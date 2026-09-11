using SparseArrays
using Statistics

function running_median(x::AbstractVector, window::Int)
    n = length(x)
    window = max(window, 1)

    # odd num window
    iseven(window) && (window += 1)

    half = window ÷ 2
    out = similar(x, Float64)

    @inbounds for i in eachindex(x)
        lo = max(firstindex(x), i - half)
        hi = min(lastindex(x), i + half)
        out[i] = median(@view x[lo:hi])
    end

    return out
end

function detection_to_liness( 
    ON, 
    OFF, 
    spec; 
    flow::Float64 = 15_000.0, 
    fhigh::Float64 = 110_000.0 
) 
    out = Vector{Vector{XY{UserUnit}}}() 
 
    fmax = maximum(spec.f) 
 
    # GUI는 y=0이 low가 아니라 화면 아래쪽이고, 
    # userunit_to_freq()에서 뒤집어서 frequency bin으로 바꿈. 
    y_low  = 1.0 - fhigh / fmax 
    y_high = 1.0 - flow  / fmax 
 
    y_low  = clamp(y_low,  0.0, 1.0) 
    y_high = clamp(y_high, 0.0, 1.0) 
 
    for (on, off) in zip(ON, OFF) 
 
        on  = clamp(on,  1, length(spec.t_sng)) 
        off = clamp(off, 1, length(spec.t_sng)) 
 
        t1 = Float64(spec.t_sng[on]) 
        t2 = Float64(spec.t_sng[off]) 
 
        push!( 
            out, 
            [ 
                XY{UserUnit}(t1, y_low), 
                XY{UserUnit}(t2, y_high) 
            ] 
        ) 
    end 
 
    sequence_check!(out) 
 
    return out 
end

function denoise_usv(sng::Matrix{Float64},f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # container for denoised spectrogram
    usv2 = deepcopy(sng)

    tind1 = findlast(x->x<15000,f)
    tind2 = findfirst(x->x>40000,f)

    noise_vals = usv2[tind1:tind2,:]
    noise_vals_vec = reshape(noise_vals,(length(noise_vals),))
    noise_sorted = sort(noise_vals_vec)

    thresh_val = noise_sorted[Int(round(length(noise_sorted)*0.95))]

    a =findall(x->x<thresh_val,usv2)

    usv2[a] .= 0

    # % Get rid of energy outside of the USV frequency range

    find1 = findall(x->x<15_000,f)
    find2 = findall(x->x>110_000,f)
    usv2[find1,:] .=0
    usv2[find2,:] .=0
   return usv2
end

# normalize a 2d array in column; sum of each column is one
function normalize_sng!(sng::Matrix{Float64},pow::BitArray{1})
    for i in 1:size(sng,2)
        sum = 0
        for j in 1:size(sng,1)
            @inbounds sum += sng[j,i]
        end
        div,id = !iszero(sum) ? (sum,false) : (1,true)
        for j in 1:size(sng,1)
            @inbounds sng[j,i] /= div
        end
        @inbounds pow[i] = id # zero
    end
    nothing
end

# ith row, jth column sdsng
function dense_sng_diff(sng::Matrix{Float64},frange::UnitRange{Int64},offset::UnitRange{Int64},i::Int64,j::Int64)
    sdsng = 0
    for n = 1:length(frange) # cacluate all in a column
        sdsng  += abs(sng[frange[n],j+1] - sng[frange[n]+offset[i],j])
    end
    return sdsng
end

function specdiscont(sng::Matrix{Float64})
    zero_ind = falses(size(sng,2))
    normalize_sng!(sng,zero_ind)

    moff = 3
    offset = -moff:moff
    frange = 1+moff:size(sng,1)-moff

    sdsng = zeros(length(offset)) # not neccessary...
    dch = zeros(Float64,size(sng,2))
    dch[zero_ind] .= 2.0
    dch[end] = 2.0

    zero_ind .= .!zero_ind # change to non-zero inds
    zero_ind[end] = 0 # make sure no calculation happens beyond the array
    zz = sparse(zero_ind)

    @inbounds for j in zz.nzind # indices for non-zero values
        for i in 1:length(offset)
             sdsng[i] = dense_sng_diff(sng,frange,offset,i,j)
        end
        dch[j] = minimum(sdsng)
    end

    return dch
end

function spectralpurity(sng::Matrix{Float64})
    specpur = zeros(size(sng,2))
    for col = 1:size(sng,2)
        temp = @view sng[:,col]
        specpur[col] = maximum(temp)/sum(temp)
    end
    nanind = findall(x->x==1,isnan.(specpur))
    specpur[nanind] .= 0
    return specpur
end

function meanfrequency(sng::Matrix{Float64},f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # sng: spectrogram
    # f: frequency range 0 to 125 khz
    mf = zeros(size(sng,2))
    @inbounds begin for col = 1:size(sng,2)
        up = 0.0
        down = 0.0
        for row = 1:size(sng,1)
             up += f[row]*sng[row,col]
             down += sng[row,col]
        end
        mf[col] = up/down
    end
    end
    nanind = findall(x->x==1,isnan.(mf))
    mf[nanind] .= 0

    return mf
end

function usv_preprocess(sng::Array{Float64,2}, dt_sng::Float64, f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})

    mf = meanfrequency(sng, f)
    sp = spectralpurity(sng)
    sd = specdiscont(sng)

    w = Int(round(0.005 / dt_sng))

    mf[2:end] .= running_median(mf, w)[2:end]
    sp[2:end] .= running_median(sp, w)[2:end]
    sd[2:end] .= running_median(sd, w)[2:end]

    return mf, sp, sd
end

function usv_ind(mf::Vector{Float64},sp::Vector{Float64},sd::Vector{Float64})
    mfthresh=20000;
    spthresh = 0.3;
    sdthresh = 0.85;

    spind = findall(x->x>spthresh,sp)
    mfind = findall(x->x>mfthresh,mf)
    sdind = findall(x->x<sdthresh,sd)

    usvind = intersect(mfind,spind,sdind)
    return usvind
end

function putative_usv(uind::Vector{Int64},t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    usv_putative = zeros(length(t))
    usv_putative[uind] .= 1
    usv_on_off = diff(usv_putative)
    usvON = findall(x->x==1,usv_on_off).-1

    if usv_putative[1] == 1
        prepend!(usvON,[1])
    end
    usvOFF = findall(x->x==-1,usv_on_off).-1
    if usv_putative[end] ==1
        append!(usvOFF,[length(t)])
    end

    return usvON, usvOFF
end

function usv_duration_check(usvON::Vector{Int64},usvOFF::Vector{Int64},dt_sng)
    keepind = Int64[];
    for z = 1:length(usvON)
        if usvOFF[z] > usvON[z] + Int(round(0.005/dt_sng))
            push!(keepind,z)
        end
    end
    usvON = usvON[keepind]
    usvOFF = usvOFF[keepind]
    return usvON, usvOFF
end

function usv_merge(usvON::Vector{Int64},usvOFF::Vector{Int64},dt_sng)
    mergeind = Int64[];
    for p = 1:length(usvOFF)-1
        if usvON[p+1] - usvOFF[p] < Int(round(0.03/dt_sng))
            push!(mergeind,p)
        end
    end
    deleteat!(usvON,mergeind.+1)
    deleteat!(usvOFF,mergeind)
    return usvON, usvOFF
end

function false_positive_usv(usvON::Vector{Int64},usvOFF::Vector{Int64},usv_putative,dt_sng)
    gapThreshold=Int(round(0.25/dt_sng))
    distanceBefore = 0
    distanceAfter = 0

    if length(usvON) > 1
        indicesToThrowOut = Int64[];
        for z = 1:length(usvON)
            if z >1
                distanceBefore=usvON[z]-usvOFF[z-1]
            else
                distanceBefore=usvON[z]
            end
            if z < length(usvON)
                distanceAfter=usvON[z+1]-usvOFF[z]
            else
                distanceAfter=length(usv_putative)-usvOFF[z]
            end
            if distanceAfter>gapThreshold;
                if distanceBefore>gapThreshold;
                    push!(indicesToThrowOut,z)
                end
            end
        end
        deleteat!(usvON,indicesToThrowOut)
        deleteat!(usvOFF,indicesToThrowOut)
    end
    return usvON, usvOFF
end

function usv_processing(uind::Vector{Int64},t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},dt_sng::Float64)
    on, off = putative_usv(uind,t)
    on,off = usv_duration_check(on,off,dt_sng)
    on,off = usv_merge(on,off,dt_sng)
    tt = Array(t)
    on,off = false_positive_usv(on,off,tt,dt_sng)

    return on,off
end

function usv_whole(sng::Matrix{Float64},dt_sng::Float64,t_sng::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # sng
    # dt_sng
    # t_sng
    # f : frequency range
    mf2,sp2,sd2 = usv_preprocess(sng,dt_sng,f)
    ui =  usv_ind(mf2,sp2,sd2)
    ON,OFF = usv_processing(ui,t_sng,dt_sng)
    return ON,OFF
end

# =====================================================================
# 2026-09 개선 추가분 (temp8.jl 전용)
#   목적:
#     (A) detection 정확도 향상
#         - amplitude/SNR gate  : 저에너지 잡음 컬럼 제거 (정밀도↑)
#         - confident false-positive : 고립된 "진짜" 콜은 유지, 고립된
#           약한 잡음만 제거 (재현율↑, 정밀도 유지)
#     (B) 세로(주파수) 방향 탐지
#         - usv_freq_bounds : 각 syllable의 실제 주파수 대역을 추정
#         - detection_to_liness_2d : box의 상/하 주파수 경계를 개별 지정
#   특징: 기존 함수는 그대로 두고 새 함수만 추가 (비파괴적).
#         Statistics/SparseArrays 외 추가 패키지 불필요.
# =====================================================================

# band index 범위를 구하는 helper (flow~fhigh Hz -> row bin 범위)
function _band_bins(f, flow::Real, fhigh::Real, nrow::Int)
    lo = findfirst(x -> x >= flow,  f)
    hi = findlast(x  -> x <= fhigh, f)
    lo === nothing && (lo = 1)
    hi === nothing && (hi = nrow)
    lo, hi = min(lo, hi), max(lo, hi)
    return lo, hi
end

# flow~fhigh 대역의 컬럼별 총 파워
function usv_band_power(sng::AbstractMatrix{<:Real}, f;
                        flow::Float64=15_000.0, fhigh::Float64=110_000.0)
    lo, hi = _band_bins(f, flow, fhigh, size(sng, 1))
    return vec(sum(@view(sng[lo:hi, :]), dims=1))
end

# robust noise floor:
#   median(nonzero) 는 USV 가 sparse 하고 denoise 가 배경을 지운 데이터에서
#   "신호 레벨"이 되어버려(잡음 아님) gate/판정 임계가 폭등 -> 전부 제거 -> 0 detection.
#   대신 nonzero 컬럼 파워의 낮은 분위수(q)를 잡음 바닥으로 사용.
function _noise_floor(colpow::AbstractVector{<:Real}; q::Float64=0.25)
    nz = colpow[colpow .> 0]
    isempty(nz) && return 0.0
    return quantile(nz, q)
end

# (A) amplitude / SNR gate : 컬럼 파워가 noise floor 의 snr_k배 초과인 컬럼만 True.
#   ※ 이 gate 는 조용한 onset/offset 을 잘라 syllable 을 truncate 하므로
#     usv_whole_2d 에서 기본 off. 필요할 때만 opt-in.
function usv_amplitude_gate(sng::AbstractMatrix{<:Real}, f;
                            flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                            snr_k::Float64=3.0, floor_q::Float64=0.25)
    p = usv_band_power(sng, f; flow=flow, fhigh=fhigh)
    floor_val = _noise_floor(p; q=floor_q)
    return p .> (snr_k * floor_val)
end

# (A) confident false-positive:
#   앞뒤로 모두 gap 이상 떨어진 "고립" 구간을, 그 자체가 충분히 강할 때만 유지.
#   기존 false_positive_usv 는 고립이면 무조건 제거 -> sparse 녹음(모든 콜이 고립)에서
#   전부 삭제되어 0 detection 이 됨. 여기서는:
#     - robust noise floor 사용 (median 대신 낮은 분위수)
#     - keep 판정을 느슨하게 (keep_snr/keep_pur)
#     - fail-safe: 전부 지워질 상황이면 아무것도 지우지 않음 (0 방지)
function false_positive_usv_confident(usvON::Vector{Int64}, usvOFF::Vector{Int64},
                                      n::Int, dt_sng::Float64,
                                      sng::AbstractMatrix{<:Real}, f;
                                      gapThreshold_s::Float64=0.25,
                                      keep_snr::Float64=1.5, keep_pur::Float64=0.25,
                                      floor_q::Float64=0.25,
                                      flow::Float64=15_000.0, fhigh::Float64=110_000.0)
    length(usvON) <= 1 && return usvON, usvOFF
    gap = Int(round(gapThreshold_s / dt_sng))

    p = usv_band_power(sng, f; flow=flow, fhigh=fhigh)
    floor_val = _noise_floor(p; q=floor_q)

    lo, hi = _band_bins(f, flow, fhigh, size(sng, 1))
    ncol = size(sng, 2)

    throwout = Int64[]
    for z in 1:length(usvON)
        db = z > 1 ? usvON[z] - usvOFF[z-1] : usvON[z]
        da = z < length(usvON) ? usvON[z+1] - usvOFF[z] : n - usvOFF[z]
        isolated = (da > gap) && (db > gap)
        isolated || continue

        a = clamp(usvON[z],  1, ncol)
        b = clamp(usvOFF[z], 1, ncol)
        seg  = @view sng[lo:hi, a:b]
        colp = vec(sum(seg, dims=1))
        amp  = any(colp .> 0) ? median(colp[colp .> 0]) : 0.0
        pur  = mean(vec(maximum(seg, dims=1)) ./ max.(colp, eps()))

        strong = (pur >= keep_pur) && (amp >= keep_snr * floor_val)
        strong || push!(throwout, z)
    end

    # fail-safe: 모두 제거될 상황이면 원본 유지 (0 detection 방지)
    if length(throwout) == length(usvON)
        return usvON, usvOFF
    end

    keep = setdiff(1:length(usvON), throwout)
    return usvON[keep], usvOFF[keep]
end

# (B) 세로 방향 탐지: 한 syllable 구간(on:off)의 실제 주파수 경계 (Hz) 추정.
#   각 컬럼의 ridge(argmax, 잡음에 강건)를 모아 분위수로 대역을 잡음.
#   컬럼 선택 방식(select):
#     :local  (기본) 각 컬럼의 peak 가 "그 컬럼 자신의 잡음"의 local_k 배를 넘으면 포함.
#             → 진폭이 변하는 스윕에서 조용한 쪽 주파수까지 포착 (편향 없음).
#     :global 세그먼트 전체에서 상위(1-col_energy_q) 로 강한 컬럼만 포함 (구버전).
#             → 조용한 부분을 버려 강한 쪽 주파수로 편향될 수 있음.
#   평탄음은 좁게, 큰 스윕은 넓게 자동 적응.
function usv_freq_bounds(sng::AbstractMatrix{<:Real}, f, on::Integer, off::Integer;
                         flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                         select::Symbol=:local, local_k::Float64=1.5,
                         col_energy_q::Float64=0.60,
                         ridge_q_lo::Float64=0.05, ridge_q_hi::Float64=0.99,
                         margin_bins::Int=4, min_height_bins::Int=6)
    ncol = size(sng, 2)
    on  = clamp(Int(on),  1, ncol)
    off = clamp(Int(off), 1, ncol)
    on, off = min(on, off), max(on, off)

    lo, hi = _band_bins(f, flow, fhigh, size(sng, 1))

    seg    = @view sng[lo:hi, on:off]
    coltot = vec(sum(seg, dims=1))
    pos    = coltot[coltot .> 0]
    isempty(pos) && return flow, fhigh

    # 전역 방식용 임계 (select==:global 일 때만 사용)
    gthr = quantile(pos, col_energy_q)

    peaks = Int[]
    @inbounds for c in 1:size(seg, 2)
        col = @view seg[:, c]
        pk  = maximum(col)
        pk <= 0 && continue

        include_col = if select === :local
            # 컬럼 자체의 잡음(nonzero 중앙값) 대비 local SNR
            locnoise = _median_nonzero(col)
            pk > local_k * locnoise
        else
            coltot[c] >= gthr
        end
        include_col || continue

        r = argmax(col)               # band 내 1-based ridge
        push!(peaks, r + lo - 1)      # 전체 스펙트럼 bin 으로 환산
    end
    isempty(peaks) && return flow, fhigh

    b_lo = floor(Int, quantile(peaks, ridge_q_lo)) - margin_bins
    b_hi = ceil(Int,  quantile(peaks, ridge_q_hi)) + margin_bins

    # 최소 높이 보장 (feature 계산이 안정적으로 되도록)
    if b_hi - b_lo < min_height_bins
        c  = (b_hi + b_lo) ÷ 2
        b_lo = c - min_height_bins ÷ 2
        b_hi = c + min_height_bins ÷ 2
    end
    b_lo = clamp(b_lo, lo, hi)
    b_hi = clamp(b_hi, lo, hi)

    return Float64(f[b_lo]), Float64(f[b_hi])
end

# 한 컬럼(벡터)의 nonzero 중앙값 (per-column 잡음 추정용)
function _median_nonzero(col::AbstractVector{<:Real})
    nz = [v_ for v_ in col if v_ > 0]
    isempty(nz) && return 0.0
    return median(nz)
end

# (B) box 변환: 시간 ON/OFF + 개별 주파수 경계로 GUI box 생성
#   좌표계: y = 1 - f/fmax  (y 작을수록 고주파, detection_to_liness 와 동일 규약)
function detection_to_liness_2d(ON, OFF, spec, sng::AbstractMatrix{<:Real};
                                flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                                kwargs...)
    out  = Vector{Vector{XY{UserUnit}}}()
    fmax = maximum(spec.f)

    for (on, off) in zip(ON, OFF)
        on_c  = clamp(on,  1, length(spec.t_sng))
        off_c = clamp(off, 1, length(spec.t_sng))
        t1 = Float64(spec.t_sng[on_c])
        t2 = Float64(spec.t_sng[off_c])

        f_lo, f_hi = usv_freq_bounds(sng, spec.f, on, off;
                                     flow=flow, fhigh=fhigh, kwargs...)

        y_hi = clamp(1.0 - f_hi / fmax, 0.0, 1.0)   # 작은 y  = 고주파
        y_lo = clamp(1.0 - f_lo / fmax, 0.0, 1.0)   # 큰 y   = 저주파

        push!(out, [XY{UserUnit}(t1, y_hi), XY{UserUnit}(t2, y_lo)])
    end

    sequence_check!(out)
    return out
end

# (A) 히스테리시스 경계 확장 (dual-threshold):
#   엄격한 core detection 이 잡은 각 구간을, 신호가 noise floor 근처로
#   떨어질 때까지 좌/우로 확장한다. 강한 core 에서만 자라기 때문에
#   순수 잡음 구간은 절대 생성되지 않는다(정밀도 유지).
#   → "가장 강한 부분만 잡히고 조용한 onset/offset 이 잘리는" 문제 해결.
#   sp 는 usv_preprocess 로 smoothing 된 것을 사용.
function extend_boundaries!(ON::Vector{Int64}, OFF::Vector{Int64},
                            sp::Vector{Float64}, colpow::Vector{Float64},
                            floor_val::Float64;
                            snr_ext::Float64=1.2, sp_lo::Float64=0.15)
    ncol = length(colpow)
    n = length(ON)
    thr = snr_ext * floor_val
    @inbounds for k in 1:n
        # 왼쪽으로 확장 (직전 구간을 넘지 않도록 bound)
        left_bound = k > 1 ? OFF[k-1] : 0
        i = ON[k]
        while i > 1 && (i - 1) > left_bound &&
              colpow[i-1] > thr && sp[i-1] > sp_lo
            i -= 1
        end
        ON[k] = i

        # 오른쪽으로 확장 (다음 구간 core 를 넘지 않도록 bound)
        right_bound = k < n ? ON[k+1] : ncol + 1
        j = OFF[k]
        while j < ncol && (j + 1) < right_bound &&
              colpow[j+1] > thr && sp[j+1] > sp_lo
            j += 1
        end
        OFF[k] = j
    end
    return ON, OFF
end

# 개선된 전체 detection (시간축 정확도 + gate + 경계 확장).
# 세로 경계는 상위에서 detection_to_liness_2d 로 붙임.
# use_*/extend 플래그로 기존 동작과 비교 가능.
#   snr_k   : core seed 용 (엄격, 정밀도 담당)
#   snr_ext : 경계 확장용 (느슨, 연속 syllable 전체 길이 회복) — snr_k 보다 낮게
function usv_whole_2d(sng::Matrix{Float64}, dt_sng::Float64,
                      t_sng::StepRangeLen, f::StepRangeLen;
                      use_amp_gate::Bool=false,   # ← 기본 off: syllable truncate/zero 유발
                      snr_k::Float64=0.8,
                      fp_mode::Symbol=:none, # :confident | :original | :none
                      extend::Bool=true,
                      snr_ext::Float64=1.3,
                      sp_lo::Float64=0.15,
                      floor_q::Float64=0.15,
                      mfthresh::Float64=20_000.0,
                      spthresh::Float64=0.1,
                      sdthresh::Float64=0.99,
                      verbose::Bool=true)
    mf2, sp2, sd2 = usv_preprocess(sng, dt_sng, f)

    spind = findall(x -> x > spthresh,  sp2)
    mfind = findall(x -> x > mfthresh,  mf2)
    sdind = findall(x -> x < sdthresh,  sd2)
    uind  = intersect(mfind, spind, sdind)
    verbose && println("  [stage] core columns (sp&mf&sd): ", length(uind))

    if use_amp_gate
        gate = usv_amplitude_gate(sng, f; snr_k=snr_k, floor_q=floor_q)
        uind = filter(i -> gate[i], uind)
        verbose && println("  [stage] after amp-gate:        ", length(uind))
    end

    on, off = putative_usv(uind, t_sng)
    on, off = usv_duration_check(on, off, dt_sng)
    on, off = usv_merge(on, off, dt_sng)
    verbose && println("  [stage] segments after merge:  ", length(on))

    tt = Array(t_sng)
    if fp_mode === :confident
        on, off = false_positive_usv_confident(on, off, length(tt), dt_sng, sng, f;
                                               keep_snr=1.5, keep_pur=0.25, floor_q=floor_q)
    elseif fp_mode === :original
        on, off = false_positive_usv(on, off, tt, dt_sng)
    elseif fp_mode === :none
        # false-positive 필터 건너뜀: 고립 검출을 지우지 않음 (재현율↑, 정밀도↓)
    else
        error("fp_mode must be :confident, :original, or :none (got $fp_mode)")
    end
    verbose && println("  [stage] segments after FP ($(fp_mode)): ", length(on))

    # 경계 확장: core 를 실제 syllable 끝까지 늘림 (구간을 늘리기만 함 → 0 유발 불가)
    if extend && !isempty(on)
        colpow = usv_band_power(sng, f)
        floor_val = _noise_floor(colpow; q=floor_q)
        extend_boundaries!(on, off, sp2, colpow, floor_val; snr_ext=snr_ext, sp_lo=sp_lo)
        on, off = usv_merge(on, off, dt_sng)  # 확장 후 붙은 구간 병합
    end
    verbose && println("  [stage] final segments:        ", length(on))

    return on, off
end