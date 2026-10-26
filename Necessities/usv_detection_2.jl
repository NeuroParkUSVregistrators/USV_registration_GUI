
using SparseArrays
using Statistics

function running_median(x::AbstractVector, window::Int)
    window = max(window, 1)
    iseven(window) && (window += 1)   # odd num window

    half = window ÷ 2
    out = similar(x, Float64)

    @inbounds for i in eachindex(x)
        lo = max(firstindex(x), i - half)
        hi = min(lastindex(x), i + half)
        out[i] = median(@view x[lo:hi])
    end
    return out
end

function denoise_usv(sng::Matrix{Float64}, f::AbstractVector{<:Real})
    usv2 = deepcopy(sng)

    tind1 = findlast(x -> x < 15000, f)
    tind2 = findfirst(x -> x > 40000, f)

    noise_vals = usv2[tind1:tind2, :]
    noise_sorted = sort(vec(noise_vals))
    thresh_val = noise_sorted[Int(round(length(noise_sorted) * 0.95))]

    usv2[findall(x -> x < thresh_val, usv2)] .= 0

    # USV 주파수 대역 밖 에너지 제거
    usv2[findall(x -> x < 15_000, f), :] .= 0
    usv2[findall(x -> x > 110_000, f), :] .= 0
    return usv2
end

# normalize a 2d array in column; sum of each column is one  (※ in-place)
function normalize_sng!(sng::Matrix{Float64}, pow::BitArray{1})
    for i in 1:size(sng, 2)
        s = 0.0
        for j in 1:size(sng, 1)
            @inbounds s += sng[j, i]
        end
        div, id = !iszero(s) ? (s, false) : (1.0, true)
        for j in 1:size(sng, 1)
            @inbounds sng[j, i] /= div
        end
        @inbounds pow[i] = id # zero column
    end
    nothing
end

# ith offset, jth column sdsng
function dense_sng_diff(sng::Matrix{Float64}, frange::UnitRange{Int64},
                        offset::UnitRange{Int64}, i::Int64, j::Int64)
    sdsng = 0.0
    for n = 1:length(frange)
        sdsng += abs(sng[frange[n], j+1] - sng[frange[n]+offset[i], j])
    end
    return sdsng
end

# ※ 인자를 in-place 로 정규화함. 원본을 보존하려면 copy 를 넘길 것 (usv_preprocess 참고).
function specdiscont(sng::Matrix{Float64})
    zero_ind = falses(size(sng, 2))
    normalize_sng!(sng, zero_ind)

    moff = 3
    offset = -moff:moff
    frange = 1+moff:size(sng, 1)-moff

    sdsng = zeros(length(offset))
    dch = zeros(Float64, size(sng, 2))
    dch[zero_ind] .= 2.0
    dch[end] = 2.0

    zero_ind .= .!zero_ind   # non-zero inds
    zero_ind[end] = 0        # 배열 밖 계산 방지
    zz = sparse(zero_ind)

    @inbounds for j in zz.nzind
        for i in 1:length(offset)
            sdsng[i] = dense_sng_diff(sng, frange, offset, i, j)
        end
        dch[j] = minimum(sdsng)
    end
    return dch
end

function spectralpurity(sng::AbstractMatrix{<:Real})
    specpur = zeros(size(sng, 2))
    for col = 1:size(sng, 2)
        temp = @view sng[:, col]
        specpur[col] = maximum(temp) / sum(temp)
    end
    specpur[isnan.(specpur)] .= 0
    return specpur
end

function meanfrequency(sng::AbstractMatrix{<:Real}, f::AbstractVector{<:Real})
    mf = zeros(size(sng, 2))
    @inbounds for col = 1:size(sng, 2)
        up = 0.0
        down = 0.0
        for row = 1:size(sng, 1)
            up += f[row] * sng[row, col]
            down += sng[row, col]
        end
        mf[col] = up / down
    end
    mf[isnan.(mf)] .= 0
    return mf
end

function usv_preprocess(sng::Matrix{Float64}, dt_sng::Float64, f::AbstractVector{<:Real})
    mf = meanfrequency(sng, f)
    sp = spectralpurity(sng)
    sd = specdiscont(copy(sng))   # ← copy: 호출자의 sng 가 정규화되지 않도록

    w = Int(round(0.005 / dt_sng))
    mf[2:end] .= running_median(mf, w)[2:end]
    sp[2:end] .= running_median(sp, w)[2:end]
    sd[2:end] .= running_median(sd, w)[2:end]

    return mf, sp, sd
end

# band index 범위를 구하는 helper (flow~fhigh Hz -> row bin 범위)
function _band_bins(f, flow::Real, fhigh::Real, nrow::Int)
    lo = findfirst(x -> x >= flow, f)
    hi = findlast(x -> x <= fhigh, f)
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

# robust noise floor: nonzero 컬럼 파워의 낮은 분위수(q)
#   (median(nonzero) 는 denoise 된 sparse 데이터에서 신호 레벨이 되어버림)
function _noise_floor(colpow::AbstractVector{<:Real}; q::Float64=0.25)
    nz = colpow[colpow .> 0]
    isempty(nz) && return 0.0
    return quantile(nz, q)
end

# 한 컬럼(벡터)의 nonzero 중앙값 (per-column 잡음 추정용)
function _median_nonzero(col::AbstractVector{<:Real})
    nz = [v_ for v_ in col if v_ > 0]
    isempty(nz) && return 0.0
    return median(nz)
end

# ON/OFF 를 [1, ncol] 로 clamp 하고 on <= off 로 정렬
function _sanitize_onoff(ON::AbstractVector{<:Integer}, OFF::AbstractVector{<:Integer}, ncol::Int)
    length(ON) == length(OFF) || throw(DimensionMismatch("ON and OFF must have equal length"))
    on  = clamp.(Int.(ON),  1, ncol)
    off = clamp.(Int.(OFF), 1, ncol)
    return min.(on, off), max.(on, off)
end

# =====================================================================
# 2. 시간축 detection
# =====================================================================

function usv_ind(mf::Vector{Float64}, sp::Vector{Float64}, sd::Vector{Float64};
                 mfthresh::Real=20_000, spthresh::Real=0.3, sdthresh::Real=0.85)
    spind = findall(x -> x > spthresh, sp)
    mfind = findall(x -> x > mfthresh, mf)
    sdind = findall(x -> x < sdthresh, sd)
    return intersect(mfind, spind, sdind)
end

# ※ 원본 MATLAB 포팅 규약 유지: ON/OFF 가 실제 전이 지점보다 1~2 bin 앞섬, ON=0 가능.
#   하위 함수들은 _sanitize_onoff / clamp 로 처리함.
function putative_usv(uind::AbstractVector{<:Integer}, t::AbstractVector)
    usv_putative = zeros(length(t))
    usv_putative[uind] .= 1
    usv_on_off = diff(usv_putative)
    usvON = findall(x -> x == 1, usv_on_off) .- 1
    usv_putative[1] == 1 && prepend!(usvON, [1])
    usvOFF = findall(x -> x == -1, usv_on_off) .- 1
    usv_putative[end] == 1 && append!(usvOFF, [length(t)])
    return usvON, usvOFF
end

function usv_duration_check(usvON::Vector{Int64}, usvOFF::Vector{Int64}, dt_sng)
    minlen = Int(round(0.005 / dt_sng))
    keepind = [z for z in 1:length(usvON) if usvOFF[z] > usvON[z] + minlen]
    return usvON[keepind], usvOFF[keepind]
end

function usv_merge(usvON::Vector{Int64}, usvOFF::Vector{Int64}, dt_sng)
    gap = Int(round(0.03 / dt_sng))
    mergeind = [p for p in 1:length(usvOFF)-1 if usvON[p+1] - usvOFF[p] < gap]
    deleteat!(usvON, mergeind .+ 1)
    deleteat!(usvOFF, mergeind)
    return usvON, usvOFF
end

function false_positive_usv(usvON::Vector{Int64}, usvOFF::Vector{Int64}, usv_putative, dt_sng)
    gapThreshold = Int(round(0.25 / dt_sng))
    if length(usvON) > 1
        indicesToThrowOut = Int64[]
        for z = 1:length(usvON)
            distanceBefore = z > 1 ? usvON[z] - usvOFF[z-1] : usvON[z]
            distanceAfter  = z < length(usvON) ? usvON[z+1] - usvOFF[z] :
                                                 length(usv_putative) - usvOFF[z]
            if distanceAfter > gapThreshold && distanceBefore > gapThreshold
                push!(indicesToThrowOut, z)
            end
        end
        deleteat!(usvON, indicesToThrowOut)
        deleteat!(usvOFF, indicesToThrowOut)
    end
    return usvON, usvOFF
end

function usv_processing(uind::AbstractVector{<:Integer}, t::AbstractVector, dt_sng::Float64)
    on, off = putative_usv(uind, t)
    on, off = usv_duration_check(on, off, dt_sng)
    on, off = usv_merge(on, off, dt_sng)
    on, off = false_positive_usv(on, off, t, dt_sng)
    return on, off
end

# 원본 detection (비교용으로 유지)
function usv_whole(sng::Matrix{Float64}, dt_sng::Float64,
                   t_sng::AbstractVector, f::AbstractVector{<:Real})
    mf2, sp2, sd2 = usv_preprocess(sng, dt_sng, f)
    ui = usv_ind(mf2, sp2, sd2)
    return usv_processing(ui, t_sng, dt_sng)
end

# (A) amplitude / SNR gate : 컬럼 파워가 noise floor 의 snr_k배 초과인 컬럼만 true.
#   ※ 조용한 onset/offset 을 잘라 syllable 을 truncate 할 수 있어 usv_whole_2d 에서 기본 off.
function usv_amplitude_gate(sng::AbstractMatrix{<:Real}, f;
                            flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                            snr_k::Float64=3.0, floor_q::Float64=0.25)
    p = usv_band_power(sng, f; flow=flow, fhigh=fhigh)
    floor_val = _noise_floor(p; q=floor_q)
    return p .> (snr_k * floor_val)
end

# (A) confident false-positive:
#   고립 구간을 그 자체가 충분히 강할 때만 유지. 전부 지워질 상황이면 아무것도 안 지움.
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
        ((da > gap) && (db > gap)) || continue   # 고립된 구간만 판정

        a = clamp(usvON[z], 1, ncol)
        b = clamp(usvOFF[z], 1, ncol)
        a, b = min(a, b), max(a, b)
        seg  = @view sng[lo:hi, a:b]
        colp = vec(sum(seg, dims=1))
        amp  = any(colp .> 0) ? median(colp[colp .> 0]) : 0.0
        pur  = mean(vec(maximum(seg, dims=1)) ./ max.(colp, eps()))

        strong = (pur >= keep_pur) && (amp >= keep_snr * floor_val)
        strong || push!(throwout, z)
    end

    length(throwout) == length(usvON) && return usvON, usvOFF   # fail-safe

    keep = setdiff(1:length(usvON), throwout)
    return usvON[keep], usvOFF[keep]
end

# (A) 히스테리시스 경계 확장 (dual-threshold):
#   core 구간을 신호가 noise floor 근처로 떨어질 때까지 좌/우로 확장.
#   core 에서만 자라므로 순수 잡음 구간은 생성되지 않음.
function extend_boundaries!(ON::Vector{Int64}, OFF::Vector{Int64},
                            sp::Vector{Float64}, colpow::Vector{Float64},
                            floor_val::Float64;
                            snr_ext::Float64=1.2, sp_lo::Float64=0.15)
    ncol = length(colpow)
    n = length(ON)
    thr = snr_ext * floor_val
    @inbounds for k in 1:n
        # 시작점을 유효 범위로 (putative_usv 는 ON=0 을 낼 수 있음)
        ON[k]  = clamp(ON[k],  1, ncol)
        OFF[k] = clamp(OFF[k], 1, ncol)

        left_bound = k > 1 ? OFF[k-1] : 0
        i = ON[k]
        while i > 1 && (i - 1) > left_bound && colpow[i-1] > thr && sp[i-1] > sp_lo
            i -= 1
        end
        ON[k] = i

        right_bound = k < n ? ON[k+1] : ncol + 1
        j = OFF[k]
        while j < ncol && (j + 1) < right_bound && colpow[j+1] > thr && sp[j+1] > sp_lo
            j += 1
        end
        OFF[k] = j
    end
    return ON, OFF
end

# 개선된 전체 detection (시간축). 세로 경계는 usv_segment_freq_bounds 로 붙임.
#   snr_k   : amp-gate core seed 용 (엄격)
#   snr_ext : 경계 확장용 (느슨)
function usv_whole_2d(sng::Matrix{Float64}, dt_sng::Float64,
                      t_sng::AbstractVector, f::AbstractVector{<:Real};
                      use_amp_gate::Bool=false,
                      snr_k::Float64=0.8,
                      fp_mode::Symbol=:none,      # :confident | :original | :none
                      extend::Bool=true,
                      snr_ext::Float64=1.3,
                      sp_lo::Float64=0.15,
                      floor_q::Float64=0.15,
                      mfthresh::Float64=20_000.0,
                      spthresh::Float64=0.1,
                      sdthresh::Float64=0.99,
                      verbose::Bool=true)
    fp_mode in (:confident, :original, :none) ||
        error("fp_mode must be :confident, :original, or :none (got $fp_mode)")

    mf2, sp2, sd2 = usv_preprocess(sng, dt_sng, f)
    uind = usv_ind(mf2, sp2, sd2; mfthresh=mfthresh, spthresh=spthresh, sdthresh=sdthresh)
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

    if fp_mode === :confident
        on, off = false_positive_usv_confident(on, off, length(t_sng), dt_sng, sng, f;
                                               keep_snr=1.5, keep_pur=0.25, floor_q=floor_q)
    elseif fp_mode === :original
        on, off = false_positive_usv(on, off, t_sng, dt_sng)
    end
    verbose && println("  [stage] segments after FP ($(fp_mode)): ", length(on))

    if extend && !isempty(on)
        colpow = usv_band_power(sng, f)
        floor_val = _noise_floor(colpow; q=floor_q)
        extend_boundaries!(on, off, sp2, colpow, floor_val; snr_ext=snr_ext, sp_lo=sp_lo)
        on, off = usv_merge(on, off, dt_sng)
    end
    verbose && println("  [stage] final segments:        ", length(on))

    return on, off
end

# =====================================================================
# 3. 주파수(세로) 경계 탐지
# =====================================================================

# --- 3-1. :ridge 방식 ------------------------------------------------
# 한 syllable(on:off)의 주파수 경계(Hz). 각 컬럼 ridge(argmax)를 모아 분위수로 대역 결정.
#   select = :local  (기본) 컬럼 peak 가 그 컬럼 자체 잡음의 local_k 배 초과면 포함
#            :global 세그먼트 전체에서 상위 (1-col_energy_q) 컬럼만 포함
function usv_freq_bounds(sng::AbstractMatrix{<:Real}, f, on::Integer, off::Integer;
                         flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                         select::Symbol=:local, local_k::Float64=1.5,
                         col_energy_q::Float64=0.60,
                         ridge_q_lo::Float64=0.05, ridge_q_hi::Float64=0.99,
                         margin_bins::Int=4, min_height_bins::Int=6)
    ncol = size(sng, 2)
    on  = clamp(Int(on), 1, ncol)
    off = clamp(Int(off), 1, ncol)
    on, off = min(on, off), max(on, off)

    lo, hi = _band_bins(f, flow, fhigh, size(sng, 1))

    seg    = @view sng[lo:hi, on:off]
    coltot = vec(sum(seg, dims=1))
    pos    = coltot[coltot .> 0]
    isempty(pos) && return flow, fhigh

    gthr = quantile(pos, col_energy_q)   # select==:global 용

    peaks = Int[]
    @inbounds for c in 1:size(seg, 2)
        col = @view seg[:, c]
        pk  = maximum(col)
        pk <= 0 && continue
        include_col = select === :local ? pk > local_k * _median_nonzero(col) :
                                          coltot[c] >= gthr
        include_col || continue
        push!(peaks, argmax(col) + lo - 1)   # 전체 스펙트럼 bin 으로 환산
    end
    isempty(peaks) && return flow, fhigh

    b_lo = floor(Int, quantile(peaks, ridge_q_lo)) - margin_bins
    b_hi = ceil(Int, quantile(peaks, ridge_q_hi)) + margin_bins
    b_lo, b_hi = _ensure_min_height(b_lo, b_hi, min_height_bins, lo, hi)

    return Float64(f[b_lo]), Float64(f[b_hi])
end

# 최소 높이 보장 후 band 로 clamp (두 방식 공통)
function _ensure_min_height(b_lo::Int, b_hi::Int, min_height_bins::Int, lo::Int, hi::Int)
    if b_hi - b_lo < min_height_bins
        c = (b_hi + b_lo) ÷ 2
        b_lo = c - min_height_bins ÷ 2
        b_hi = c + min_height_bins ÷ 2
    end
    return clamp(b_lo, lo, hi), clamp(b_hi, lo, hi)
end

# --- 3-2. :threshold 방식 (구 TimeBinDetector) ------------------------
# 각 time bin 마다 threshold = max(absolute_threshold,
#                                  median + noise_k·1.4826·MAD,
#                                  peak_fraction·peak)
# 을 넘는 최저/최고 bin 을 찾고, 세그먼트 전체의 min/max 를 경계로 사용.
# 반환: 세그먼트별 NamedTuple
#   (on, off, time_bins, bin_min, bin_max, freq_min, freq_max,   ← time bin 별 contour
#    min_bin, max_bin, min_freq, max_freq, valid)                ← 세그먼트 전체 경계
function detect_frequency_bounds(
    sng::AbstractMatrix{<:Real}, f::AbstractVector{<:Real},
    ON::AbstractVector{<:Integer}, OFF::AbstractVector{<:Integer};
    flow::Real=15_000, fhigh::Real=110_000, noise_k::Real=3.0,
    peak_fraction::Real=0.05, absolute_threshold::Real=0.0,
    margin_bins::Integer=0,
)
    Base.require_one_based_indexing(sng, f, ON, OFF)
    length(f) == size(sng, 1) || throw(DimensionMismatch("f must match spectrogram rows"))
    length(ON) == length(OFF) || throw(DimensionMismatch("ON and OFF must have equal length"))
    isempty(f) && throw(ArgumentError("frequency axis must not be empty"))
    all(isfinite, f) && all(diff(f) .> 0) ||
        throw(ArgumentError("frequencies must be finite and strictly increasing"))
    isfinite(flow) && isfinite(fhigh) && 0 <= flow <= fhigh ||
        throw(ArgumentError("require 0 <= flow <= fhigh, both finite"))
    isfinite(noise_k) && noise_k >= 0 || throw(ArgumentError("noise_k must be finite and nonnegative"))
    isfinite(peak_fraction) && 0 <= peak_fraction < 1 ||
        throw(ArgumentError("peak_fraction must be in [0, 1)"))
    isfinite(absolute_threshold) && absolute_threshold >= 0 ||
        throw(ArgumentError("absolute_threshold must be finite and nonnegative"))
    margin_bins >= 0 || throw(ArgumentError("margin_bins must be nonnegative"))
    lo, hi = searchsortedfirst(f, flow), searchsortedlast(f, fhigh)
    lo <= hi || throw(ArgumentError("no frequency bins in the requested band"))
    for (on, off) in zip(ON, OFF)
        1 <= on <= off <= size(sng, 2) ||
            throw(ArgumentError("require 1 <= ON <= OFF <= number of time bins"))
    end

    results = NamedTuple[]
    for (on, off) in zip(ON, OFF)
        time_bins = Int(on):Int(off)
        bin_min = Vector{Union{Missing,Int}}(missing, length(time_bins))
        bin_max = copy(bin_min)
        freq_min = Vector{Union{Missing,Float64}}(missing, length(time_bins))
        freq_max = copy(freq_min)
        overall_lo, overall_hi = hi + 1, lo - 1
        for (j, t) in enumerate(time_bins)
            col = @view sng[lo:hi, t]
            values = Float64[v for v in col if isfinite(v)]
            isempty(values) && continue
            minimum(values) >= 0 || throw(ArgumentError("sng must contain nonnegative linear values"))
            peak = maximum(values)
            peak > 0 || continue
            baseline = median(values)
            mad = median(abs.(values .- baseline))
            # noise_k=0 explicitly disables the noise-floor criterion.
            noise_threshold = iszero(noise_k) ? 0.0 : baseline + noise_k * 1.4826 * mad
            threshold = max(absolute_threshold, noise_threshold, peak_fraction * peak)
            first_bin, last_bin = hi + 1, lo - 1
            for r in lo:hi
                v = sng[r, t]
                if isfinite(v) && v > threshold
                    first_bin = min(first_bin, r)
                    last_bin = max(last_bin, r)
                end
            end
            first_bin <= last_bin || continue
            bin_min[j], bin_max[j] = first_bin, last_bin
            freq_min[j], freq_max[j] = Float64(f[first_bin]), Float64(f[last_bin])
            overall_lo = min(overall_lo, first_bin)
            overall_hi = max(overall_hi, last_bin)
        end
        valid = overall_lo <= overall_hi
        min_bin = valid ? overall_lo - min(margin_bins, overall_lo - lo) : missing
        max_bin = valid ? overall_hi + min(margin_bins, hi - overall_hi) : missing
        push!(results, (; on=Int(on), off=Int(off), time_bins, bin_min, bin_max,
                          freq_min, freq_max, min_bin, max_bin,
                          min_freq=valid ? Float64(f[min_bin]) : missing,
                          max_freq=valid ? Float64(f[max_bin]) : missing, valid))
    end
    return results
end

# detect_frequency_bounds 를 detection 출력(ON/OFF)에 바로 쓰는 adapter.
#   - ON/OFF 를 sanitize 해서 넘김
#   - invalid 세그먼트는 (flow, fhigh) 로 fallback (ridge 방식과 동일 규약)
#   - min_height_bins 보장 (ridge 방식과 동일)
#   - min_col_peak_frac: 세그먼트 최대 peak 대비 이 비율 미만인 time bin 은 경계 계산에서 제외.
#       detect_frequency_bounds 는 모든 time bin 의 min/max 를 합치므로, 경계 확장으로
#       섞여 들어온 약한 잡음 컬럼 하나가 box 를 전체 대역으로 벌릴 수 있음 → 이를 방지.
#       0.0 이면 원래 detect_frequency_bounds 결과와 동일.
# 반환: Vector{Tuple{Float64,Float64}}  (f_lo, f_hi) per segment
function usv_freq_bounds_threshold(sng::AbstractMatrix{<:Real}, f::AbstractVector{<:Real},
                                   ON::AbstractVector{<:Integer}, OFF::AbstractVector{<:Integer};
                                   flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                                   margin_bins::Int=4, min_height_bins::Int=6,
                                   noise_k::Real=3.0, peak_fraction::Real=0.05,
                                   absolute_threshold::Real=0.0,
                                   min_col_peak_frac::Real=0.2)
    isempty(ON) && return Tuple{Float64,Float64}[]
    on, off = _sanitize_onoff(ON, OFF, size(sng, 2))
    res = detect_frequency_bounds(sng, f, on, off;
                                  flow=flow, fhigh=fhigh, noise_k=noise_k,
                                  peak_fraction=peak_fraction,
                                  absolute_threshold=absolute_threshold,
                                  margin_bins=0)          # margin 은 아래에서 직접 적용
    lo, hi = searchsortedfirst(f, flow), searchsortedlast(f, fhigh)
    out = Vector{Tuple{Float64,Float64}}(undef, length(res))
    for (k, r) in enumerate(res)
        if !r.valid
            out[k] = (flow, fhigh)
            continue
        end
        # time bin 별 band 내 peak
        colpeak = [maximum(x -> isfinite(x) ? Float64(x) : 0.0, @view sng[lo:hi, t])
                   for t in r.time_bins]
        cut = min_col_peak_frac * maximum(colpeak)
        b_lo, b_hi = hi + 1, lo - 1
        for j in eachindex(r.time_bins)
            (ismissing(r.bin_min[j]) || colpeak[j] < cut) && continue
            b_lo = min(b_lo, r.bin_min[j])
            b_hi = max(b_hi, r.bin_max[j])
        end
        if b_lo > b_hi                                    # 강한 컬럼이 임계를 못 넘은 경우 fallback
            b_lo, b_hi = r.min_bin, r.max_bin
        end
        b_lo, b_hi = _ensure_min_height(b_lo - margin_bins, b_hi + margin_bins,
                                        min_height_bins, lo, hi)
        out[k] = (Float64(f[b_lo]), Float64(f[b_hi]))
    end
    return out
end

# --- 3-3. 공통 인터페이스 ---------------------------------------------
# fmethod = :ridge | :threshold. kwargs 는 선택한 방식의 함수로 그대로 전달.
# 반환: Vector{Tuple{Float64,Float64}}  (f_lo, f_hi) per segment
function usv_segment_freq_bounds(sng::AbstractMatrix{<:Real}, f::AbstractVector{<:Real},
                                 ON::AbstractVector{<:Integer}, OFF::AbstractVector{<:Integer};
                                 fmethod::Symbol=:ridge,
                                 flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                                 kwargs...)
    if fmethod === :ridge
        on, off = _sanitize_onoff(ON, OFF, size(sng, 2))
        return Tuple{Float64,Float64}[usv_freq_bounds(sng, f, a, b; flow=flow, fhigh=fhigh, kwargs...)
                                      for (a, b) in zip(on, off)]
    elseif fmethod === :threshold
        return usv_freq_bounds_threshold(sng, f, ON, OFF; flow=flow, fhigh=fhigh, kwargs...)
    else
        error("fmethod must be :ridge or :threshold (got $fmethod)")
    end
end

# =====================================================================
# 4. GUI 변환  (XY, UserUnit, sequence_check! 는 GUI 쪽 정의 필요)
#   좌표계: y = 1 - f/fmax  (y 작을수록 고주파)
# =====================================================================

# 고정 대역 box (기존)
function detection_to_liness(ON, OFF, spec;
                             flow::Float64=15_000.0, fhigh::Float64=110_000.0)
    out = Vector{Vector{XY{UserUnit}}}()
    fmax = maximum(spec.f)

    y_low  = clamp(1.0 - fhigh / fmax, 0.0, 1.0)
    y_high = clamp(1.0 - flow / fmax, 0.0, 1.0)

    for (on, off) in zip(ON, OFF)
        on  = clamp(on, 1, length(spec.t_sng))
        off = clamp(off, 1, length(spec.t_sng))
        t1 = Float64(spec.t_sng[on])
        t2 = Float64(spec.t_sng[off])
        push!(out, [XY{UserUnit}(t1, y_low), XY{UserUnit}(t2, y_high)])
    end

    sequence_check!(out)
    return out
end

# 이미 계산된 (f_lo, f_hi) 로 box 생성
function bounds_to_liness(ON, OFF, bounds::AbstractVector{<:Tuple}, spec)
    length(bounds) == length(ON) == length(OFF) ||
        throw(DimensionMismatch("ON, OFF, bounds must have equal length"))
    out  = Vector{Vector{XY{UserUnit}}}()
    fmax = maximum(spec.f)
    nt   = length(spec.t_sng)

    for ((on, off), (f_lo, f_hi)) in zip(zip(ON, OFF), bounds)
        a, b = minmax(clamp(Int(on), 1, nt), clamp(Int(off), 1, nt))
        t1 = Float64(spec.t_sng[a])
        t2 = Float64(spec.t_sng[b])
        y_hi = clamp(1.0 - f_hi / fmax, 0.0, 1.0)   # 작은 y = 고주파
        y_lo = clamp(1.0 - f_lo / fmax, 0.0, 1.0)   # 큰 y  = 저주파
        push!(out, [XY{UserUnit}(t1, y_hi), XY{UserUnit}(t2, y_lo)])
    end

    sequence_check!(out)
    return out
end

# 시간 ON/OFF + 개별 주파수 경계 box. fmethod 로 :ridge / :threshold 선택.
function detection_to_liness_2d(ON, OFF, spec, sng::AbstractMatrix{<:Real};
                                fmethod::Symbol=:ridge,
                                flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                                kwargs...)
    bounds = usv_segment_freq_bounds(sng, spec.f, ON, OFF;
                                     fmethod=fmethod, flow=flow, fhigh=fhigh, kwargs...)
    return bounds_to_liness(ON, OFF, bounds, spec)
end

# =====================================================================
# 5. Pipeline: 시간축 detection + 주파수 경계를 한 번에 (GUI 불필요)
#   detect_kw : usv_whole_2d 에 전달할 NamedTuple   예) (; fp_mode=:confident)
#   freq_kw   : 주파수 방식에 전달할 NamedTuple      예) (; noise_k=2.0)
#   반환: (ON, OFF, bounds)   bounds[k] = (f_lo, f_hi) Hz
# =====================================================================
function usv_detect_2d(sng::Matrix{Float64}, dt_sng::Float64,
                       t_sng::AbstractVector, f::AbstractVector{<:Real};
                       fmethod::Symbol=:ridge,
                       flow::Float64=15_000.0, fhigh::Float64=110_000.0,
                       detect_kw::NamedTuple=(;), freq_kw::NamedTuple=(;))
    ON, OFF = usv_whole_2d(sng, dt_sng, t_sng, f; detect_kw...)
    bounds = usv_segment_freq_bounds(sng, f, ON, OFF;
                                     fmethod=fmethod, flow=flow, fhigh=fhigh, freq_kw...)
    return (; ON, OFF, bounds)
end
