include("GUItype_temp2.jl")

# =====================================================================
# gui_function_temp.jl
#   의존: usv_detection.jl (userunit_to_freq_bin / freq_bin_to_userunit /
#         usv_harmonic_analysis). 메인 파일에서 함께 include 되어 있어야 함
#         (include 순서는 상관없음 — 호출 시점에만 필요).
#
#   2026-10 수정
#     - 주파수 좌표 규약을 usv_detection.jl 과 하나로 통일
#         y(bin) = (nfreq + 1 - bin) / nfreq ,  bin(y) = nfreq + 1 - round(y·nfreq)
#       기존 userunit_to_freq 공식(257 고정)과 같은 식이라 저장된 수동 box 와 호환.
#       → 자동 탐지 box 를 읽어 들일 때 생기던 ±1 bin 어긋남 제거, 257 하드코딩 제거.
#     - harmonic 을 포함한 box 에서 feature(purity, mean F, pitch var)가
#       harmonic 때문에 왜곡되지 않도록, 기본적으로 fundamental 대역에서만 계산
#       (FEATURE_BAND[] = :box 로 바꾸면 이전 동작).
#     - find_the_time / time_to_ind 타입 제한 완화 (AbstractVector).
#     - spec 에서 주파수축을 찾는 spec_freq_axis 추가 (f / f_sng / fs 순).
#     - load_file: 확장자 판정을 endswith 로, 알 수 없는 확장자는 명확한 에러.
#     - normalize_liness 를 이 파일로 통합 (옛 normalize_to_gui_liness 대체).
#         · load_file 에서 바로 nfreq 와 함께 호출 → bin 단위로 저장된 y 가 실제로 변환됨
#           (이전엔 normalize_to_gui_liness 가 먼저 XY 로 바꿔서 bin 변환이 한 번도 안 돌았음)
#         · y 단위(ratio/bin)를 box 마다가 아니라 파일 전체로 판정, 저장된 "y_unit" 우선
#         · x1 > x2 box 정렬, NaN/Inf box 제거, 항상 새 벡터 반환
# =====================================================================

# ---------------------------------------------------------------------
# spec 축 helper
# ---------------------------------------------------------------------

# spec_data 에서 주파수축(Hz)을 찾음. GUI / 탐지 / feature 모두 이 함수를 쓸 것.
function spec_freq_axis(spec)
    nf = size(spec.fout, 1)
    for name in (:f, :f_sng)
        if hasproperty(spec, name)
            fv = getproperty(spec, name)
            length(fv) == nf && return fv
        end
    end
    hasproperty(spec, :fs) && return range(0.0, getproperty(spec, :fs) / 2, length=nf)
    error("spec 에서 주파수축을 찾을 수 없음 (f, f_sng 또는 fs 필드 필요)")
end

# ---------------------------------------------------------------------
# 파일 로딩
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# liness 정규화 (JLD 에서 읽은 어떤 형식이든 → Vector{Vector{XY{UserUnit}}}, y 는 ratio)
#   ※ 다른 파일에 있던 normalize_liness / 이 파일의 옛 normalize_to_gui_liness,
#     point_xy, box_bounds_any 는 이 구현으로 대체. 옛 normalize_liness 정의는 지울 것.
# ---------------------------------------------------------------------

# UserUnit 같은 wrapper(.val) 와 일반 숫자를 모두 Float64 로
_num(v) = v isa Real ? Float64(v) :
          hasproperty(v, :val) ? Float64(getproperty(v, :val)) : Float64(v)

function _point(p)
    if hasproperty(p, :x) && hasproperty(p, :y)
        return _num(getproperty(p, :x)), _num(getproperty(p, :y))
    elseif (p isa Tuple || p isa AbstractVector) && length(p) >= 2
        return _num(p[1]), _num(p[2])
    end
    error("Unsupported point type: $(typeof(p))")
end

# 한 box → (x1, x2, y1, y2). USVBox 류(x1,y1,x2,y2 필드) 또는 2-point 벡터/튜플.
function _box_coords(b)
    if hasproperty(b, :x1) && hasproperty(b, :x2) && hasproperty(b, :y1) && hasproperty(b, :y2)
        return _num(getproperty(b, :x1)), _num(getproperty(b, :x2)),
               _num(getproperty(b, :y1)), _num(getproperty(b, :y2))
    elseif (b isa AbstractVector || b isa Tuple) && length(b) >= 2
        x1, y1 = _point(b[1])
        x2, y2 = _point(b[2])
        return x1, x2, y1, y2
    end
    error("Unsupported liness element: $(typeof(b))")
end

# y 단위 판정: box 하나씩이 아니라 파일 전체 기준.
#   ratio : 0..1 (캔버스 밖으로 드래그하면 약간 벗어날 수 있어 ratio_max 까지 허용)
#   bin   : 1..nfreq
function _liness_y_unit(ys::AbstractVector{<:Real}, nfreq::Int; ratio_max::Real=2.0)
    isempty(ys) && return :ratio
    maximum(abs, ys) <= ratio_max && return :ratio
    all(y -> 0.5 <= y <= nfreq + 0.5, ys) && return :bin
    error("liness 의 y 범위를 ratio 로도 bin 으로도 해석할 수 없음: " *
          "extrema=$(extrema(ys)), nfreq=$nfreq")
end

# y_unit = :auto | :ratio | :bin   (저장 시 "y_unit" 을 같이 저장했다면 그 값을 넘길 것)
function normalize_liness(raw; nfreq::Int=257, y_unit::Symbol=:auto)
    out = Vector{Vector{XY{UserUnit}}}()
    raw === nothing && return out
    raw isa Observable && (raw = raw[])
    isempty(raw) && return out

    coords = [_box_coords(b) for b in raw]
    unit = y_unit === :auto ? _liness_y_unit([c[k] for c in coords for k in (3, 4)], nfreq) : y_unit
    unit in (:ratio, :bin) || throw(ArgumentError("y_unit must be :auto, :ratio or :bin (got $y_unit)"))

    nbad = 0
    for (x1, x2, y1, y2) in coords
        if unit === :bin    # freq_bin_to_userunit 와 같은 식 (소수 bin 도 허용)
            y1 = (nfreq + 1 - y1) / nfreq
            y2 = (nfreq + 1 - y2) / nfreq
        end
        if !(isfinite(x1) && isfinite(x2) && isfinite(y1) && isfinite(y2))
            nbad += 1
            continue
        end
        push!(out, [XY{UserUnit}(min(x1, x2), y1), XY{UserUnit}(max(x1, x2), y2)])
    end
    nbad > 0 && @warn "normalize_liness: 좌표가 NaN/Inf 인 box $nbad 개를 버림"

    sequence_check!(out)
    return out
end

# 옛 이름 호환 (nfreq 를 반드시 넘길 것 — 기본 257 은 fout 행 수가 다르면 틀림)
normalize_to_gui_liness(raw; nfreq::Int=257, y_unit::Symbol=:auto) =
    normalize_liness(raw; nfreq=nfreq, y_unit=y_unit)

function load_liness(nm::AbstractString)
    d = JLD.load(nm)
    if haskey(d, "liness")
        return d["liness"]
    elseif !isempty(keys(d))
        k = first(keys(d))
        @warn "No 'liness' key in JLD. Using first key = $k"
        return d[k]
    else
        @warn "JLD has no keys : $nm"
        return nothing
    end
end

# 저장할 때 "y_unit" => "ratio" 를 같이 넣어두면 자동 판정 대신 그 값을 사용
function load_liness_unit(nm::AbstractString)
    d = JLD.load(nm)
    haskey(d, "y_unit") || return :auto
    u = Symbol(d["y_unit"])
    return u in (:ratio, :bin) ? u : :auto
end

function load_file(nm::String, idx::Int)
    lname = lowercase(nm)
    if endswith(lname, ".h5")              # first time
        nmh5 = nm
    elseif endswith(lname, ".jld")         # second time, already created liness
        nmh5 = string(nm[1:findlast('.', nm)-1], ".h5")
        isfile(nmh5) || error("짝이 되는 h5 파일이 없음: $nmh5")
    else
        error("지원하지 않는 파일 형식 (.h5 또는 .jld): $nm")
    end

    anl, dig = loadNIDAQdataFromh5(nmh5)
    sp = stft_param()
    audio = anl[idx]
    spec = spec_data(audio, sp)

    # spec 을 먼저 만들어야 nfreq 를 알 수 있음 → 그다음 liness 정규화
    liness_gui = if endswith(lname, ".jld")
        normalize_liness(load_liness(nm); nfreq=size(spec.fout, 1), y_unit=load_liness_unit(nm))
    else
        Vector{Vector{XY{UserUnit}}}()
    end

    return spec, Observable(liness_gui)
end

# ---------------------------------------------------------------------
# feature 계산
# ---------------------------------------------------------------------

# feature 를 계산할 주파수 대역
#   :fundamental (기본) harmonic syllable 이면 box 안에서 fundamental 대역만 사용
#                        (harmonic 이 없는 syllable 은 box 전체 → 이전과 동일)
#   :box                box 전체 (이전 동작)
const FEATURE_BAND = Ref(:fundamental)
# usv_harmonic_analysis 에 넘길 추가 옵션 (예: (; peak_rel=0.05)). harmonic_box 는 넣지 말 것.
const FEATURE_HARM_KW = Ref{NamedTuple}((;))

# box (x1:x2, ylo:yhi) 에서 feature 계산에 쓸 row 범위
function _feature_rows(spec, x1::Int, x2::Int, ylo::Int, yhi::Int; band::Symbol=FEATURE_BAND[])
    band === :box && return ylo, yhi
    band === :fundamental || throw(ArgumentError("band must be :fundamental or :box (got $band)"))
    yhi - ylo < 2 && return ylo, yhi

    f = spec_freq_axis(spec)
    a = usv_harmonic_analysis(spec.fout, f, x1, x2;
                              flow=Float64(f[ylo]), fhigh=Float64(f[yhi]),
                              FEATURE_HARM_KW[]..., harmonic_box=:fundamental)
    a.is_harmonic || return ylo, yhi

    lo = clamp(searchsortedfirst(f, a.fund_lo), ylo, yhi)
    hi = clamp(searchsortedlast(f, a.fund_hi), ylo, yhi)
    return lo <= hi ? (lo, hi) : (ylo, yhi)
end

# 한 syllable 의 (dur[ms], dB, purity, meanF, pitchV)
function _syllable_features(spec::spec_data, box::Vector{XY{UserUnit}};
                            room_noise::Float64=0.16, band::Symbol=FEATURE_BAND[])
    t1, t2 = box[1].x.val, box[2].x.val
    dur = (t2 - t1) * 1000   # in ms

    x1, x2 = minmax(time_to_ind([t1, t2], spec.t_sng)...)
    ylo, yhi = userunit_to_freq(box, size(spec.fout, 1))
    ylo, yhi = _feature_rows(spec, x1, x2, ylo, yhi; band=band)

    syllable_sng = @view spec.fout[ylo:yhi, x1:x2]
    purity, amp = spectral(syllable_sng)
    dB = mean(power_to_decibels.(amp, room_noise))
    meanF, pitchV = frequency(syllable_sng, spec, ylo)

    return dur, dB, mean(purity), meanF, pitchV
end

function gui_duration!(dur::Vector{Float64}, lns::Vector{Vector{XY{UserUnit}}})
    for i in eachindex(lns)
        dur[i] = (lns[i][2].x.val - lns[i][1].x.val) * 1000 # in ms
    end
    nothing
end

function gui_usvfeatures!(spec::spec_data,
                          lns::Vector{Vector{XY{UserUnit}}},
                          dur::Vector{Float64},
                          dB::Vector{Float64}, spur::Vector{Float64},
                          mf::Vector{Float64}, pv::Vector{Float64};
                          room_noise::Float64=0.16, band::Symbol=FEATURE_BAND[])
    for i in eachindex(lns)
        dur[i], dB[i], spur[i], mf[i], pv[i] =
            _syllable_features(spec, lns[i]; room_noise=room_noise, band=band)
    end
    nothing
end

function single_usvfeatures(spec::spec_data,
                            lns::Vector{Vector{XY{UserUnit}}},
                            usvind::Int64;
                            room_noise::Float64=0.16, band::Symbol=FEATURE_BAND[])
    return _syllable_features(spec, lns[usvind]; room_noise=room_noise, band=band)
end

# ---------------------------------------------------------------------
# 좌표 변환
# ---------------------------------------------------------------------

function time_to_ind(time::Vector{Float64}, ref::AbstractVector{<:Real})
    ind, = find_the_time(time, ref)
    return ind
end

# box 의 y(UserUnit) → 주파수 bin [low, high]
#   nfreq 는 size(spec.fout, 1). 기본값 257 은 기존 호출 호환용.
#   공식은 usv_detection.jl 의 userunit_to_freq_bin 과 공유 (기존 257 공식과 동일).
function userunit_to_freq(user::Vector{XY{UserUnit}}, nfreq::Int=257)
    a = userunit_to_freq_bin(user[1].y.val, nfreq)
    b = userunit_to_freq_bin(user[2].y.val, nfreq)
    return [min(a, b), max(a, b)]
end

function userunit_to_freq(users::Vector{Vector{XY{UserUnit}}}, nfreq::Int=257)
    [userunit_to_freq(us, nfreq) for us in users]
end

function unit_to_time(unit::XY{UserUnit}, ctime::Float64, tspan::Float64)
    XY{UserUnit}(ctime + unit.x.val * tspan, unit.y.val)
end

function unit_to_time(unit::Vector{XY{UserUnit}}, ctime::Float64, tspan::Float64)
    # time = current time + userunit*timespan
    return [XY{UserUnit}(ctime + unit[1].x.val * tspan, unit[1].y.val),
            XY{UserUnit}(ctime + unit[2].x.val * tspan, unit[2].y.val)]
end

function time_to_unit(time::Vector{XY{UserUnit}}, ctime::Float64, tspan::Float64)
    # current position = (time - current time) / timespan
    return [XY{UserUnit}((time[1].x.val - ctime) / tspan, time[1].y.val),
            XY{UserUnit}((time[2].x.val - ctime) / tspan, time[2].y.val)]
end

function time_to_unit(times::Vector{Vector{XY{UserUnit}}}, ctime::Float64, tspan::Float64)
    [time_to_unit(ts, ctime, tspan) for ts in times]
end

# sometimes usv syllable onset pairs can be mixed, so to make sure they are in time sequence orders
function sequence_check!(ut::Vector{Vector{XY{UserUnit}}})
    sort!(ut, by=x -> (x[1].x.val))
    nothing
end

# find the index of each time point
#   tp = time point, ta = time array (정렬되어 있어야 함). 반환 (index, time)
#   ※ 메인 GUI 파일에 있는 같은 이름의 StepRangeLen 전용 정의는 지울 것
#     (더 좁은 타입이라 그쪽이 우선 호출되고, 메서드 덮어쓰기 경고가 남).
function find_the_time(tp::Real, ta::AbstractVector{<:Real})
    isempty(ta) && throw(ArgumentError("time axis is empty"))
    i0 = firstindex(ta)
    isfinite(tp) || return i0, ta[i0]
    ti = clamp(searchsortedfirst(ta, tp), i0, lastindex(ta))
    return ti, ta[ti]
end

function find_the_time(temp::AbstractVector{<:Real}, ta::AbstractVector{<:Real})
    a = zeros(Int, length(temp))
    b = zeros(Float64, length(temp))
    for (i, tp) in enumerate(temp)
        a[i], b[i] = find_the_time(tp, ta)
    end
    return a, b
end
