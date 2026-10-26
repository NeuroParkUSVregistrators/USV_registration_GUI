# sdts_filter.jl
# Short-Duration Transient Suppression (SDTS) filter
#
# Based on:
#   "Spectrogram denoising and automated extraction of the fundamental
#    frequency variation of dolphin whistles"
#
# 입력 행렬 레이아웃: S[freq, time]  (spec.fout 과 동일; size(S,1) = nfreq)
#
#   out = alpha * S + (1 - alpha) * ( max(H, D1, D2) - V )
#
#   H  : 시간축으로 길쭉한 가우시안 (tonal / USV 성분 강조)
#   D1 : +45° 대각선으로 길쭉한 가우시안 (FM sweep 강조)
#   D2 : -45° 대각선으로 길쭉한 가우시안
#   V  : 주파수축으로 길쭉한 가우시안 (짧은 broadband transient = click 등)
#
# Python 원본의 _gaussian_filter / _diagonal_gaussian_filter 는
# (sigma_x, sigma_y) = (a, 6a) / (6a, a) 조합을 쓰는데, 여기서는 의미가
# 헷갈리지 않도록 "시간축 σ / 주파수축 σ / 회전각" 으로 직접 지정한다.
# 4개 커널을 한 번의 루프에서 동시에 계산하므로 중간 행렬을 만들지 않는다
# (긴 녹음에서도 메모리 = 출력 행렬 1개).

"""
    sdts_kernel(m, σ_t, σ_f; θ=0.0) -> Matrix{Float64}

m×m 정규화 가우시안 커널. 행 = 주파수 오프셋, 열 = 시간 오프셋.
θ (rad) 만큼 회전된 좌표계에서 시간방향 σ_t, 주파수방향 σ_f.
"""
function sdts_kernel(m::Int, σ_t::Real, σ_f::Real; θ::Real=0.0)
    h = m ÷ 2
    K = Matrix{Float64}(undef, m, m)
    c, s = cos(θ), sin(θ)
    for (jj, dt) in enumerate(-h:h), (ii, df) in enumerate(-h:h)
        u = c * dt + s * df      # 회전된 "시간" 축
        v = -s * dt + c * df     # 회전된 "주파수" 축
        K[ii, jj] = exp(-(u^2 / (2 * σ_t^2) + v^2 / (2 * σ_f^2)))
    end
    K ./= sum(K)
    return K
end

"""
    sdts_filter(S; alpha=0.5, m=3, a=nothing, clamp_negative=true) -> Matrix{Float64}

SDTS 필터 (out-of-place). 원본 `S`는 바뀌지 않는다.

- `alpha` : 원본 스펙트로그램을 얼마나 남길지 (0..1). 0 이면 필터 결과만.
- `m`     : 커널 크기 (홀수).
- `a`     : 커널 σ 기준값. 기본 m/10 (Python 원본과 동일).
- `clamp_negative` : `max(...) - V` 가 음수가 되는 픽셀을 0으로 자름.
  (이후 detection 이 음수 파워를 가정하지 않도록 기본 true)

경계는 replicate(가장자리 값 반복) 패딩.
"""
function sdts_filter(S::AbstractMatrix{<:Real};
    alpha::Real=0.5,
    m::Int=3,
    a::Union{Nothing,Real}=nothing,
    clamp_negative::Bool=true,
)
    isodd(m) || throw(ArgumentError("m must be an odd number (got $m)"))
    0 <= alpha <= 1 || throw(ArgumentError("alpha must be in [0, 1] (got $alpha)"))
    a = a === nothing ? m / 10 : float(a)

    Kh  = sdts_kernel(m, 6a, a)              # horizontal (time-long)
    Kv  = sdts_kernel(m, a, 6a)              # vertical   (freq-long)
    Kd1 = sdts_kernel(m, 6a, a; θ=π / 4)     # diagonal +45°
    Kd2 = sdts_kernel(m, 6a, a; θ=-π / 4)    # diagonal -45°

    nf, nt = size(S)
    h = m ÷ 2
    α = float(alpha)
    out = Matrix{Float64}(undef, nf, nt)

    # 커널은 180° 회전 대칭이라 correlation == convolution
    Threads.@threads for j in 1:nt
        @inbounds for i in 1:nf
            sh = 0.0; sv = 0.0; s1 = 0.0; s2 = 0.0
            for dj in -h:h
                jj = clamp(j + dj, 1, nt)
                kj = dj + h + 1
                for di in -h:h
                    ii = clamp(i + di, 1, nf)
                    ki = di + h + 1
                    x = Float64(S[ii, jj])
                    sh += Kh[ki, kj] * x
                    sv += Kv[ki, kj] * x
                    s1 += Kd1[ki, kj] * x
                    s2 += Kd2[ki, kj] * x
                end
            end
            y = α * Float64(S[i, j]) + (1 - α) * (max(sh, s1, s2) - sv)
            out[i, j] = clamp_negative ? max(y, 0.0) : y
        end
    end
    return out
end

"""
    sdts_filter!(S; kwargs...)

Python 원본처럼 in-place 로 쓰고 싶을 때. `S` 를 결과로 덮어쓴다.
(내부적으로는 출력 버퍼 하나를 만든 뒤 copy!)
"""
function sdts_filter!(S::AbstractMatrix{Float64}; kwargs...)
    copy!(S, sdts_filter(S; kwargs...))
    return S
end

# ---- GUI 에서 쓰는 설정 (REPL에서 바꿔가며 튜닝 가능) ----
const SDTS_ENABLED = Ref(true)
const SDTS_ALPHA   = Ref(0.5)
const SDTS_M       = Ref(7)
const SDTS_A       = Ref{Union{Nothing,Float64}}(nothing)   # nothing → m/10
