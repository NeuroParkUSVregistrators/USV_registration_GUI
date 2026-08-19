# JPEG image saving
# 2026-02-23

# Collecting Full-band range
# White noise padding (Padval removed)

using Pkg
using JLD
using HDF5
#Pkg.add("Images")
#Pkg.add("ImageIO")
#Pkg.add("FileIO")
#Pkg.add("ImageMagick")
#Pkg.add("ImageTransformations") 
#Pkg.add("Interpolations")
using Images
using ImageMagick
using FileIO
using Printf
using Statistics
using ImageTransformations  # imresize
using Interpolations

include("wavesurfer.jl")
include("acoustic.jl")
include("gui_function_temp.jl")
include("GUItype.jl")
include("Stretcher.jl")

function sequence_check_any!(ut)
    sort!(ut, by=b -> box_bounds(b)[1])  # tstart
    return nothing
end

# Try to turn whatever "point" is into (x::Float64, y::Float64)
getval(v) = v isa Number ? Float64(v) :
            (hasproperty(v, :val) ? Float64(getproperty(v, :val)) : Float64(v))

function point_xy(p)::Tuple{Float64,Float64}
    if p isa Tuple && length(p) >= 2
        return (Float64(p[1]), Float64(p[2]))
    elseif p isa AbstractVector && length(p) >= 2
        return (Float64(p[1]), Float64(p[2]))
    elseif hasproperty(p, :x) && hasproperty(p, :y)
        return (getval(getproperty(p, :x)), getval(getproperty(p, :y)))
    else
        error("Unsupported point type: $(typeof(p))")
    end
end

# USVBox (x1,y1,x2,y2)
function box_bounds(box)::Tuple{Float64,Float64,Float64,Float64}
    if hasproperty(box, :x1) && hasproperty(box, :x2) && hasproperty(box, :y1) && hasproperty(box, :y2)
        x1 = getval(getproperty(box, :x1))
        x2 = getval(getproperty(box, :x2))
        y1 = getval(getproperty(box, :y1))
        y2 = getval(getproperty(box, :y2))
        return (min(x1, x2), max(x1, x2), y1, y2)
    end
    if hasproperty(box, :t1) && hasproperty(box, :t2) && hasproperty(box, :y1) && hasproperty(box, :y2)
        x1 = getval(getproperty(box, :t1))
        x2 = getval(getproperty(box, :t2))
        y1 = getval(getproperty(box, :y1))
        y2 = getval(getproperty(box, :y2))
        return (min(x1, x2), max(x1, x2), y1, y2)
    end
    if hasproperty(box, :start) && hasproperty(box, :stop) && hasproperty(box, :low) && hasproperty(box, :high)
        x1 = getval(getproperty(box, :start))
        x2 = getval(getproperty(box, :stop))
        y1 = getval(getproperty(box, :low))
        y2 = getval(getproperty(box, :high))
        return (min(x1, x2), max(x1, x2), y1, y2)
    end
    #error("Unsupported box type: $(typeof(box)); fields=$(fieldnames(typeof(box)))")
end

# ratio(0..1) -> bin (1..nfreq)
@inline function _ratio_to_bin(y::Float64, nfreq::Int)
    a = clamp(y, 0.0, 1.0)
    b = -Int(round(a * nfreq)) + (nfreq + 1)  # y=0 -> nfreq+1 -> clamp -> nfreq, y=1 -> 1
    return clamp(b, 1, nfreq)
end

# Hz -> bin (freq_axis가 있으면 그 축에 맞춰 매핑)
@inline function _hz_to_bin(yhz::Float64, freq_axis::AbstractVector{<:Real})
    i = searchsortedfirst(freq_axis, yhz)
    return clamp(i, firstindex(freq_axis), lastindex(freq_axis))
end

function y_to_bins_ratio_clamped(y1::Real, y2::Real, nfreq::Int)
    a = _ratio_to_bin(y1, nfreq)
    b = _ratio_to_bin(y2, nfreq)
    tmp = a
    a > b && (a=b, b=tmp)
    return (a, b)
end

# normalization / image converting
@inline function to_gray_img(mat; normalize::Symbol, gmin::Float64=0.0, gmax::Float64=1.0)
    m = Float32.(mat) # for JPEG saving
    if normalize === :local
        mn = minimum(m)
        mx = maximum(m)
        if mx <= mn
            fill!(m, 0f0)
        else
            @inbounds m .= (m .- mn) ./ (mx - mn)
        end
    elseif normalize === :global
        if gmax <= gmin
            fill!(m, 0f0)
        else
            @inbounds m .= (m .- Float32(gmin)) ./ Float32(gmax - gmin)
        end
        @inbounds m .= clamp.(m, 0f0, 1f0)
    else
        throw(ArgumentError("normalize must be :local or :global"))
    end
    return Gray.(m)
end

function normalize01(mat::AbstractArray{<:Real};
    contrast::Symbol=:local,
    gmin::Float64=0.0,
    gmax::Float64=1.0,
    qlo::Float64=0.02,
    qhi::Float64=0.98,
    gamma::Float64=1.0
)

    if isempty(mat)
        return fill(0f0, size(mat))
    end

    m = Float32.(mat)

    if contrast === :quantile
        v = vec(Float64.(m))
        vf = filter(isfinite, v)
        isempty(vf) && return fill(0f0, size(m))
        qlo2 = clamp(qlo, 0.0, 1.0)
        qhi2 = clamp(qhi, 0.0, 1.0)
        qhi2 <= qlo2 && (qhi2 = min(1.0, qlo2 + 1e-3))
        lo = Float32(quantile(vf, qlo2))
        hi = Float32(quantile(vf, qhi2))
    elseif contrast === :local
        lo = Float32(minimum(m))
        hi = Float32(maximum(m))
    elseif contrast === :global
        lo = Float32(gmin)
        hi = Float32(gmax)
    else
        throw(ArgumentError("contrast must be :local, :global, or :quantile"))
    end

    if !(hi > lo)
        fill!(m, 0f0)
    else
        @inbounds m .= (m .- lo) ./ (hi - lo)
        @inbounds m .= clamp.(m, 0f0, 1f0)
    end

    if gamma != 1.0
        invg = Float32(1.0 / gamma)
        @inbounds m .= m .^ invg
    end

    return m  # Float32, 0..1
end

function export_usv_jpegs(nmjld::String;
    outdir::String=SAVE_DIR,
    pad_ms::Float64=0.0, # time padding (ms)
    pad_freqbins::Int=12, # freq padding (bin)
    normalize::Symbol,
    quality::Int=92,
    flip_freq::Bool=true, #GUI: reverse!(dims=1), so flipping again
    min_w::Int=2,
    min_h::Int=2
)
    mkpath(outdir)

    d = load(nmjld)
    @assert haskey(d, "liness") "No 'liness' key in JLD file"
    lns = d["liness"]
    sequence_check_any!(lns)

    dotloc = findlast('.', nmjld)
    nmh5 = string(nmjld[1:dotloc-1], ".h5")

    anl, dig = loadNIDAQdataFromh5(nmh5)
    sp = stft_param()
    audio = anl[1]                 # 1/2
    spec = spec_data(audio, sp)    # spec.fout, spec.t_sng

    fout = spec.fout
    nfreq, ntime = size(fout)

    # global min/max
    gmin, gmax = 0.0, 1.0
    if normalize === :global
        gmin = minimum(fout)
        gmax = maximum(fout)
    end

    # pad_ms to time index padding
    dt = step(spec.t_sng)
    pad_t = pad_ms <= 0 ? 0 : Int(ceil((pad_ms / 1000) / dt))

    println("Boxes: ", length(lns))
    println("fout size: ", size(fout))
    println("dt: ", dt, " sec   pad_t: ", pad_t, " frames")
    println("Export -> ", outdir)

    # Crop out each boxes
    for (i, box) in enumerate(lns)
        (tstart, tend, yy1, yy2) = try
            box_bounds(box)
        catch e
            @warn "skip box $i: box_bounds failed" exception = (e, catch_backtrace())
            continue
        end

        # time index (x)
        x = time_to_ind([tstart, tend], spec.t_sng)
        x1 = clamp(x[1] - pad_t, 1, ntime)
        x2 = clamp(x[2] + pad_t, 1, ntime)
        tmp = x1
        x1 > x2 && (x1=x2, x2=tmp)

        # freq index (y)
        (y1, y2) = y_to_bins_ratio_clamped(yy1, yy2, nfreq)
        y1 = clamp(y1 - pad_freqbins, 1, nfreq)
        y2 = clamp(y2 + pad_freqbins, 1, nfreq)
        tmp = y1
        y1 > y2 && (y1=y2, y2=tmp)

        # minimum size ensurance
        if (x2 - x1 + 1) < min_w
            x2 = min(ntime, x1 + (min_w - 1))
        end
        if (y2 - y1 + 1) < min_h
            c = (y1 + y2) ÷ 2
            half = (min_h - 1) ÷ 2
            y1 = clamp(c - half, 1, nfreq)
            y2 = clamp(y1 + min_h - 1, 1, nfreq)
        end

        # crop (view)
        cut = @view fout[y1:y2, x1:x2]

        # flipping freq axis (goes w/ update_sng)
        mat = flip_freq ? reverse(cut; dims=1) : Array(cut)
        img = to_gray_img(mat; normalize=normalize, gmin=gmin, gmax=gmax)

        # index + time + freqbin range
        fname = joinpath(outdir, @sprintf("%04d_%.4f-%.4f_y%03d-%03d.jpeg", i, tstart, tend, y1, y2))

        try
            Images.save(fname, img; quality=quality)
        catch
            Images.save(fname, img)
        end
    end

    println("Done.")
    return nothing
end


function export_usv_jpegs_fullband(
    nmjld::String;
    outdir::String,
    pad_ms::Float64=0.0,
    flip_freq::Bool=true,
    quality::Int=92,

    # output
    outH::Int=128,
    outW::Int=128,
    padval::Float32=0f0,
    keep_aspect::Bool=true,
    min_w::Int=2,

    # contrast
    contrast::Symbol=:quantile,
    qlo::Float64=0.02,
    qhi::Float64=0.98,
    gamma::Float64=1.0, idx::Int=1
)
    mkpath(outdir)

    d = load(nmjld)
    @assert haskey(d, "liness") "No 'liness' key in JLD file"
    lns = d["liness"]
    sort!(lns, by=b -> min(getfield(b, :x1), getfield(b, :x2)))  # time sort

    dotloc = findlast('.', nmjld)
    nmh5 = string(nmjld[1:dotloc-1], ".h5")

    anl, dig = loadNIDAQdataFromh5(nmh5)
    sp = stft_param()
    audio = anl[idx]
    spec = spec_data(audio, sp)

    fout = spec.fout
    nfreq, ntime = size(fout)

    # global min/max (contrast=:global)
    gmin = minimum(fout)
    gmax = maximum(fout)

    dt = step(spec.t_sng)
    pad_t = pad_ms <= 0 ? 0 : Int(ceil((pad_ms / 1000) / dt))

    println("Boxes: ", length(lns))
    println("fout size: ", size(fout))
    println("dt: ", dt, " sec   pad_t: ", pad_t, " frames")
    println("Export -> ", outdir)

    for (i, box) in enumerate(lns)
        # USVBox: x1,x2
        tstart = min(getfield(box, :x1), getfield(box, :x2))
        tend = max(getfield(box, :x1), getfield(box, :x2))

        # time index
        x = time_to_ind([tstart, tend], spec.t_sng)
        x1 = clamp(x[1] - pad_t, 1, ntime)
        x2 = clamp(x[2] + pad_t, 1, ntime)
        tmp = x1
        x1 > x2 && (x1=x2, x2=tmp)

        if (x2 - x1 + 1) < min_w
            x2 = min(ntime, x1 + (min_w - 1))
        end

        if x1 > x2
            @warn "skip box $i: invalid crop range" x1 = x1 x2 = x2
            continue
        end

        if (x2 - x1 + 1) < 2
            x1 = clamp(x1, 1, ntime - 1)
            x2 = x1 + 1
        end

        cut = @view fout[:, x1:x2]
        mat = flip_freq ? reverse(cut; dims=1) : Array(cut)

        # normalize/contrast -> 0..1
        mat01 = normalize01(mat; contrast=contrast, gmin=gmin, gmax=gmax, qlo=qlo, qhi=qhi, gamma=gamma)

        # resize + pad to 128x128
        out01 = resize_pad(mat01; outH=outH, outW=outW, qlo=qlo, qhi=qhi, keep_aspect=keep_aspect)

        img = to_gray(out01)

        fname = joinpath(outdir, @sprintf("%04d_%.4f-%.4f_fullband.jpeg", i, tstart, tend))
        try
            Images.save(fname, img; quality=quality)
        catch
            Images.save(fname, img)
        end
    end

    println("Done.")
    return nothing
end

#Cutoff 30-110kHz, according to the Goffinet et al. / stretch 적용 전. (짧은 것만 tmax로 늘렸다고 서술되어 있음.: slight gain reported)
function export_usv_jpegs_cutoff(
    nmjld::String;
    outdir::String,
    pad_ms::Float64=0.0,
    flip_freq::Bool=true,
    quality::Int=92,

    # output
    outH::Int=128,
    outW::Int=128,
    padval::Float32=0f0,
    keep_aspect::Bool=true,
    min_w::Int=2,

    # contrast
    contrast::Symbol=:quantile,
    qlo::Float64=0.02,
    qhi::Float64=0.98,
    gamma::Float64=1.0, idx::Int=1
)
    mkpath(outdir)

    d = load(nmjld)
    @assert haskey(d, "liness") "No 'liness' key in JLD file"
    lns = d["liness"]
    sort!(lns, by=b -> min(getfield(b, :x1), getfield(b, :x2)))  # time sort

    dotloc = findlast('.', nmjld)
    nmh5 = string(nmjld[1:dotloc-1], ".h5")

    anl, dig = loadNIDAQdataFromh5(nmh5)
    sp = stft_param()
    audio = anl[idx]
    spec = spec_data(audio, sp)

    fout = spec.fout
    freq_axis = collect(spec.f)

    f_lo = 30_000.0
    f_hi = 120_000.0
    #Goffinet et al. 에서는 30-110kHz, 짧은 syllable stretch, zero padding

    y1_band = clamp(searchsortedfirst(freq_axis, f_lo), 1, length(freq_axis))
    y2_band = clamp(searchsortedlast(freq_axis, f_hi), 1, length(freq_axis))
    nfreq, ntime = size(fout)

    # global min/max (contrast=:global)
    gmin = minimum(fout)
    gmax = maximum(fout)

    dt = step(spec.t_sng)
    pad_t = pad_ms <= 0 ? 0 : Int(ceil((pad_ms / 1000) / dt))

    println("Boxes: ", length(lns))
    println("fout size: ", size(fout))
    println("dt: ", dt, " sec   pad_t: ", pad_t, " frames")
    println("Export -> ", outdir)

    for (i, box) in enumerate(lns)
        # USVBox: x1,x2
        tstart = min(getfield(box, :x1), getfield(box, :x2))
        tend = max(getfield(box, :x1), getfield(box, :x2))

        # time index
        x = time_to_ind([tstart, tend], spec.t_sng)
        x1 = clamp(x[1] - pad_t, 1, ntime)
        x2 = clamp(x[2] + pad_t, 1, ntime)
        tmp = x1
        x1 > x2 && (x1=x2, x2=tmp)

        if (x2 - x1 + 1) < min_w
            x2 = min(ntime, x1 + (min_w - 1))
        end

        if x1 > x2
            @warn "skip box $i: invalid crop range" x1 = x1 x2 = x2
            continue
        end

        if (x2 - x1 + 1) < 2
            x1 = clamp(x1, 1, ntime - 1)
            x2 = x1 + 1
        end

        cut = @view fout[y1_band:y2_band, x1:x2]
        mat = flip_freq ? reverse(cut; dims=1) : Array(cut)

        # normalize/contrast -> 0..1
        mat01 = normalize01(mat; contrast=contrast, gmin=gmin, gmax=gmax, qlo=qlo, qhi=qhi, gamma=gamma)

        # resize + pad to 128x128
        out01 = resize_pad(mat01; outH=outH, outW=outW, qlo=qlo, qhi=qhi, keep_aspect=keep_aspect)

        img = to_gray(out01)

        fname = joinpath(outdir, @sprintf("%04d_%.4f-%.4f_cutoff.jpeg", i, tstart, tend))
        try
            Images.save(fname, img; quality=quality)
        catch
            Images.save(fname, img)
        end
    end

    println("Done.")
    return nothing
end

function export_usv_jpegs_lowcut(
    nmjld::String;
    outdir::String,
    pad_ms::Float64=0.0,
    flip_freq::Bool=true,
    quality::Int=92,

    # output
    outH::Int=128,
    outW::Int=128,
    padval::Float32=0f0,
    keep_aspect::Bool=true,
    min_w::Int=2,

    # contrast
    contrast::Symbol=:quantile,
    qlo::Float64=0.02,
    qhi::Float64=0.98,
    gamma::Float64=1.0, idx::Int=1
)
    mkpath(outdir)

    d = load(nmjld)
    @assert haskey(d, "liness") "No 'liness' key in JLD file"
    lns = d["liness"]
    sort!(lns, by=b -> min(getfield(b, :x1), getfield(b, :x2)))  # time sort

    dotloc = findlast('.', nmjld)
    nmh5 = string(nmjld[1:dotloc-1], ".h5")

    anl, dig = loadNIDAQdataFromh5(nmh5)
    sp = stft_param()
    audio = anl[idx]
    spec = spec_data(audio, sp)

    fout = spec.fout
    freq_axis = collect(spec.f)

    if length(spec.f) != size(fout, 1)
        fout = permutedims(fout)
        freq_axis = freq_axis
    end

    f_lo = 30_000.0
    f_hi = 125_000.0
    #Goffinet et al. 에서는 30-110kHz, 짧은 syllable stretch, zero padding

    y2_band = clamp(searchsortedlast(freq_axis, f_hi), 1, length(freq_axis))
    y1_band = clamp(searchsortedfirst(freq_axis, f_lo), 1, length(freq_axis))
    nfreq, ntime = size(fout)

    # global min/max (contrast=:global)
    gmin = minimum(fout)
    gmax = maximum(fout)

    dt = step(spec.t_sng)
    pad_t = pad_ms <= 0 ? 0 : Int(ceil((pad_ms / 1000) / dt))

    println("Boxes: ", length(lns))
    println("fout size: ", size(fout))
    println("dt: ", dt, " sec   pad_t: ", pad_t, " frames")
    println("Export -> ", outdir)

    for (i, box) in enumerate(lns)
        # USVBox: x1,x2
        tstart = min(getfield(box, :x1), getfield(box, :x2))
        tend = max(getfield(box, :x1), getfield(box, :x2))

        # time index
        x = time_to_ind([tstart, tend], spec.t_sng)
        x1 = clamp(x[1] - pad_t, 1, ntime)
        x2 = clamp(x[2] + pad_t, 1, ntime)
        tmp = x1
        x1 > x2 && (x1=x2, x2=tmp)

        if (x2 - x1 + 1) < min_w
            x2 = min(ntime, x1 + (min_w - 1))
        end

        if x1 > x2
            @warn "skip box $i: invalid crop range" x1 = x1 x2 = x2
            continue
        end

        if (x2 - x1 + 1) < 2
            x1 = clamp(x1, 1, ntime - 1)
            x2 = x1 + 1
        end

        if y1_band > y2_band
            y1_band, y2_band = y2_band, y1_band
        end
        cut = @view fout[y1_band:y2_band, x1:x2]
        mat = flip_freq ? reverse(cut; dims=1) : Array(cut)

        # normalize/contrast -> 0..1
        mat01 = normalize01(mat; contrast=contrast, gmin=gmin, gmax=gmax, qlo=qlo, qhi=qhi, gamma=gamma)

        # resize + pad to 128x128
        out01 = resize_pad(mat01, dt; outH=outH, outW=outW, qlo=qlo, qhi=qhi, keep_aspect=keep_aspect)

        img = to_gray(out01)

        fname = joinpath(outdir, @sprintf("%04d_%.4f-%.4f_lowcut.jpeg", i, tstart, tend))
        try
            Images.save(fname, img; quality=quality)
        catch
            Images.save(fname, img)
        end
    end

    println("Done.")
    return nothing
end
# 30kHz 이상만 trim; 2026-03-30
# Stretcher로 0.2초 이하인 signal sqrt(dur*dur_max)로 늘림


function resize_pad(
    mat01::AbstractArray{<:Real},
    dt::Real;
    outH::Int=128,
    outW::Int=128,
    qlo::Float64=0.02,
    qhi::Float64=0.98,
    keep_aspect::Bool=true,
)
    H, W = size(mat01)

    # degenerate safety
    if H < 1 || W < 1
        return clamp.(noise_sigma .* randn(Float32, outH, outW) .+ 0.1f0, 0f0, 0.5f0)
    end

    if keep_aspect
        s = min(outH / H, outW / W)
        newH = max(1, Int(round(H * s)))
        newW = max(1, Int(round(W * s)))
    else
        s = min(outH / H, outW / W)
        newH = max(1, Int(round(H * s)))

        val = stretcher(W, dt * 200)
        if !isfinite(val) || val <= 0
            val = W
        end
        newW = clamp(max(1, Int(round(val))), 1, outW)
    end

    src = Float32.(mat01)

    method = if H == 1 && W == 1
        (NoInterp(), NoInterp())
    elseif H == 1
        (NoInterp(), BSpline(Linear()))
    elseif W == 1
        (BSpline(Linear()), NoInterp())
    else
        BSpline(Linear())
    end

    # imresize expects array, returns Float-like
    resized = imresize(src, (newH, newW); method=method)

    v = vec(Float64.(resized))
    vf = filter(isfinite, v)
    bg = isempty(vf) ? qlo : Float32(quantile(vf, 0.30)) #background (30%)
    #hiN = isempty(vf) ? qhi : Float32(quantile(vf, 0.55)) #cap_hi (55%)

    noise_sigma = Float32(quantile(vf, 0.65)) #0.08f0 이전 값, 어두움

    σ = noise_sigma * max(qhi - bg, 1f-3)
    out = bg .+ σ .* randn(Float32, outH, outW)
    out .= clamp.(out, 0f0, qhi)

    y0 = (outH - newH) ÷ 2 + 1
    x0 = (outW - newW) ÷ 2 + 1
    out[y0:y0+newH-1, x0:x0+newW-1] .= resized
    return out
end
to_gray(m01::AbstractArray{<:Real}) = Gray.(Float32.(m01))
