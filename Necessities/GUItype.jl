
# temporary

struct USVBox
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
end


function restore_gui_types(boxes)
    [[
        XY{UserUnit}(b.x1, b.y1),
        XY{UserUnit}(b.x2, b.y2)
    ] for b in boxes]
end

function strip_gui_types(lines::Vector{Vector{XY{UserUnit}}})
    [
        USVBox(
            Float64(l[1].x.val),
            Float64(l[1].y.val),
            Float64(l[2].x.val),
            Float64(l[2].y.val)) for l in lines
    ]
end

"""
function normalize_liness(x)
    if x isa Observable
        return normalize_liness(x[])  # unwrap
    elseif x isa AbstractVector{USVBox}
        return Observable(restore_gui_types(x))
    elseif x isa AbstractVector && (isempty(x) || first(x) isa AbstractVector{XY{UserUnit}})
        # already GUI type: Vector{Vector{XY{UserUnit}}}}
        return Observable(x)
    else
        error("Unsupported liness type: (typeof(x))")
    end
end
"""

# USVBox -> Vector{XY{UserUnit}} convert
# + y가 0..1 ratio & 1..nfreq bin case both
function normalize_liness(raw; nfreq::Int=257)
    raw === nothing && return Vector{Vector{XY{UserUnit}}}()
    raw isa Observable && (raw = raw[])
    isempty(raw) && return Vector{Vector{XY{UserUnit}}}()

    b0 = first(raw)

    # Case 1) already [ [XY,XY], [XY,XY], ... ]
    if b0 isa AbstractVector && length(b0) >= 2 && hasproperty(b0[1], :x) && hasproperty(b0[1], :y)
        return raw
    end

    # Case 2) USVBox-like: has x1,y1,x2,y2
    fn = fieldnames(typeof(b0))
    if (:x1 in fn) && (:y1 in fn) && (:x2 in fn) && (:y2 in fn)
        out = Vector{Vector{XY{UserUnit}}}(undef, length(raw))
        @inbounds for i in eachindex(raw)
            b = raw[i]
            x1 = Float64(getfield(b, :x1))
            x2 = Float64(getfield(b, :x2))
            y1 = Float64(getfield(b, :y1))
            y2 = Float64(getfield(b, :y2))

            if !(abs(y1) <= 1.1 && abs(y2) <= 1.1) && (0.5 <= y1 <= nfreq + 0.5) && (0.5 <= y2 <= nfreq + 0.5)
                y1 = (nfreq + 1 - y1) / nfreq
                y2 = (nfreq + 1 - y2) / nfreq
            end

            out[i] = [XY{UserUnit}(x1, y1), XY{UserUnit}(x2, y2)]
        end
        sequence_check!(out)
        return out
    end

    error("normalize_liness: unsupported liness element type=$(typeof(b0)), fields=$(fn)")
end
