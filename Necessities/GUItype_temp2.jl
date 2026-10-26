
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
