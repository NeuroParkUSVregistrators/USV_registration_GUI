using Pkg

function stretcher(dur, max_dur)
    #x1, x2, y1, y2 = box_bounds(box)
    #dur = x2 - x1
    if dur < max_dur
        dur = sqrt(dur * max_dur)
    end
    return dur
end # stretching time

#how to ...