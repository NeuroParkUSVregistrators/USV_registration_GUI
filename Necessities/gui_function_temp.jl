include("GUItype.jl")
include("wavesurfer.jl")

function point_xy(p)::Tuple{Float64, Float64}
    if p isa Tuple && length(p) >= 2
        return (p[1], p[2])
    elseif p isa AbstractVector && length(p) >= 2
        return (p[1], p[2])
    elseif hasproperty(p,:x) && hasproperty(p,:y)
        return (getproperty(p,:x), getproperty(p,:y))
    else
        error("Unsupported point type: $(typeof(p))")
    end
end

function box_bounds_any(box)::Tuple{Float64,Float64,Float64,Float64} 
    # USVBox
    if hasproperty(box, :x1) && hasproperty(box, :x2) && hasproperty(box, :y1) && hasproperty(box, :y2) 
        x1 = getproperty(box, :x1)
        x2 = getproperty(box, :x2)
        y1 = getproperty(box, :y1)
        y2 = getproperty(box, :y2)
        return (min(x1,x2), max(x1,x2), y1, y2) 
    end 
    error("Unsupported liness element: $(typeof(box))") 
end

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
        return d
    end
end

function normalize_to_gui_liness(raw) 
    out = Vector{Vector{XY{UserUnit}}}() 
 
    raw = raw isa Observable ? raw[] : raw 
 
    isempty(raw) && return out 
 
    # 1) Vector{Vector{XY{UserUnit}}}
    if raw isa Vector{<:Vector{<:Any}} 
        for b in raw 
            if length(b) >= 2 && (b[1] isa XY{UserUnit}) && (b[2] isa XY{UserUnit}) 
                push!(out, [b[1], b[2]]) 
            else 
                x1,x2,y1,y2 = box_bounds_any(b) 
                push!(out, [XY{UserUnit}(x1,y1), XY{UserUnit}(x2,y2)]) 
            end 
        end 
        return out 
    end 
 
    # 2) Vector{USVBox}
    if raw isa AbstractVector 
        for box in raw 
            x1,x2,y1,y2 = box_bounds_any(box) 
            push!(out, [XY{UserUnit}(x1,y1), XY{UserUnit}(x2,y2)]) 
        end 
        return out 
    end 
 
    error("Unsupported liness container: $(typeof(raw))") 
end 

function load_file(nm::String, idx::Int)
        # fout and liness
        if occursin(".h5",nm) # first time
            # fout
            anl,dig = loadNIDAQdataFromh5(nm)
            sp = stft_param()
            #audio = anl[2]
            audio = anl[idx]
            spec = spec_data(audio,sp)
            
            # create a new liness as an observable
            liness_gui = Vector{Vector{XY{UserUnit}}}(undef,0)
    
        elseif occursin(".jld",nm) # second time, already created liness
            # fout
            dotloc = findlast('.',nm)
            nmh5 = string(nm[1:dotloc-1],".h5")
            anl,dig = loadNIDAQdataFromh5(nmh5)
            sp = stft_param()
            #audio = anl[2]
            audio = anl[idx]
            spec = spec_data(audio,sp)

            # bounding box
            #exlines = JLD.load(nm)
            #liness = exlines["liness"]
            
            raw_liness = load_liness(nm)
            liness_gui = normalize_to_gui_liness(raw_liness)

        end

    return spec, Observable(liness_gui)
end

function gui_duration!(dur::Vector{Float64},lns::Vector{Vector{XY{UserUnit}}})
    for i in 1:length(lns)
        dur[i] = (lns[i][2].x - lns[i][1].x) * 1000 # in ms
    end
    nothing
end

function gui_usvfeatures!(spec::spec_data,
                            lns::Vector{Vector{XY{UserUnit}}},
                            dur::Vector{Float64},
                            dB::Vector{Float64},spur::Vector{Float64},
                            mf::Vector{Float64},pv::Vector{Float64},
                            ;room_noise::Float64=0.16)
    # duration
    for i in 1:length(lns)
        # duration
        dur[i] = (lns[i][2].x - lns[i][1].x) * 1000 # in ms

        # rest
        # x and y for boundaries
        x = time_to_ind([lns[i][1].x.val,lns[i][2].x.val],spec.t_sng)
    
        y = userunit_to_freq(lns[i])

        # single sng 
        syllable_sng = @view spec.fout[y[1]:y[2],x[1]:x[2]]

        purity, amp = spectral(syllable_sng)
        spur[i] = mean(purity)
        dB[i] = mean(power_to_decibels.(amp,room_noise))
        meanF,pitchV = frequency(syllable_sng,spec,y[1])

        mf[i] = meanF
        pv[i] = pitchV
        
    end
    nothing
end

function single_usvfeatures(spec::spec_data,
                            lns::Vector{Vector{XY{UserUnit}}},
                            usvind::Int64
                            ;room_noise::Float64=0.16)
        # duration
        dur = (lns[usvind][2].x - lns[usvind][1].x) * 1000 # in ms

        # rest
        # x and y for boundaries
        x = time_to_ind([lns[usvind][1].x.val,lns[usvind][2].x.val],spec.t_sng)
    
        y = userunit_to_freq(lns[usvind])

        # single sng 
        syllable_sng = @view spec.fout[y[1]:y[2],x[1]:x[2]]

        purity, amp = spectral(syllable_sng)
        spur = mean(purity)
        dB = mean(power_to_decibels.(amp,room_noise))
        meanF,pitchV = frequency(syllable_sng,spec,y[1])

    return dur, dB, spur, meanF, pitchV
end

function time_to_ind(time::Vector{Float64},ref::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64})
    ind, = find_the_time(time,ref)
    return ind
end
    
function userunit_to_freq(user::Vector{XY{UserUnit}})
    
    fq = 257
    # value in ratio
    low = Int(round(user[1].y.val*fq))
    high = Int(round(user[2].y.val*fq))

    low *=-1
    low += 258

    high *=-1
    high +=258
    if low > high
        low,high = high,low
    end

    return [Int(clamp(low,1,fq)),Int(clamp(high,1,fq))]
end

function userunit_to_freq(users::Vector{Vector{XY{UserUnit}}})
    [userunit_to_freq(us) for us in users]
end 


function unit_to_time(unit::XY{UserUnit},ctime::Float64,tspan::Float64)
    # use the same
    y1 = unit.y.val
    x1 = ctime + unit.x.val*tspan

    XY{UserUnit}(x1,y1)
end

function unit_to_time(unit::Vector{XY{UserUnit}},ctime::Float64,tspan::Float64)
    
    # use the same values
    y1 = unit[1].y.val
    y2 = unit[2].y.val

    # time = current time + userunit*timespan
    x1 = ctime + unit[1].x.val*tspan
    x2 = ctime + unit[2].x.val*tspan
    
    first = XY{UserUnit}(x1,y1)
    last = XY{UserUnit}(x2,y2)
    return [first,last]
end


function time_to_unit(time::Vector{XY{UserUnit}},ctime::Float64,tspan::Float64)
    
    # use the same values
    y1 = time[1].y.val
    y2 = time[2].y.val
    
    # current position = (time - current time) / timespan
    x1 = (time[1].x.val-ctime)/tspan
    x2 = (time[2].x.val-ctime)/tspan
        
    first = XY{UserUnit}(x1,y1)
    last = XY{UserUnit}(x2,y2)
    return [first,last]
end

function time_to_unit(times::Vector{Vector{XY{UserUnit}}},ctime::Float64,tspan::Float64)
    [time_to_unit(ts,ctime,tspan) for ts in times]
end

# sometimes usv syllabel onset pairs can be mixed, so to make sure they are in time sequence orders
function sequence_check!(ut::Vector{Vector{XY{UserUnit}}}) # Int64 for inds, Float64 for times
    sort!(ut,by = x->(x[1].x.val)) # x[1][1][1]: an usv onset
    nothing
end

"""
# find the indices of each time point
function find_the_time(tp::Float64,ta::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # tp = time point
    # ta = time array e.g. NIDAQ OR RWD OR INTAN
    # ti = time index
    ti = findfirst(x->x>=tp,ta)
    ti,ta[ti] # (index,time)
end

function find_the_time(temp::Vector{Float64,},ta::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    b = zeros(length(temp))
    a = zeros(Int64,length(temp))
    for i in 1:length(temp)
        a[i],b[i] = find_the_time(temp[i],ta)
    end
    return a,b
end"""

# find the indices of each time point
# 2026-01-30 temporary  
function find_the_time(tp::Float64, ta::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
    # tp = time point
    # ta = time array e.g. NIDAQ OR RWD OR INTAN
    # ti = time index
    if !isfinite(tp) || isempty(ta)
        return 1, ta[1]
    end

    #ti = findfirst(x -> x >= tp, ta)
    ti = searchsortedfirst(ta, tp)
    ti = clamp(ti, firstindex(ta), lastindex(ta))
    return ti, ta[ti] # (index,time)
end

function find_the_time(temp::Vector{Float64,}, ta::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
    b = zeros(length(temp))
    a = zeros(Int64, length(temp))
    for i in 1:length(temp)
        a[i], b[i] = find_the_time(temp[i], ta)
    end
    return a, b
end
