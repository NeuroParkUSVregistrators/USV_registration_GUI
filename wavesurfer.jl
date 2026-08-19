# functions for wavesurfer reading

# Scaling for the analog singals
function scaledDoubleAnalogData(data::Matrix{Int16},Scales::Matrix{Float64},sc::Matrix{Float64})

    inverseChannelScales=1.0./Scales;
    nScans,nChannels = size(data);
    nCoefficients = size(sc,1) ;
    scaledData=zeros(nScans,nChannels) ;

    if nScans>0 && nChannels>0 && nCoefficients>0
        for j = 1:nChannels
            @inbounds for i = 1:nScans
                datumAsADCCounts = Float64(data[i,j]) ;
                datumAsADCVoltage = sc[nCoefficients,j] ;
                @inbounds for k = (nCoefficients-1):-1:1
                    datumAsADCVoltage = sc[k,j] + datumAsADCCounts*datumAsADCVoltage ;
                end
                scaledData[i,j] = inverseChannelScales[j] * datumAsADCVoltage ;
            end
        end
    end
    return scaledData
end

function matrix_to_vector_vector(data::Matrix{T}) where T
    r,c = size(data)
    data_vec_vec = [data[:,i] for i in 1:c]
end

# load .h5 file from wavesurfer
function loadNIDAQdataFromh5(name::String)
    fid = h5open(name,"r")
    sweep = name[findlast("_",name)[1]+1:findlast(".",name)[1]-1] # check the sweeps

    # matrix
    timestamps = read(fid[string("sweep_",sweep,"/digitalScans")]) # TTL inputs e.g. laer, camera, fiberphotometry
    analogs = read(fid[string("sweep_",sweep,"/analogScans")])
    channelScales = read(fid["header/AIChannelScales"])
    scalingCoefficients = read(fid["header/AIScalingCoefficients"])

    # scale factors
    datas = scaledDoubleAnalogData(analogs,channelScales,scalingCoefficients);

    # convert Matrix{T} to Vector{Vector{T}}
    data = matrix_to_vector_vector(datas)
    ttl = matrix_to_vector_vector(timestamps)

    return data,ttl # Matrix for analog data and TTL inputs
end

# downsampling for analog signals and create a new timevector for it
function analog_down(anal::Vector{Float64},sample::Int64,down::Int64)
    # sample = sampling rate of the original signal e.g. 250khz = 250_000
    # down = the downsampled rate e.g. 1khz = 1000
    rate = down/sample
    anal_resample = resample(anal,rate)
    len = Int(length(anal)*rate)
    time_resample = range(0,step=1/down,length=len)

    return anal_resample,time_resample
end

function analog_down(anal::Vector{Vector{Float64}},sample::Int64,down::Int64)
    rate = down/sample
    new_anal = [resample(a,rate) for a in anal]
    
    len = Int(length(anal[1])*rate)
    time_resample = range(0,step=1/down,length=len)
    
    return new_anal, time_resample
end

# downsampling for ttl signals

function ttl_down(ttl::Vector{UInt8},sample::Int64,down::Int64)

    rate  = Int(sample/down)
    new_ttl_length = Int(length(ttl)/rate)
    new_ttl = zeros(UInt8,new_ttl_length,)

    for n = 1:new_ttl_length
        new_ttl[n] = ttl[rate*(n-1)+1]
    end
    return new_ttl
end

# for V{V{T}}
function td(ttl::Vector{UInt8},rate::Int64,ntl::Int64)
    [ttl[rate*(n-1)+1] for n in 1:ntl]
end

function ttl_down(tl::Vector{Vector{UInt8}},sample::Int64,down::Int64)
    rate  = Int(sample/down)
    new_ttl_length = Int(length(tl[1])/rate)
    new_ttl = [td(t,rate,new_ttl_length) for t in tl]
    
end

# here raw laser ind as UInt8
function digital_only(name::String)
    # load the data
    d,t =  loadNIDAQdataFromh5(name)

    # laser 
    ## TTL downsampling
    tslaser = ttl_down(bit.(t[1],digital_laser),250_000,1000) # 250kHz to 1000Hz

    # camera start time with respect to NIDAQ recording to align laser stim
    tscamera = findfirst(x->x==1,bit.(t[1],digital_camera1))
    camera_start = (tscamera-1)/samplingfrequency

    # time
    time = range(0,step=1/1000,length=length(tslaser))

    # tslaser::Vector{UInt8} time::range  camera_start = Float64
    return tslaser,time,camera_start

end

function find_laser(ls::Vector{UInt8})
    
    ls = Int.(ls)
    diff_ls = diff(ls)
    on_ind = findall(x->x==1,diff_ls) .+1
    off_ind = findall(x->x==-1,diff_ls) 
    ind = [[_a,_b] for (_a,_b) in zip(on_ind,off_ind)]
   return ind
end

# return laser as laser ind 
function digital_only_lsind(name::String)
    # load the data
    d,t =  loadNIDAQdataFromh5(name)

    # laser 
    ## TTL downsampling
    tslaser = ttl_down(bit.(t[1],digital_laser),250_000,1000) # 250kHz to 
    lsind = find_laser(tslaser)

    # camera start time with respect to NIDAQ recording to align laser stim
    tscamera = findfirst(x->x==1,bit.(t[1],digital_camera1))
    camera_start = (tscamera-1)/samplingfrequency

    # time
    time = range(0,step=1/1000,length=length(tslaser))

    # tslaser::Vector{UInt8} time::range  camera_start = Float64
    return lsind,time,camera_start

end
