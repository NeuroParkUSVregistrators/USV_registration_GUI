# requring 
# FFTW DSP LinearAlgebra


mutable struct stft_param
    nfft::Int64
    skip::Int64
    nout::Int64
    window::Vector{Float64}
    fin::Vector{Float64}
    tmp::Vector{ComplexF64}
    plan::FFTW.rFFTWPlan{Float64, -1, false, 1, Tuple{Int64}}

    function stft_param(;nfft::Int64=512,skip::Int64=128)
        nout = (nfft >> 1)+1
        window = hanning(nfft)
        fin = zeros(Float64,nfft)
        tmp = zeros(ComplexF64,nout)
        plan = plan_rfft(fin)

        sp = new(nfft,skip,nout,window,fin,tmp,plan)
        return sp
    end
end

mutable struct spec_data
    # param
    sp::stft_param

    # data
    audio::Vector{Float64}
    fout::Matrix{Float64}

    # dependency
    npts::Int64
    fs::Int64
    f::StepRangeLen

    t::StepRangeLen
    t_sng::StepRangeLen

    dt::Float64
    dt_sng::Float64

    function spec_data(data::Vector{Float64},sp::stft_param;fs::Int64=250_000)
    
        fout = Stft(data,sp)
        npts = length(data)
        f = range(0,fs/2,length=sp.nout)
        dt = 1/fs
        t = range(0,step=dt,length=npts) 
        dt_sng = sp.skip/fs
        nblocks = size(fout,2)
        t_sng = range(0,step=dt_sng,length=nblocks)

        new(sp,data,fout,npts,fs,f,t,t_sng,dt,dt_sng)
    end
end

function nblock(npts::Int64,sp::stft_param)
    Int(floor((npts-sp.nfft)/sp.skip))+1  
end

# apply hanning window filter to smoothe out the ends
function sig!(data::Vector{Float64},f::Vector{Float64},ui::UnitRange{Int64},wa::Vector{Float64})
    # data = whole audio 1 to end
    # f = filtered data mutated
    # ui = unit range 1:512 513:1024 ....
    # wa = hanning window
    of =1
    for i = ui
        @inbounds f[of] = data[i]*wa[of]
        of+=1
    end
    nothing
end

# convert into real number
function fout!(fout::Matrix{Float64},tmp::Vector{ComplexF64},offset::Int64)
    for i = 1:length(tmp)
        @inbounds fout[offset+i] = abs2(tmp[i])
    end
    nothing
end

# low level one
function Stft(data::Vector{Float64},
                        nfft::Int64,skip::Int64,nout::Int64,
                        window::Vector{Float64},
                        fin::Vector{Float64},tmp::Vector{ComplexF64},
                        plan::FFTW.rFFTWPlan{Float64, -1, false, 1, Tuple{Int64}})
    npts = length(data) # for spectrogram not fixed length inputs
    nblocks = Int(floor((npts-nfft)/skip))+1
    fout = zeros(Float64,nout,nblocks)
    offset = 0
    @inbounds begin for k = 1:nblocks
            period = (1:nfft) .+ (k-1)*skip
            sig!(data,fin,period,window)
            mul!(tmp,plan,fin)
            fout!(fout,tmp,offset)
            offset += nout
        end
    end
    return fout
end
# with stft_param
function Stft(data::Vector{Float64},sp::stft_param)
    Stft(data,sp.nfft,sp.skip,sp.nout,sp.window,sp.fin,sp.tmp,sp.plan)
end

# for @view version
function spectral(sng::SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false})
    sbl_sum = sum(sng,dims=1)
    sbl_max = maximum(sng,dims=1)
    
    sp = sbl_max./sbl_sum
    nanind = findall(x->x==1,isnan.(sp))
    sp[nanind] .= 0
    
    return sp, sbl_sum
end

function spectral(sng::Matrix{Float64})
    sbl_sum = sum(sng,dims=1)
    sbl_max = maximum(sng,dims=1)
    
    sp = sbl_max./sbl_sum
    nanind = findall(x->x==1,isnan.(sp))
    sp[nanind] .= 0
    
    return sp, sbl_sum
end

function frequency(sng::SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false},
                    sp::spec_data,lowbound::Int64)
        
    n = argmax(sng,dims=1)    
    freq = zeros(Float64,size(sng,2))
    
    nran = length(n)
    for j in 1:nran
        freq[j]= sp.f[lowbound+ n[j][1]-1]
    end
    
    meanfreq = mean(freq)
    pitch = (std(freq))^2
    if pitch == 0 || pitch == NaN
        println("something wrong: may want to press delete button if no orange box is visible")
        pitch = 100_000
    end
    outpitch = log10(pitch)
    # double check whether pitch is NaN
    if isnan(outpitch)
        outpitch = 5.0
    end
    return meanfreq,outpitch
end

function frequency(sng::Matrix{Float64},
                    sp::spec_data,lowbound::Int64)
        
    n = argmax(sng,dims=1)    
    freq = zeros(Float64,size(sng,2))
    
    nran = length(n)
    for j in 1:nran
        freq[j]= sp.f[lowbound+ n[j][1]-1]
    end
    
    meanfreq = mean(freq)
    pitch = (std(freq))^2
    if pitch == 0
        println("something wrong")
        pitch = 100_000
    end
    return meanfreq,log10(pitch)
end



function power_to_decibels(x::Float64,x_ref=10e-12)
    10log10(x/x_ref)
end
