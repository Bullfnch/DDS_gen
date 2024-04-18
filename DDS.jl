#Функция генерирует сигнал методом DDS используя таблицу сигнала

using FFTW
using Plots
using DSP

plotly()
SYS_CLK=90E6
F0 = 30E3
T=1/F0
t = 0:1/SYS_CLK:T

function sinus_dds()
    tbl=UInt16[]
    for i in 1:256
        s = round(UInt16,(sin(2*pi*(i-1)/256)+1)*(4095+1)/2,RoundDown)
        push!(tbl,s)
    end
    return tbl
end
sinus = sinus_dds()
sinus_inv = sinus.*(-1).+4096

function sinus_index(time,f0)
    stps=0
    index_zero=0
    for i in 2:1:length(time)
        sin_val=sin(2*pi*time[i]*f0)
        if  abs(sin_val) < 2E-15
            if stps==1
                index_zero=i-1
                break
            end
            stps=stps+1
        end
    end
    return index_zero
end

ind=sinus_index(t,F0)
# plot(t[1:ind],sin.(2*pi*F0*t[1:ind]))

function generate_sinus(index,time,f0)
    sinus=Float64[]
    N=length(time)
    for i in 1:N
        if i>index && i!=N
            cycle=round(Int,i/index,RoundDown)*index
            if cycle == i
                cycle=i-1
            end
            sample=sin(2*pi*f0*time[i-cycle])
            push!(sinus,sample)
        else
            sample=sin(2*pi*f0*time[i])
            push!(sinus,sample)
        end
        
    end
    return sinus
end

# sinus_timed=generate_sinus(ind,t,F0)


function pulse()
    tbl=UInt16[]
    for i in 1:256
        if i>127
            s=0
            push!(tbl,s)
        else
            s=255
            push!(tbl,s)
        end
    end
    return tbl
end

pls=pulse()


function DDS(SYS_CLK,F0,t,tb)

    phase_accum = UInt32(0x0001)
    incr = round(UInt32,2^32*F0/SYS_CLK)
    println(incr)
    tb_val = tb
    signal = typeof(tb[1])[]
    for i in t
        addr=phase_accum>>29
        addr=addr+0x1
        if addr>0x006
            phase_accum=0x0001
            addr=0x001
        end
        # println(addr)
        signal = vcat(signal,tb_val[addr])
        phase_accum = phase_accum + incr
    end
    return signal

end

seq=[1,0,1,0,0,1]

cod=DDS(SYS_CLK,F0,t,seq)

function DDS_PSK2(SYS_CLK,F0,t,tb1,tb2,sequence)
    phase_accum = UInt32(0x0001)
    incr = round(UInt32,2^32*F0/SYS_CLK)
    tb_val_1 = tb1
    tb_val_2 = tb2
    tb_val = tb1
    signal = typeof(tb_val[1])[]
    seq=typeof(sequence[1])[]
    stps=0
    ind=1
    while length(seq)<length(t)
        if stps == round(Int,SYS_CLK/F0*2)
            stps = 0
            ind = ind+1
            if ind>length(sequence)
                ind=1
            end
            push!(seq,sequence[ind])
        else
            if ind>length(sequence)
                ind=1
            end
            push!(seq,sequence[ind])
            stps=stps+1
        end
    end
    for i in 1:length(t)-1
        
        if (seq[i+1]-seq[i])!=0 && seq[i+1]==1
            addr=phase_accum>>24
            tb_val=tb_val_1
            signal = vcat(signal,tb_val[addr+0x0001])
            phase_accum = phase_accum + incr
        elseif  (seq[i+1]-seq[i])!=0 && seq[i+1]==0
            addr=phase_accum>>24
            tb_val=tb_val_2
            signal = vcat(signal,tb_val[addr+0x0001])
            phase_accum = phase_accum + incr
        else
            addr=phase_accum>>24
            signal = vcat(signal,tb_val[addr+0x0001])
            phase_accum = phase_accum + incr
        end
    end
    return signal
end


function DDS_saw(SYS_CLK,F0,t)
    phase_accum = UInt16(0)
    incr = UInt16(0)
    incr = UInt16(round(2^16*F0/SYS_CLK))
    signal=typeof(phase_accum)[]
    for i in t
        signal = vcat(signal,phase_accum>>8)
        phase_accum = phase_accum + incr
    end
    return signal
end


# Блок генерации сигналов

# bpsk=DDS_PSK2(SYS_CLK,F0,t,sinus,sinus_inv,[1,0,1,0,0,1])

# saw=DDS_saw(SYS_CLK,F0,t)
# saw=(saw*3.3/typemax(UInt8)).-3.3/2
# sig=DDS(SYS_CLK,F0,t,sinus)
# sig = (sig*3.3/4096).-3.3/2

# meandr=DDS(SYS_CLK,F0,t,pls)

# mnd=[1,0]
# mnd=repeat(mnd,round(Int,t[end]/T))
# mnd=repeat(mnd,inner=round(Int,length(meandr)/length(mnd)))

#y=fft(sig)
#G=2*abs.(y)/length(y)
#plot(sig)
#plot(10*log.(10,G))
##
# response = Lowpass(F0, fs=SYS_CLK)
# designmethod = FIRWindow(hamming(length(sig)))
# LowFilter = digitalfilter(response, designmethod)

# sig=filt(sig,LowFilter)
function upsample(x,order)
    ind_x=1:length(x)
    ind_y=1:length(x)*order
    ind_y_x=1:order:length(x)*order
    y=typeof(x[1])[]
    id=1
    stps=0
    push!(y,x[id])
    while length(y)<length(x)*order
        if stps == order-1
        push!(y,x[id])
        stps = 0
        id=id+1
        else
            push!(y,0)
            stps = stps+1
        end
    end
    return y
end

# saw=upsample(saw,10)

# signal=3.3/2*sin.(2*pi*F0*t)

# plot(sig)

# plot(sig)
# plot!(3.3/2*sinus_timed)

# plot(meandr/255)
# plot!(mnd)

# Y1=fft(signal)
# Y2=fft(sig)
# G1=abs.(Y1)/length(Y1)
# G2=abs.(Y2)/length(Y2)
# plot(10*log.(G1))
# plot!(10*log.(G2))

# plot(t[1:end-1],bpsk/4096)