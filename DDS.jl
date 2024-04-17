#Функция генерирует сигнал методом DDS используя таблицу сигнала

using FFTW
using Plots
using DSP

plotly()
SYS_CLK=90E6
F0 = 30E3
t = 0:1/SYS_CLK:1/F0*120

function sinus_dds()
    tbl=UInt16[]
    for i in 1:256
        s = round(UInt16,(sin(2*pi*(i-1)/256)+1)*(4095+1)/2,RoundDown)
        push!(tbl,s)
    end
    return tbl
end
sinus = sinus_dds()

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
    phase_accum = UInt16(0x0001)
    incr = round(UInt16,2^16*F0/SYS_CLK,RoundUp)
    tb_val = tb
    signal = typeof(tb[1])[]
    for i in t
        addr=phase_accum>>8
        addr=addr+0x001
        if addr == 0x0000
            addr = 0x0001
        end
        signal = vcat(signal,tb_val[addr])
        phase_accum = phase_accum + incr
    end
    # println(incr)
    return signal
end

function DDS_PSK2(SYS_CLK,F0,t,tb,sequence)
    phase_accum = UInt16(0)
    incr = UInt16(0)
    incr = UInt16(round(2^16*F0/SYS_CLK))
    tb_val = tb
    signal = typeof(tb_val[1])[]
    seq=repeat(sequence,inner=round(Int,length(t)/length(sequence)))
    seq=seq[1:(length(seq)-(length(seq)-length(t)))]
    flag=false;

    for i in seq
        if i == 1
            addr=phase_accum>>8
            if addr == 0x0000
                addr = 0x0001
            end
            signal = vcat(signal,tb_val[addr])
            phase_accum = phase_accum + incr
            flag = true
        else

            if flag == true
                phase_accum = phase_accum + round(UInt16,typemax(typeof(phase_accum))/2)
                flag == false
            else
                addr=phase_accum>>8
                # if addr == 0x0000
                #     addr = 0x0001
                # end
                signal = vcat(signal,tb_val[addr])
                phase_accum = phase_accum + incr
            end

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
seq=UInt16[1,0,1,1]

saw=DDS_saw(SYS_CLK,F0,t)
saw=(saw*3.3/typemax(UInt8)).-3.3/2
sig=DDS(SYS_CLK,F0,t,sinus)
sig = (sig*3.3/4096).-3.3/2
# meandr=DDS(SYS_CLK,F0,t,pls)

#y=fft(sig)
#G=2*abs.(y)/length(y)
#plot(sig)
#plot(10*log.(10,G))
##
response = Lowpass(F0, fs=SYS_CLK)
designmethod = FIRWindow(hamming(length(sig)))
LowFilter = digitalfilter(response, designmethod)

sig=filt(sig,LowFilter)
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

signal=3.3/2*sin.(2*pi*F0*t)
plot(sig)
plot!(signal)

# Y1=fft(signal)
# Y2=fft(sig)
# G1=abs.(Y1)
# G2=abs.(Y2)
# plot(10*log.(G1))
# plot!(10*log.(G2))