CLK=90E6
F0 = 30E3

t=0:1/CLK:10*1/F0
using Plots
plotly()

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
plot(t[1:ind],sin.(2*pi*F0*t[1:ind]))

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

sinus=generate_sinus(ind,t,F0)
plot(sinus)