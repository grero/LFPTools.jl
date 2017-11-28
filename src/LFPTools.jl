module LFPTools
using DSP
using LsqFit


sinfunc(x, f0, p) = p[1] + p[2]*sin.(2*pi*f0*x+p[3])
function sinfunc!(y, x, f0, p)
    for i in 1:length(y)
        y[i] = sinfunc(x[i], f0, p)
    end
    y
end

"""
Estimate the parameters of a sinusoid with frequency `f0` by band-pass filtering the signal in `Y` between `f0 - 30.0 and f0 + 30.0` Hz.
"""
function estimate_sinusoid(Y::AbstractVector{T}, f0=50.0,fs=30_000.0, x=(0.0:1/fs:(length(Y)-1)/fs)) where T <: Real
    ff2 = digitalfilter(Bandpass(f0-30.0, f0+30.0;fs=fs), Butterworth(4))
    Y2 = filtfilt(ff2, Y)
    YY = zeros(Y2)
    func(x,p) = sinfunc!(YY, x, f0, p)
    pp = curve_fit(func, x, Y2, [0.0, 1.0, 0.0])
    pp.param
end

"""
Remove line noise at frequency `f0` Hz from the signal `Y`.
"""
function remove_linenoise(Y, f0=50.0, fs=30_000.0, x=(0.0:1/fs:(length(Y)-1)/fs))
    pp = estimate_sinusoid(Y, f0, fs, x)
    Y2 = copy(Y)
    remove_linenoise!(Y2, pp, f0, fs)
    Y2, pp
end

"""
Remove line noise at frequency `f0` in-place using the fitted parameters `pp`.
"""
function remove_linenoise!(Y, pp, f0=50.0, fs=30_000.0)
    dt = 1/fs
    x = 0.0
    for i in 1:length(Y)
       Y[i] -= sinfunc(x,f0, pp)
       x += dt
    end
end

end#module
