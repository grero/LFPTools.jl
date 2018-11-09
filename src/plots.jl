using PyPlot
using Colors
using DSP

function plot_lfp(Xβ::Matrix{Float64}, Xγ::Matrix{Float64}, x;labels=["0.1 - 10 Hz", "20 - 40 Hz"])
    fig = plt[:figure]()
    ax1 = fig[:add_subplot](211)
    ax2 = fig[:add_subplot](212)
    ntrials = size(Xγ,2)
    _cc = parse(Colorant, "#1f77b4")
    for (ax,X) in zip([ax1, ax2], [Xβ, Xγ])
        μ = mean(X, 2)[:]
        σ = std(X,2)[:]
        se = σ/sqrt(ntrials)
        ax[:plot](x, μ)
        ax[:fill_between](x, μ-2*se, μ+2*se;color=(_cc.r, _cc.g, _cc.b, 0.8))
        ax[:axvline](0.0; color="r", label="Fixation start")
        ax[:axvline](0.35; color="k", label="Target on")
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
    end
    ax1[:set_xticklabels]([])
    ax1[:set_title](labels[1])
    ax2[:set_title](labels[2])
    plt[:legend]()
    fig
end

function plot_spectrogram(Y::Union{Array{T}, DSP.Periodograms.Spectrogram}, args...;kvs...) where T <: Real
    fig = plt[:figure]()
    ax = fig[:add_subplot](111)
    plot_spectrogram(ax, Y, args...;kvs...)
    fig
end

"""
Plot the time-frequency power spectrum of the signal `Y` sampled at sampling rate `fs` using the matplotlib axes object `ax`. The parameter `t0` indicates the time in seconds of the first point in `Y`. The optional parameeter `logscale` indicates wheter to plot the log of the power, and `interpolation` specifies an interpolation method (defaults to "sinc") to be applied to the image.
"""
function plot_spectrogram(ax, Y::Vector{T}, t0::Real, fs=1000;kvs...) where T <: Real
    PP =  spectrogram(Y, fs=fs)
    plot_spectrogram(ax, PP, t0, fs;kvs...)
end

function plot_spectrogram(ax, Y::Array{T,2}, t0::Real, fs=1000;kvs...) where T <: Real
    nbins, ntrials = size(Y)
    PP = spectrogram(Y[:,1], fs=fs)
    μ = PP.power
    σ² = PP.power.*PP.power
    for ti in 2:ntrials
        PP = spectrogram(Y[:,ti], fs=fs)
        pp = PP.power
        μ .+= pp
        σ² .+= pp.*pp
    end
    μ ./= ntrials
    σ² .= sqrt.(σ²/ntrials - μ.*μ)
    μ, σ², PP.freq, PP.time
end

function plot_spectrogram(ax, PP::DSP.Periodograms.Spectrogram, t0::Real, fs=1000;fmax=fs/2, logscale=false, interpolation="sinc",events=Float64[])
    fidx = find(PP.freq .<= fmax)
    if logscale
        ppower = log10.(PP.power)
        clabel = "log Power"
    else
        ppower = PP.power
        clabel= "Power"
    end
    II = ax[:imshow](ppower[fidx,:];aspect="auto", extent=(t0, t0 + PP.time[end], PP.freq[1], PP.freq[fidx[end]]), origin="lower",interpolation=interpolation)
    for ee in events
        ax[:axvline](ee;color="w")
    end
    if length(ax[:figure][:axes])>1
        plt[:colorbar](II,cax=ax[:figure][:axes][2];label=clabel)
    else
        plt[:colorbar](II;label=clabel)
    end
    ax[:set_xlabel]("Time [s]")
    ax[:set_ylabel]("Frequency [Hz]")
    ax[:spines]["right"][:set_visible](false)
    ax[:spines]["top"][:set_visible](false)
end
