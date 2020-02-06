"""
    fit_ir(x::AbstractVector, fs, nfft, max_decay)
    fit_ir(x::SampleBuf, nfft, max_decay)

Analyze the given impulse response `x` at sample rate `fs` by curve fitting in
each frequency band. The impulse response audio should be relatively clean (not
too much nonstationary background noise).

`nfft` gives the number of frequency bands that will be used.

`max_decay` should be the longest time (in seconds) that it might take for
the impulse to decay to the noise floor. Setting this too high may affect
convergence, and setting it too low may effect accuracy.

Returns a vector of named tuples for each band containing:
    - `freq` - the band center frequency (in Hz)
    - `noise_pow` - noise floor (in dB)
    - `peak_snr` - Peak SNR (in dB)
    - `rt60` - RT60 for each band (in seconds)
    - `drr` - Direct-to-Reverberant Ratio for each band (in dB)

If the curve fit does not converge the RT60 and DRR may be `nothing`.
"""
fit_ir(x::SampleBuf, nfft, max_decay) =
    fit_ir(x.data, samplerate(x), nfft, max_decay)

function fit_ir(x::AbstractVector, fs, nfft, max_decay)
    maxidx = findmax(abs2.(x))[2]
    # now we take the STFT to generate our frequency bands
    hop = nfft÷2
    if maxidx <= 2hop
        throw(ErrorException("IR Fitting needs at least $(2hop) samples (2 hops) prior to peak (at $maxidx). Provide more signal or decrease your hop size"))
    end
    # set things up so that the 2nd STFT frame is centered on the peak.
    # this also ensures that no padding is present in that window
    spec = stft2(x[maxidx-2hop:end], nfft, hop; window=hanning)
    specpow = abs2.(spec)
    logspecpow = 10log10.(specpow)

    # TODO: we could probably re-structure this to also allow for nmax
    # starting too low, but then it seems less guaranteed that it will
    # converge
    nmax_init = round(Int, max_decay*fs/hop)+npeak

    noise = 10log10.(mean(specpow[:, nmax_init:end], dims=2))

    map(1:size(specpow, 1)) do k
        # The peak index should be 2, but in practice it seems delayed for some
        # frequencies. All the following analysis is performed relative to the
        # peak position
        npeak = findmax(logspecpow[k, :])[2]
        if npeak == 1
            throw(ErrorException("Bin $k has a peak that's too early, this function should be modified to add more space before the time-domain peak"))
        end
        nmax = nmax_init
        # from Traer they say the decay starts after about 20ms. We want to make
        # sure we're in the exponential decay portion, so we wayt 50ms  after
        # the peak:
        # nstart = ceil(Int, 0.05 * fs / hop) + npeak
        # unfortunately this process somewhat sensitive to how much of the early
        # IR we include - if I try to skip the first part as in the
        # commented-out code above then it seems to give much more shallow fits,
        # which correspond to much larger RT60s. It's not super clear from Traer
        # whether he included the early IR in the fit.
        nstart = npeak

        # let's try a linear fit from nstart to nmax. We'll then use that to
        # figure out where the decay hits the noise floor. Initially the fitted
        # data will include the noise that's after the decay, so we'll
        # iteratively shorten the region we're fitting over until it converged.
        converged = false
        local L

        while true
            a = [nstart:nmax ones(nmax-nstart+1)]
            b = logspecpow[k, nstart:nmax]
            # aL = b, so a \ b = L
            L = a \ b
            # positive slope - not going to converge
            L[1] >= 0 && break
            # starting below noise floor - not going to converge
            L[2] < noise[k] && break
            nmax_prev = nmax
            nmax = round(Int, (noise[k] - L[2]) / L[1])
            nmax >= nmax_init && break
            # fit region is too small to be reliable
            nmax <= nstart + 10 && break
            noise[k] = 10log10.(mean(specpow[k, nmax:end]))
            # I think it's possible that this could get into an oscillation,
            # so we'll consider any increase in the noise intersect point to be
            # convergence
            if nmax >= nmax_prev
                converged = true
                break
            end
        end

        # we'll estimate the peak power with a quadratic interpolation. The
        # point we interpolate at is based on the relative magnitudes of the
        # neighboring points.
        # details at:
        # https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html

        peakpow, peakpos = let
            α, β, γ = specpow[k, (-1:1).+npeak]
            p = 0.5 * (α-γ)/(α-2β+γ)
            pow = β - 0.25 * (α-γ)*p
            (10log10(pow), p+npeak)
        end

         # the peak is in the 2nd STFT frame, so that's the part of the linear
         # fit that we compare to for the peak SNR
        (freq=(k-1)/nfft*fs,
         noise_pow=noise[k], # noise power
         peak_snr=peakpow - noise[k], # peak SNR
         drr=converged ? peakpow - (L[2] + peakpos*L[1]) :
                         nothing, # DRR
         rt60=converged ? -60/L[1]*hop/fs : nothing) # RT60
    end
end
