"""
    circ_smooth(v, L=501)

Smooth a signal with a gaussian window with support over about `L`. Assumes the
signal is periodic.
"""
function circ_smooth(v, L=501)
    N = length(v)
    win = gaussian(N, 0.15*L/N; zerophase=true)
    win ./= sum(win)
    circ_conv(v, win)
end

"""
Perform circular convolution using the FFT algorithm. The result length will be the
greater of the arguments.
"""
function circ_conv(u, v)
    lu = length(u)
    lv = length(v)
    nfft = nextfastfft(max(lu, lv))
    upad = similar(u, nfft)
    copyto!(upad, u)
    for i in lu+1:nfft
        @inbounds upad[i] = zero(eltype(u))
    end
    vpad = similar(v, nfft)
    copyto!(vpad, v)
    for i in lv+1:nfft
        @inbounds vpad[i] = zero(eltype(v))
    end
    p = plan_rfft(upad)
    U = p*upad
    V = p*vpad
    @. U *= V
    # I think the arguments for this are reversed in the FFTW docs
    ldiv!(upad, p, U)
    resize!(upad, max(lu,lv))
end

"""
    circ_resample(v, rate)

Performs time-domain resampling of `v` by the ratio `rate`, treating `v` as
periodic. When downsampling the signal is lowpassed before resampling.
"""
function circ_resample(v, rate)
    L = length(v)
    s = rate < 1 ? circ_smooth(v, 2/rate) : v
    map(0:round(Int, L*rate)-1) do i
        cubic(s, i/rate+1, :periodic)
    end
end
