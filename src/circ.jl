"""
    circ_smooth(v, L=501)

Smooth a signal with a gaussian window with support over about `L`. Assumes the
signal is periodic.
"""
function circ_smooth(v, L=501)
    N = length(v)
    win = gaussian(N, 0.15*L/N; zerophase=true)
    win ./= sum(win)
    irfft(rfft(v) .* rfft(win), N)
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
