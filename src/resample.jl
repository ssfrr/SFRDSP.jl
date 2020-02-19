"""
    resample2(v, rate)

Performs time-domain resampling of `v` by the ratio `rate` using cubic interpolation.
When downsampling the signal is lowpassed before resampling. When up-sampling it is
lowpassed afterwards to remove interpolation noise.

Similarities to `DSP.resample`

- The signals are aligned such that the first sample of `v` corresponds to the first
  sample of `resample2(v, rate)`.

Differences from `DSP.resample`:

- this guarantees that `length(resample2(v, rate)) == round(Int, length(v)*rate)`
- `resample2` uses :reflect_at boundary conditions to avoid ringing when the signal
  starts far from 0.
"""
function resample2(v, rate, lpf=resample_filter2(rate))
    L = length(v)
    Lf = length(lpf)
    @assert isodd(Lf)
    fcomp = Lf÷2
    # reflecting adds Lf-1 samples
    # convolution adds Lf-1 samples
    if rate <= 1
        vf = conv(lpf, pad(v, :reflect_at, fcomp))
        map(0:round(Int, L*rate)-1) do i
            cubic(vf, i/rate+Lf, :reflect_at)
        end
    else
        out = map(0:round(Int, L*rate)-1) do i
            cubic(v, i/rate+1, :reflect_at)
        end
        conv(lpf, pad(out, :reflect_at, fcomp))[Lf:end-Lf+1]
    end
end

function pad(v, bc, pre, post=pre)
    # TODO: support padding longer than the signal, which is useful for large resampling
    # factors
    prepad, postpad = if bc == :reflect_at
        reverse(v[2:pre+1]), reverse(v[end-post:end-1])
    else
        throw(ArgumentError("Unsupported boundary condition type"))
    end

    [prepad; v; postpad]
end

function resample_filter2(rate)
    # DSP.jl doesn't force FIRWindow to be odd-length, so we do that here
    # see https://github.com/JuliaDSP/DSP.jl/issues/345
    cutoff = rate <= 1 ? rate : 1/rate
    ord, α = kaiserord(0.1*cutoff)
    iseven(ord) && (ord += 1)
    digitalfilter(Lowpass(cutoff), FIRWindow(kaiser(ord, α)))
end

# generate data without information too close to nyquist
# x = resample2(rand(100).+10, 2)
# r = 0.5
# plot(0:length(x)-1, x)
# x2 = resample2(x, r)
# plot!((0:length(x2)-1)/r, x2)
#
# x = randn(100000)
# @benchmark resample(x, 0.05)
# @benchmark resample2(x, 0.05)
# currently `resample2` is about 4-5x faster than `resample` for this case
