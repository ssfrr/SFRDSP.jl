const stft_docs = """
The short-time fourier transform can be thought of as a filterbank, with each
filter centered at each fft bin frequency. The window function determines the
frequency response of each bandpass filter (each filter is an FIR filter of the
window function modulated by the center frequency).

This implementation has a few differences from the `stft` implementation in
DSP.jl that I think make it more consistent with the filterbank interpretation,
and make the output easier to reason about:

1. The first window is centered at the first sample of `x`, and in general the
   nth window is centered at `(n-1)*hop`. This also ensures that no input samples
   are multiplied by zero, which would make them unrecoverable.
2. `x` is zero-padded to ensure that the center of the last window is no
   earlier than the last sample of `x`. Similarly this ensures that no information
   from `x` is lost, and also makes the width of the result easier to reason about.
   (it is only a function of the hop size, not the window size).
3. Each window is centered at index 1 before taking the FFT, avoiding a linear
   phase shift in each frame. One benefit of this is that when the same frequency
   component present in `x` is picked up by multiple filters, the phase will be
   coherent across the bands (assuming it is not demodulated).
4. It provides the `demod` option to work with demodulated filter bank outputs.
   This is useful because the phase progresses much more slowly (proportional to
   the difference between the signal frequency and the center frequency, rather
   than at the signal frequency itself).
5. An `istft2` implementation is provided. Note that `istft2` is currently not
   fully compatible with DSP.jl's `stft` function. In order to make it compatible
   we'd need to match the window offset and also perform an `fftshift`.

`stft2` and `istft2` should satisfy the property that for any real-valued signal
`x` and `n`:

    istft2(stft2(x, n), n)[1:length(x)] ≈ x

Note that the round-trip may introduce some zero-padding at the end.
"""

# TODO: support windows longer than nfft
# TODO: add options for boundary conditions so you don't get big bursts
# in the STFT at the beginning and end for truncated signals
"""
    stft2(x, nfft, hop=nfft÷2;
          window=rect,
          demod=false)

Perform a short-time fourier transform on a real-valued time signal `x`. The
result is a complex-valued matrix of dimension:

    (nfft÷2+1, cld(length(x)-1, hop)+1)

- `nfft` is the FFT size.
- `hop` is the hop size from frame to frame. You can think of this as the
   sampling period of the continuous filterbank.
- `window` adds an analysis window. It can be an `AbstractArray`
  of length `nfft`, or a function that can be called as
  `f(nfft; zerophase=true)` to generate such a window. Note that generally the
  window function should be periodic and centered around index `1`, which is
  what the `zerophase` option in DSP.jl's window functions does. See the DSP.jl
  documentation for more details.
- `demod` should be set to `true` if the STFT phase should be demodulated (i.e
  a "sliding window" stft). This is equivilent to de-modulating the output of
  each filter by its center frequency.

$stft_docs
"""
function stft2(x, nfft, hop=nfft÷2;
               window=rect,
               demod=false)
    # we could handle this, but it's simpler not to, and this is the common case
    @assert iseven(nfft)
    @assert hop > 0

    L = length(x)
    N = nfft÷2+1
    # the first window will be centered at the first sample. The last window
    # should be centered at the last sample or later, so that COLA will apply,
    # and we don't need to worry about numerical issues because some samples
    # got scaled way down by the window function.
    M = cld(L-1, hop)+1

    X = zeros(Complex{eltype(x)}, N, M)

    win = if window isa Function
        window(nfft; zerophase=true)
    elseif window isa AbstractArray
        @assert length(window) == nfft
        window
    else
        throw(ArgumentError("window must be an Array or Function"))
    end

    xbuf = zeros(eltype(x), nfft)
    Xbuf = similar(X, N)

    p = plan_rfft(xbuf)

    # each time we copy into the buffer for the FFT, we have to copy the
    # positive part (including the window center) and the negative part
    # separately. We want the window center to be at the beginning of the fft
    # buffer, or else we add a half-buffer circular delay to each frame

    # we could simplify this (at the expense of more copying) by just copying
    # the data in one block and calling `ifftshift`.

    xidx = 1
    nfft2 = nfft÷2

    for c in 1:M
        fill!(xbuf, zero(eltype(xbuf)))
        # copy the positive half to the beginning of the fft buffer
        # note the center of the window could be past the last sample
        # of `x`, in which case `npos` should be 0.
        npos = max(0, min(nfft2, L-xidx+1))
        copyto!(xbuf, 1, x, xidx, npos)
        # copy the negative half to the end of the fft buffer
        nneg = min(nfft2, xidx-1)
        # note that if xidx is past the end of `x`, we don't actually have
        # nneg samples available to copy. it's still useful to compute the
        # indices for the copying though
        copystart = xidx-nneg
        ncopy = max(0, min(nneg, L-copystart+1))
        copyto!(xbuf, nfft-nneg+1, x, copystart, ncopy)
        # zero-out any parts of the buffer that weren't filled. This should
        # only happen at the beginning and end of the process
        # xbuf[npos+1:end-nneg] .= zero(eltype(xbuf))
        if win !== nothing
            xbuf .*= win
        end
        mul!(Xbuf, p, xbuf)
        copyto!(X, (c-1)*N+1, Xbuf, 1, N)
        xidx += hop
    end

    if demod
        # we want to convert to a sliding-window STFT  by multiplying by
        # exp(-jnω). This is equivilent to demodulating each band of the STFT
        # filter bank
        n = (0:M-1)*hop
        ω = range(0, π, length=N)
        # TODO we could do this without allocating a whole copy of X by doing
        # the multipliation frame-wise
        X *= exp.(-im .* n' .* ω)
    end

    X
end


"""
    istft2(X, nfft, hop=nfft÷2;
           out_nfft=nfft,
           out_hop=Int(out_nfft * hop / nfft),
           window=rect,
           demod=false)
Perform an inverse short-time fourier transform on the complex-valued matrix
`X`.

- `nfft` is the FFT size that was used to compute the STFT.
- `hop` is the hop size of the frames in `X`
- `out_nfft` is the desired iFFT size to use for reconstruction.
- `out_hop` is the hop size for the resynthesized windows. You can set this
  to something different from `hop` for basic time stretching (though you will
  likely also need to manually update the phase values).
- `window` adds a resynthesis window. It can be an `AbstractArray`
  of length `nfft`, or a function that can be called as
  `f(nfft; zerophase=true)` to generate such a window. We currently only handle
  the painless case where the round-trip window (elementwise product of the
  analysis and synthesis windows) has the constant-overlap-add property. Note
  that generally the window function should be periodic and centered around
  index `1`, which is what the `zerophase` option in DSP.jl's window functions
  does. See the DSP.jl documentation for more details.
- `demod` should be set to `true` if the STFT phase has been demodulated (i.e
  a "sliding window" stft).

$stft_docs
"""
function istft2(X, nfft, hop=nfft÷2;
                out_nfft=nfft,
                # out_hop will throw an error if the given sizes aren't
                # compatible. TODO: give a better error message or handle this
                # more gracefully
                out_hop=Int(out_nfft * hop / nfft),
                window=rect,
                demod=false)
    # we could handle this case, but it's simpler not to
    @assert iseven(nfft)
    @assert iseven(out_nfft)
    @assert hop > 0
    @assert out_hop > 0

    N = size(X, 1)
    M = size(X, 2)

    # the output signal starts in the center of the first window, and
    # ends no later than the center of the last window
    L = (M-1) * out_hop + 1
    out = zeros(real(eltype(X)), L)

    # TODO: if the element type is e.g. Float32, we should use the same for
    # the window
    # TODO: scaling should probably be handled more automatically, and
    # also we should handle the non-COLA case. Also we're not scaling
    # the windows above because we're assuming they're COLA and add to 1
    # ALSO there are still some artifacts at the very beginning and end,
    # where the samples don't get the right number of windows adding up,
    # so the scale factor is off.
    win = if window isa Function
        window(out_nfft, zerophase=true)
    elseif window isa AbstractArray
        @assert length(window) == out_nfft
        window
    else

    end

    if demod
        # we want to convert to a sliding-signal STFT (which is what most stft
        # implementations do natively) by multiplying by exp(jnω). This is
        # equivilent to modulating the STFT filter bank by the band's center
        # frequency.
        n = (0:M-1)*hop
        ω = range(0, π, length=N)
        # TODO we could do this without allocating a whole copy of X by doing
        # the multipliation frame-wise
        X *= exp.(im .* n' .* ω)
    end

    # this is probably pretty close to supporting nfft != N, but not really
    # tested in that case
    Xbuf = similar(X, out_nfft÷2+1)

    # fill with zeros, so that if we partially fill it later, it ends up being
    # zero-padded
    fill!(Xbuf, zero(eltype(X)))
    outbuf = similar(out, out_nfft)
    p = plan_irfft(Xbuf, out_nfft)
    copylen = min(length(Xbuf), N)

    nfft2 = out_nfft÷2

    outoffset = 0
    for c in 1:M
        copyto!(Xbuf, 1, X, (c-1)*N+1, copylen)
        mul!(outbuf, p, Xbuf)
        outbuf .*= win
        npos = min(nfft2, L-outoffset)
        @simd for i in 1:npos
            @inbounds out[outoffset+i] += outbuf[i]
        end
        nneg = min(nfft2, outoffset)
        @simd for i in 1:nneg
            @inbounds out[outoffset-nneg+i] += outbuf[end-nneg+i]
        end
        outoffset += out_hop
    end

    out
end
