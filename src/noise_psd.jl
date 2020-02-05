"""
    noise_psd(x::Vector{<:Real}, smoothing, samplerate, nfft, hop=nfft÷2)
    noise_psd(X::Matrix{<:Complex}, smoothing, framerate)

Estimate the noise power-spectral-density of a noisy time-domain signal `x`,
or an STFT `X` using minimum statistics. This is similar to the approach in [1], but assumes a constant noise PSD. To reduce variance the signal will be
smoothed with a gaussian window `smoothing` seconds long. Genrally longer windows will provide a less-biased estimate, but it must be short enough that there will be some regions of the signal where the the smoothing window has none of the target signal in it.

`samplerate` is used for time-domain signals, vs `framerate` for STFTs. You can
convert between them with `framerate = samplerate / hop`.

The result will be of size `nfft÷2+1` for a time signal `x`, or the same height
as an STFT `X`.


[1] Martin, R., Spectral subtraction based on minimum statistics, 1994
"""
function noise_psd(x::AbstractVector{<:Real},
                   smoothing, samplerate,
                   nfft, hop=nfft÷2)
    noise_psd(stft2(x, nfft, hop), smoothing, samplerate/hop)
end

function noise_psd(X::AbstractMatrix{<:Complex}, smoothing, framerate)
    L = round(Int, smoothing * framerate)
    if iseven(L)
        L += 1
    end
    win = gaussian(L, 0.14)
    win ./= sum(x->x^2, win)
    map(1:size(X, 1)) do band
        sm = conv(abs2.(X[band, :]), win)
        minimum(sm[L:end-L+1])
    end
end
