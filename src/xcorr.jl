"""
    circ_xcorr2(u,v; norm=true, center=true)

Perform circular cross-correlation using the FFT algorithm, shifting `v`
relative to `u`, so delaying `u` will shift the result to the right. Index 1
corresponds to zero lag.

Another interpretation is that index of the peak of this cross correlation
(minus 1) represents the delay of `u` with respect to `v`, or alternatively the
position of `v` within `u`.
"""
function circ_xcorr2(u,v; norm=true, center=true)
    L = length(u)
    @assert length(v) == L
    if center
        u .-= mean(u)
        v .-= mean(v)
    end
    p = plan_rfft(u)
    U = p*u
    V = p*v
    U .*= conj.(V)
    result = irfft(U, L)

    if norm
        result ./= sqrt(sum(x->x^2, u)*sum(x->x^2, v))
    end

    result
end

"""
    xcorr2(u,v;
           lagbounds=(-length(v)+1, length(u)-1),
           unbiased=true,
           norm=true,
           center=true)

Perform linear cross-correlation using the FFT algorithm, shifting `v` relative
to `u`, so delaying `u` will shift the result to the right. Index
`-lagbounds[1]+1` (default `length(v)`) corresponds to zero lag.

If `unbiased` is true (the default) then this will compensate for lost energy
due to lags where the signals don't entirely overlap.
"""
function xcorr2(u,v;
                lagbounds=(-length(v)+1, length(u)-1),
                unbiased=true,
                norm=true,
                center=true)
    # TODO make sure it works for zero-length vectors
    # TODO: handle minlag == maxlag
    # TODO: really `unbiased` is a form of normalization, where we normalize
    # each lag by the number of positions that overlap. This also suggests yet
    # another normalization approach where we normalize each lag by geometric
    # mean (similar to the current `norm=true` but we only include the parts
    # that overlap)
    if center
        u .-= mean(u)
        v .-= mean(v)
    end
    su = length(u)
    sv = length(v)
    minlag, maxlag = lagbounds
    minlag <= maxlag || throw(ArgumentError(
        "maximum lag must be greater or equal to minimum lag"))
    # if `u` and `v` are different lengths and the lags don't cover the whole
    # signal then we might not need all of the signals. E.g. if `u` is longer
    # than `v` by more than `minlag` then the end of the `u` is never used,
    # so we can chop it off.
    su, sv = if su > sv
        clamp(sv + maxlag, sv, su), sv
    elseif su < sv
        su, clamp(su - minlag, su, sv)
    else
        su, sv
    end
    nlags = maxlag - minlag + 1
    # TODO: we should use `nextfastfft`
    nfft = max(su - minlag, sv + maxlag)
    upad = zeros(eltype(u), nfft)
    vpad = zeros(eltype(v), nfft)

    copyto!(upad, 1, u, 1, su)
    # pre-shift `v` to correspond to the `minlag` position
    if minlag < 0
        copyto!(vpad, 1, v, -minlag+1, sv+minlag)
        copyto!(vpad, nfft+minlag+1, v, 1, -minlag)
    else
        copyto!(vpad, minlag+1, v, 1, sv)
    end

    p = plan_rfft(upad)
    uspec = p * upad
    vspec = p * vpad
    uspec .*= conj.(vspec)

    result = resize!(irfft(uspec, nfft), nlags)

    if unbiased
        # maximum overlap when lag == 0
        maxoverlap = min(su, sv)
        for (i, lag) in enumerate(UnitRange(lagbounds...))
            overlap = if lag < 0
                max(0, maxoverlap+lag)
            elseif lag < su - sv
                maxoverlap
            else
                max(0, maxoverlap - lag + su - sv)
            end
            if overlap > 0
                result[i] *= maxoverlap / overlap
            end
        end
    end

    if norm
        # TODO: this could be a little more efficient by combining with the
        # unbiasing and/or only applying it to the shortest of u, v, and result
        result ./= sqrt(sum(x->x^2, u)*sum(x->x^2, v))
    end
    result
end
