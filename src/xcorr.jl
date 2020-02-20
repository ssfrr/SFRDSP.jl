"""
    circ_xcorr2(u,v; norm=true, center=true, phat=false)

Perform circular cross-correlation using the FFT algorithm, shifting `v`
relative to `u`, so delaying `u` will shift the result to the right. Index 1
corresponds to zero lag.

Another interpretation is that index of the peak of this cross correlation
(minus 1) represents the delay of `u` with respect to `v`, or alternatively the
position of `v` within `u`.

If `phat=true`, this implements GCC-PHAT, where the amplitudes are set to one
and only the phases are used.
"""
function circ_xcorr2(u,v; norm=true, center=true, phat=false)
    L = length(u)
    @assert length(v) == L
    if center
        u .-= mean(u)
        v .-= mean(v)
    end
    p = plan_rfft(u)
    U = p*u
    V = p*v
    if phat
        @. U = U * conj(V) / (abs(U) * abs(V))
    else
        @. U *= conj(V)
    end
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
           center=true,
           phat=false,
           plantype=FFTW.ESTIMATE)

Perform linear cross-correlation using the FFT algorithm, shifting `v` relative
to `u`, so delaying `u` will shift the result to the right. Index `1` corresponds to
zero lag, with negative lag wrapped to the end of the result vector.

If `unbiased` is true (the default) then this will compensate for lost energy
due to lags where the signals don't entirely overlap. This is generally appropriate
when the target signals being correlated occupy a substantial fraction of `u` and `v`.

If `phat=true`, this implements GCC-PHAT, where the amplitudes are set to one
and only the phases are used.

`plantype` can be set to one of:
    - `FFTW.ESTIMATE` (the default)
    - `FFTW.MEASURE`
    - `FFTW.PATIENT`
    - `FFTW.EXHAUSTIVE`

See the FFTW documentation for more information on these settings. The plan information
is saved internally within FFTW on each session, so you only pay the cost of planning
the first time `xcorr2` is called.

!!! WARNING: Currently (Feb 2020) I get very long pauses and sometimes hangs with any
plan other than `FFTW.ESTIMATE`
"""
function xcorr2(u,v;
                lagbounds=(-length(v)+1, length(u)-1),
                unbiased=true,
                norm=true,
                center=true,
                phat=false,
                plantype=FFTW.ESTIMATE)
    # TODO make sure it works for zero-length vectors
    # TODO: handle minlag == maxlag
    # TODO: really `unbiased` is a form of normalization, where we normalize
    # each lag by the number of positions that overlap. This also suggests yet
    # another normalization approach where we normalize each lag by geometric
    # mean (similar to the current `norm=true` but we only include the parts
    # that overlap)
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
    nfft = nextfastfft(max(su - minlag, sv + maxlag))
    upad = zeros(float(eltype(u)), nfft)
    vpad = zeros(float(eltype(v)), nfft)

    copyto!(upad, 1, u, 1, su)
    copyto!(vpad, 1, v, 1, sv)

    if center
        upad[1:su] .-= mean(u)
        vpad[1:sv] .-= mean(v)
    end
    # we don't normalize if we're using GCC-PHAT because we set all the amplitudes
    # to 1
    normval = if norm && !phat
        sqrt(sum(abs2, upad)*sum(abs2, vpad))
    else
        one(eltype(upad))
    end

    # if unbiased
    #     minlag_norm =
    #     maxlag_norm =
    #     biasnorm = map(minlag:maxlag) do lag
    #     end
    # end

    p = plan_rfft(upad; flags=plantype | FFTW.DESTROY_INPUT)
    uspec = p * upad
    vspec = p * vpad
    magthresh = sqrt(eps())
    for i in eachindex(uspec)
        @inbounds uspec[i] *= conj(vspec[i])
        if phat
            @inbounds mag = sqrt(abs2(uspec[i]) * abs2(vspec[i]))
            mag > magthresh && (uspec[i] /= norm)
        end
    end

    result = inv(p) * uspec
    # remove any padding and aliased data from the middle of the buffer
    # for lagbounds=(-3,4) the result would look like:
    # [0PPPPNNN] (0 is 0-lag, P for positive, N for negative)
    # TODO: support fully positive and fully negative lagbounds
    @assert minlag <= 0
    @assert maxlag >= 0
    @assert length(result) == nfft
    copyto!(result, maxlag+2, result, nfft+minlag+1, -minlag)
    resize!(result, nlags)

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
                # this power is empirical - I haven't throught through whether
                # there's a theoretical reason, or whether it is signal-dependent.
                # this is probably another good reason to normalize each lag based
                # on the actual energy in the overlap instead of this estimate
                result[i] *= (maxoverlap / overlap)^0.5
            end
        end
    end

    if norm
        result ./= normval
    end
    result
end
