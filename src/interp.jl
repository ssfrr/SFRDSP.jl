"""
    cubic(v, i, bc=:constant, constval=zero(eltype(v)))

Evaluate the vector `v` at the interpolated index `i` using cubic interpolation.
`bc` sets the boundary condition, which is used for extrapolation of the signal,
when `i` is outside the original support of `v`.

`bc` can be set to one of:
- :constant     - The function is `constval` outside the given indices
- :nearest      - The first and last values are extended
- :linear       - The first two and last two values are used to compute a slope
- :periodic     - The function is assumed to repeat infinitely
- :reflect_at   - The function is assumed to reflect at the first and last
                  points (so they don't repeat). This is equivalent to scipy's
                  `mirror` mode.
- :reflect_past - The function is assumed to reflect at a point 1/2 sample past
                  the first and last points, so they repeat. This is equivalent
                  to scipy's `reflect` mode.
"""
function cubic(v, i, bc=:linear, constval=zero(eltype(v)))
    # TODO: add :quadratic and :cubic boundary conditions that set a constant
    # 2nd and 3rd derivatives
    # TODO: for efficiency we also may want a `cubic_inbounds` that doesn't
    # check `i` and assumes we're within the signal. Definitely set up some
    # profiling and comparison to existing interpolation first.
    L = length(v)
    idxs = floor(Int, i) .+ (-1, 0, 1, 2)
    # generate y valuels, extrapolating if we're near the edge of the signal
    ys = if bc == :constant
        map(idx->1 <= idx <= L ? v[idx] : constval, idxs)
    elseif bc == :nearest
        map(x->v[clamp(x, 1, L)], idxs)
    elseif bc == :linear
        map(idxs) do idx
            if idx < 1
                v[1] + (v[2]-v[1]) * (idx-1)
            elseif idx > L
                v[L] + (v[L]-v[L-1]) * (idx-L)
            else
                v[idx]
            end
        end
    elseif bc == :periodic
        map(x->v[mod1(x, L)], idxs)
    elseif bc == :reflect_at # reflect at the edge
        # ABCDCB|ABCD|CBABCD
        map(idxs) do idx
            idx = mod1(idx, 2L-2)
            idx > L ? v[2L-idx] : v[idx]
        end
    elseif bc == :reflect_past # reflect 1/2 sample past the edge
        # DCBA|ABCD|DCBA
        map(idxs) do idx
            idx = mod1(idx, 2L)
            idx > L ? v[2L+1-idx] : v[idx]
        end
    else
        throw(ArgumentError("Unrecognized boundary condition $bc"))
    end

    # coeficients from http://paulbourke.net/miscellaneous/interpolation/
    # this is the Paul Breeuwsma / Catmull-Rom variant
    a0 = -0.5*ys[1] + 1.5*ys[2] - 1.5*ys[3] + 0.5*ys[4]
    a1 = ys[1] - 2.5*ys[2] + 2*ys[3] - 0.5*ys[4]
    a2 = -0.5*ys[1] + 0.5*ys[3]
    a3 = ys[2]

    i -= floor(i)
    a0*i^3 + a1*i^2 + a2*i + a3
end
