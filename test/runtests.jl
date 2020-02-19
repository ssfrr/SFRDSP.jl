using SFRDSP
using Test
using DSP

@testset "SFRDSP.jl" begin
    @testset "stft" begin
        @testset "basic round-trip (L=$L)" for L in 10:17
            x = randn(L)
            X = stft2(x, 8, 8) # window centers at 1, 9, 17
            @test size(X) == (5,3)
            x2 = istft2(X, 8, 8)
            @test length(x2) == 17
            @test x2 ≈ [x; zeros(17-L)]
        end

        # TODO: overlap-add still doesn't work quite right with non-1/2 hop
        # sizes (need to handle the endpoints better)
        @testset "overlap-add (windows=($awin, $swin))" for (awin, swin, scale) in
            [(rect, rect, 2),
             (cosine, cosine, 1),
             (hanning, rect, 1),
             (rect, hanning, 1)]
            x = randn(32)
            X = stft2(x, 8, 4, window=awin)
            x2 = istft2(X, 8, 4, window=swin)
            @test x2 ≈ scale * [x; zeros(length(x2)-length(x))]
        end
    end

    @testset "xcorr2" begin
        # setting the length of `u` and `v` to 6 forces the desired nfft to be 11,
        # which gets bumped to 12 because of `nextfastfft`, so these test that all
        # the index twiddling still holds.
        @testset "basic cross-correlation" begin
            @test xcorr2([1; zeros(5)],
                         [1; zeros(5)];
                         center=false,
                         unbiased=false) ≈ [zeros(5); 1; zeros(5)]
            @test xcorr2([zeros(5); 1],
                         [1; zeros(5)];
                         center=false,
                         unbiased=false) ≈ [zeros(10); 1]
            @test xcorr2([1; zeros(5)],
                         [zeros(5); 1];
                         center=false,
                         unbiased=false) ≈ [1; zeros(10)]
        end
        @testset "xcorr with lagbounds" begin
            u = [1; zeros(5)]
            v = [1; zeros(5)]
            @test xcorr2(u,v; center=false, lagbounds=(-1,1)) ≈ [0, 1, 0]
        end
        @testset "xcorr with short u" begin
            u = [1, 0]
            v = [1, 0, 0]
            @test xcorr2(u,v; center=false, unbiased=false) ≈ [0, 0, 1, 0]
        end
        @testset "xcorr with short v" begin
            u = [1, 0, 0]
            v = [1, 0]
            @test xcorr2(u,v; center=false, unbiased=false) ≈ [0, 1, 0, 0]
        end
    end

end
