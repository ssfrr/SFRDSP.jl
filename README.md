# SFRDSP

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ssfrr.github.io/SFRDSP.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ssfrr.github.io/SFRDSP.jl/dev)
[![Build Status](https://travis-ci.org/ssfrr/SFRDSP.jl.svg?branch=master)](https://travis-ci.org/ssfrr/SFRDSP.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ssfrr/SFRDSP.jl?svg=true)](https://ci.appveyor.com/project/ssfrr/SFRDSP-jl)
[![Codecov](https://codecov.io/gh/ssfrr/SFRDSP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ssfrr/SFRDSP.jl)

SFRDSP.jl is a collection of DSP routines I've written, mostly for my personal use. They are generally not well-tested well-documented, or guaranteed to be stable, but you are welcome to use them if they are useful to you.

Some of these may make it into DSP.jl eventually, but the benefit for me of having my own routines is that I can write the minimal code I need for my use-cases and avoid getting too bogged down in generalizing.

Often when functionality from DSP.jl is duplicated here (usually because I have a different design or different priorities) I frequently affix a `2`, e.g. `stft2` for my short-time fourier transform implementation.

Feel free to submit PRs or issues if you find any bugs, but I will likely not spend very much time reviewing contributions that aren't useful to me personally.
