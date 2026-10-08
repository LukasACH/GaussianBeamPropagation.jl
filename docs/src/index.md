# GaussianBeamPropagation.jl

!!! warning

    This package is very much a work in progress.
    The final API is not stable yet, and may change drastically over time until I am happy with at and stabilise it with a 1.0 release.

The goal of this package is to provide a framework to propagate Gaussian beams, e.g., Gaussian laser beams, through an optical system supporting general astigmatism.
It focusses on a fast propagation over perfect accuracy to enable numerical optimisation strategies to, e.g., minimise astigmatism in optical systems.

It implements extended ABCD matrices used with the [ray transfer matrix analysis](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis) that supports both [non-rotationally symmetric elements and third-order/oblique astigmatism](https://en.wikipedia.org/wiki/Astigmatism_optical_systems#Forms_of_astigmatism).
