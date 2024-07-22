segmentation
==

An R package implementing the piecewise constant function segmentation algorithms from
Nilsen et al. (https://doi.org/10.1186/1471-2164-13-591).

These algorithms are already available through the BioConductor `copynumber` package,
but have been reimplemented here because:
  - The intensive parts of the algorithm are now in C++ rather than R, which
  makes them faster and have lower memory requirements.
  - Parts of the algorithm were not implemented in the BioConductor package, mainly
  the `kmin` parameter (controlling the minimum size of a segment) was not available
  in "multipcf", but is provided here.
