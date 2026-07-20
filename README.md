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

## Release v1.1
This release adds PELT optimisation to exact PCF (Pruned Exact Linear Time). This version has the same worst-case quadratic time complexity as Exact PCF, but approaches linear time in the best case. It does this by discovering positions that can never improve the overall segmentation score, and removing them from future iterations.
