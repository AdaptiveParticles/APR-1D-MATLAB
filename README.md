# APR_1D_demo
This repository gives a didactic introduction to the Adaptive Particle Representation with simple 1D examples.

For an introduction or citation for the Adaptive Particle Representation please see: Cheeseman et al. 2018 (Forget Pixels: Adaptive Particle Representation of Fluorescence Microscopy Images)

## Usage:

Please include demo_source into your Matlab path (including subfolders).

Note: Some files require the use of simulink. 

## Demos Descriptions:

### demo_basic
*demo_apr_1D_naive_scheme_comparison: Computes the APR, using both a naive approach using Particle Cells, and the Pulling Scheme.

*demo_apr_1D_reconstruction: Computes the APR, comparing different reconstruction methods, and showing the satisfaction of the Reconstruction Bound.

*demo_equivalence_optimization: Compares forming the APR, with, and without, the Equivalence Optimization for the Pulling Scheme.

### demo_reconstruction_condition
*demo_apr_1D_change_E: Produces the APR for different relative error values E, and compares different reconstruction methods with the Reconstruction Bound.

### demo_symbolic_vs_numeric
*demo_apr_1D_numerical: Compares the solution of the APR using numerical estimates or the gradient, vs. analytic function

### demo_compare_continuous
*demo_compare_continuous: Compares the Implied Resolution Function R^*(y), Local Resolution Estimate L(y), with continuous optimal solutions, and Particle Cells. Reproduces plots from Figure2 of Cheeseman et al. 2018.

## Notes:

*The code assumes you are fimiliar with the naming concepts introduced in Cheeseman et al. 2018.

*Feel free to play with the functions, domain size, and bounds, and relative error E.

*Conceptually little changes in higher dimensions, with the exception of complexity of data-structures and access.






