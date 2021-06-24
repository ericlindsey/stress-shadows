Shadows: A stress-constrained inversion method for geodetic data 
------
Eric Lindsey, Earth Observatory of Singapore

Version 1.0: October 2019

The stress-shadows repository contains code and examples for running stress-constrained inversions for interseismic slip rate deficit (kinematic coupling) on faults. It is particularly applicable to offshore megathrusts, where near-trench resolution is low.

The method is based on a general-purpose inversion toolbox for solving problems of the type G\*m = d, subject to regularization L\*m = k, constraints A\*m >= 0, and bounds m1 <= m <= m2. The program is implemented in Object-oriented Matlab and currently contains modules for fault sources (based on the Unicycle toolkit) and geodetic data. Extending the toolbox to handle additional problem types can be accomplished by adding matlab functions that create the specified matrices.

The basic use of the toolbox involves creating one parameter file (matlab .m file) listing the names of source and dataset objects, and desired inversion parameters. An arbitrary number of sources and datasets may be specified together. 

See the examples for simplified instruction on the use of the toolbox.

**References:**

If you use this code, please cite our paper:

  Lindsey, E. O., R. Mallick, J. A. Hubbard, K. E. Bradley, R. Almeida, J. D. P. Moore, R. Burgmann, and E. M. Hill, Slip rate deficit and earthquake potential on shallow megathrusts, Nature Geoscience, [doi:10.1038/s41561-021-00736-x](https://doi.org/10.1038/s41561-021-00736-x), 2021.

Depending on which elastic deformation greens functions you use (rectangular or triangular), also cite:

 (rectangular) Okada, Y., 1992. Internal deformation due to shear and tensile faults in a half-space, Bull. Seism. Soc. Am., 82, 1018–1040.

 (triangular) Nikkhoo, M., and T. R. Walter (2015), Triangular dislocation: an analytical, artefact-free solution, Geophys. J. Int., 201(2), 1119–1141, https://doi.org/10.1093/gji/ggv035.
