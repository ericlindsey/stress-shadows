% Unicycle (Unified Cycle of Earthquakes) is a framework for the simulation
% of fault slip using rate- and state-dependent friction with the radiation
% damping approximation.
%
% Geometry
%   unicycle/geometry/coseismicPatch         - coseismic slip distribution
%   unicycle/geometry/patch                  - fault patch
%   unicycle/geometry/passiveReceiver        - passive stressed faults (no motion)
%   unicycle/geometry/source                 - moving fault loading receiver faults
%   unicycle/geometry/receiver               - receiver fault with slip evolution
%   unicycle/geometry/segment                - fault segment grouping patches
%   unicycle/geometry/triangle               - triangular fault patch
%   unicycle/geometry/triangleReceiver       - triangular receiver w/ slip evolution
%
% Green's function
%   unicycle/greens/stress                   - stress on receiver faults
%   unicycle/greens/stressKernels            - Okada (1992) stress kernels
%   unicycle/greens/triangleStressKernels    - Meade (2007) stress kernels
%   unicycle/greens/okada85                  - displacement at surface
%
%   hmmvp/unihmmvp                           - Green's functions with Hierarchical Matrices
%   
% Manifold
%   unicycle/manifold/gps                    - create a gps object
%
% Ordinary Differential Equation
%   unicycle/ode/evolution                   - models of fault slip evolution
%   unicycle/ode/rateandstate                - rate-and-state friction
%   unicycle/ode/rateandstatedamping         - rate-and-state friction w/ radiation damping
%   unicycle/ode/ratestrengthening           - rate-strengthening approximation
%   unicycle/ode/ratestrengthening_prestress - rate-strengthening friction w/ pre-stress
%
% Input/Output and formats
%   unicycle/export/exportflt_rfaults        - compatible w/ Relax, EDCMP, Gamra
%   unicycle/export/exportvtk_rfaults        - 3D visualization w/ Paraview
%   unicycle/export/exportxyz_rfaults        - GMT ASCII format
%   unicycle/export/grdread                  - read GMT GRD format
%
% Optimization
%   unicycle/optim/sim_anl                   - simulated annealing
%   unicycle/optim/mh                        - Metropolis Hastings
%   unicycle/optim/laplacian                 - smoothing matrix operator
%
% Author:
%   Sylvain Barbot, Earth Observatory of Singapore, May 24, 2013
