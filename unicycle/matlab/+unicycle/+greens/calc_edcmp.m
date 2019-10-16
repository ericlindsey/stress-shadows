% Simulation of displacements due to fault slip in a layered half space
% using the semi-analytic solution of Wang et al. (2003).
%
% CALC_EDCMP('init',prefix)
%
% loads the point-source Green's function calculated with EDGRN. prefix is
% base name of the Green's functions, for example:
%
%     calc_edcmp('init','grnfcts/prem');
%
% loads the files prem.ss, prem.ds and prem.cl from the grnfcts directory.
% Files prem.ss, prem.ds and prem.cl cam be obtained using the edgrn 
% program.
%
% [ux,uy,uz]=CALC_EDCMP('displacement',slip,x1,x2,x3,...
%                  length,width,strike,dip,rake,x,y);
%
% computes the displacement components ux, uy, uz due to a fault defined
% using the Relax and edcmp conventions:
%
%                 N (x1,y)
%                /
%               /| Strike
%   x1,x2,x3 ->@------------------------      (x2,x)
%              |\        p .            \ W
%              :-\      i .              \ i
%              |  \    l .                \ d
%              :90 \  S .                  \ t
%              |-Dip\  .                    \ h
%              :     \. | Rake               \
%              |      -------------------------
%              :             L e n g t h
%              Z (x3,-z)
%
% displacements are calculated at positions (x, y) for east and north,
% respectively. x1 and y are oriented north. x2 and x are oriented east. x3
% and z are pointing in opposite directions: x3 is positive down and z is
% positive up.
%
% edcmp core calculations are implemented using C++, Fortran and OpenMP 
% (see edcmp.cpp and the edcmp library).
