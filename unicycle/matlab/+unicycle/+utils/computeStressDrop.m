function dtau=computeStressDrop(sourceFilename,G,nu)
% function dtau = COMPUTESTRESSDROP(sourceFilename,G,nu) computes the stress drop
% dtau from a coseismic slip distribution using the stress interactions
% described by Okada (1992) following the formula
%
%         /          
%         | tau * s dA
%         /          
%   dtau = --------------
%             /
%             | s dA
%             /
%
% where tau is the amplitude of shear stress on the fault, s is the
% coseismic slip on the fault.
%
% example:
%
%    unicycle.utils.computeStressDrop('faults/hill+14_0089.flt',30e3,0.25)
%
% AUTHOR:
% Sylvain Barbot (April 10, 2014), Earth Observatory of Singapore
%
% SEE ALSO: unicycle, hmmvp.computeStressDrop
 
import unicycle.greens.computeStressKernelsOkada92
import unicycle.greens.okada92
import unicycle.geometry.source
 
% homogeneous elastic half space, Poisson's solid
earthModel=okada92(G,nu);
 
% source and receiver faults
flt = source(sourceFilename,earthModel);
           
[Kss,Ksd,Kds,Kdd] = computeStressKernelsOkada92(flt,flt,G,nu);
 
% shear stress change
ts=Kss*(cosd(flt.rake).*flt.slip)+Kds*(sind(flt.rake).*flt.slip);
td=Ksd*(cosd(flt.rake).*flt.slip)+Kdd*(sind(flt.rake).*flt.slip);
tau=(ts.*cosd(flt.rake)+td.*sind(flt.rake));

% slip-averaged stress drop
dtau=-sum(flt.slip.*tau)/sum(flt.slip);

end
