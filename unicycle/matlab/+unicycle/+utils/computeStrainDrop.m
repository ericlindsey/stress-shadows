function depsilon=computeStrainDrop(sourceFilename)
% function depsilon = COMPUTESTRAINDROP(sourceFilename) computes the strain drop 
% depsilon from a coseismic slip distribution using the expression
%
%         /           
%         | n . E . s dA
%         /           
%   dE  = --------------
%             /
%             | s dA
%             /
%
% where n is the fault normal vector, E is the strain tensor, s is the
% coseismic slip on the fault.
%
% example:
%
%    unicycle.utils.computeStrainDrop('faults/hill+14_0089.flt')
%
% AUTHOR:
% Sylvain Barbot (April 10, 2014), Earth Observatory of Singapore
%
% SEE ALSO: unicycle

import unicycle.greens.computeStressKernelsOkada92
import unicycle.greens.okada92
import unicycle.geometry.source

% homogeneous elastic half space, Poisson's solid
earthModel=okada92(1/2,0);

% source and receiver faults
flt=source(sourceFilename,earthModel);
            
[Kss,Ksd,Kds,Kdd] = computeStressKernelsOkada92(flt,flt,1/2,0);

% stretch nhat . E . shat
nDotEStrike=Kss*(cosd(flt.rake).*flt.slip)+Kds*(sind(flt.rake).*flt.slip);
nDotEDip   =Ksd*(cosd(flt.rake).*flt.slip)+Kdd*(sind(flt.rake).*flt.slip);
stretch=nDotEStrike.*cosd(flt.rake)+nDotEDip.*sind(flt.rake);

% slip-averaged stress drop
depsilon=-sum(stretch.*flt.slip)/sum(flt.slip);

end
