function [L11,L12,L13,L22,L23,L33]=computeDisplacementKernelsVerticalShearZone(shz,nu,x,vecsize)
% COMPUTEDISPLACEMENTKERNELSVERTICALSHEARZONE computes inversion kernels
% (Green's functions) to connect fault slip to surface deformation using 
% the analytic solution of Barbot, Moore and Lambert (2016).
%
% INTERFACE:
%
%   [L11,L12,L13,L22,L23,L33] = ...
%            unicycle.greens.computeDisplacementKernelsVerticalShearZone(src,nu,x,vecsize)
%
% INPUT:
%
% shz     geometry.shearZone object
% nu      Poisson's ratio
% x       coordinates of observation point (east, north, up)
% vecsize number of component for displacement vector
%
%
% DATA LAYOUT:
%
%            strain eij
%       /e1               \
%       |n1               |
%       |u1               |
%       |e2               |
%       |n2               |
%       |u2               |
% Lij = | .               |
%       | .               |
%       | .               |
%       | .               |
%       |en               |
%       |nn               |
%       \un               /
%
% Sylvain Barbot and James D. P. Moore, 07/16/2016
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle, greens.shearZone16
import unicycle.utils.*

% number of GPS stations
D=size(x,1);

% Green's function matrix
L=cell(6,1);

% individual strain components
e=eye(6);

textprogressbar('# shear zone displacement kernels: ');

for j=1:6

    G=zeros(vecsize*D,shz.N);
    
    for i=1:shz.N
        
        if 0==mod(i-1,2)
            textprogressbar((i/shz.N+(j-1))*100/6);
        end
        [u1,u2,u3]=unicycle.greens.computeDisplacementVerticalShearZone( ...
            x(:,2),x(:,1),-x(:,3), ...
            shz.x(i,2),shz.x(i,1),-shz.x(i,3), ...
            shz.L(i),shz.T(i),shz.W(i),shz.strike(i), ...
            e(j,1),e(j,2),e(j,3),e(j,4),e(j,5),e(j,6),1,nu);
            
        switch(vecsize)
            case 2
                u=[u2,u1]';
            case 3
                u=[u2,u1,-u3]';
            otherwise
                error('greens:copmuteDisplacementKernelsOkada85','invalid degrees of freedom (%d)\n',...
                    unicycle.greens.model.dgf);
        end
       
        G(:,i)=u(:);
        
    end
    
    L{j,1}=G;
    
    
end

% distribute kernels into two separate variables.
[L11,L12,L13,L22,L23,L33]=deal(L{:});

textprogressbar(100);
textprogressbar('');

end



