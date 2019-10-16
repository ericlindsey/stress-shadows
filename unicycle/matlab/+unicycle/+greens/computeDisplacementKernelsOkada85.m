function [Ks,Kd]=computeDisplacementKernelsOkada85(patch,nu,x,vecsize)
% COMPUTEDISPLACEMENTKERNELSOKADA85 computes inversion kernels
% (Green's functions) to connect fault slip to surface deformation using 
% the analytic solution of Okada (1985) for a half space.
%
% INTERFACE:
%
%   [Ks,Kd] = unicycle.greens.computeDisplacementKernelsOkada85(src,nu,x,vecsize)
%
% INPUT:
%
% src     geometry.patch object
% nu      Poisson's ratio
% x       coordinates of observation point (east, north, up)
%         *** but point should be at the surface ***
% vecsize number of component for displacement vector
%
%
% DATA LAYOUT:
%
%          strike slip                dip slip 
%      /e1               \      /e1               \
%      |n1               |      |n1               |
%      |u1               |      |u1               |
%      |e2               |      |e2               |
%      |n2               |      |n2               |
%      |u2               |      |u2               |
% Ks = | .               | Kd = | .               |
%      | .               |      | .               |
%      | .               |      | .               |
%      | .               |      | .               |
%      |en               |      |en               |
%      |nn               |      |nn               |
%      \un               /      \un               /
%
% SEE ALSO: unicycle, greens.okada92

import unicycle.greens.*
import unicycle.utils.*

assert(unicycle.greens.model.dgf==2,'greens.okada92.G assume 2=dgf');

% number of GPS stations
D=size(x,1);

% Green's function matrix
K=cell(2,1);

textprogressbar('# rectangle displacement kernels: ');

for j=1:2
    
    G=zeros(vecsize*D,patch.N);
    
    for i=1:patch.N
        % east coordinate of fault upper corner
        xg=patch.x(i,1);
        % north coordinate of fault upper corner
        yg=patch.x(i,2);
        % depth coordinate of fault upper corner
        zg=-patch.x(i,3);
        
        % length of fault path (in strike direction)
        L=patch.L(i);
        % width of fault path (in dip direction)
        W=patch.W(i);
        
        % orientation
        strike=patch.strike(i)/180*pi;
        dip=patch.dip(i)/180*pi;
        
        % relative distance between source and receiver
        xd=x(:,1)-xg;
        yd=x(:,2)-yg;
        
        [ux,uy,uz]=computeDisplacementOkada85(1,xd,yd,nu,dip,zg,L,W,j,strike);
        
        switch(vecsize)
            case 2
                u=[ux,uy]';
            case 3
                u=[ux,uy,uz]';
            otherwise
                error('greens:computeDisplacementKernelsOkada85','invalid degrees of freedom (%d)\n',...
                    unicycle.greens.model.dgf);
        end
       
        G(:,i)=u(:);
        
    end
    
    K{j}=G;
    
end

% distribute kernels into two separate variables.
[Ks,Kd]=deal(K{:});

textprogressbar(100);
textprogressbar('');

end



