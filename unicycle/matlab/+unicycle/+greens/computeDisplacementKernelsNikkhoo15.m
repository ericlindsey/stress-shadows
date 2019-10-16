function [Gs,Gd]=computeDisplacementKernelsNikkhoo15(triangle,nu,x,vecsize)
% COMPUTEDISPLACEMENTKERNELSNIKKHOO15 computes inversion kernels
% (Green's functions) to connect fault slip to surface deformation using
% the analytic solution of Nikkhoo et al. (2015) for triangular
% dislocations in a half space.
%
% INTERFACE:
%
%   G = unicycle.greens.computeDisplacementKernelsNikkhoo15(triangle,nu,x,vecsize)
%
% INPUT:
%
% triangle geometry.triangle object
% nu       Poisson's ratio
% x        coordinates of observation point (east, north, up)
%
% vecsize number of component for displacement vector.
%
%
% DATA LAYOUT:
%
%      --- strike slip ---; --- dip slip ---
%     /e1                                   \
%     |n1                                   |
%     |u1                                   |
%     |e2                                   |
%     |n2                                   |
%     |u2                                   |
% G = | .                                   |
%     | .                                   |
%     | .                                   |
%     | .                                   |
%     |en                                   |
%     |nn                                   |
%     \un                                   /
%
% SEE ALSO: unicycle, geometry.triangle
import unicycle.utils.*

assert(unicycle.greens.model.dgf==2,'greens.nikkhoo15.G assume 2=dgf');

% number of GPS stations
D=size(x,1);

% Green's function matrix


strikeSlip=[1,0,0];
dipSlip=[0,-1,0];
opening=[0,0,1];

KO=cell(2,1);

textprogressbar('# triangle displacement kernels: ');

% interlace east, north and up forward models
for j=1:2
    
    G=zeros(vecsize*D,triangle.N);
    
    for i=1:triangle.N
        
        if 0==mod(i-1,2)
            textprogressbar(i/triangle.N*50+50*(j-1));
        end
        
        % grouping coordinates of triangular dislocation vertices
        P1=triangle.x(triangle.vertices(i,1),:);
        P2=triangle.x(triangle.vertices(i,2),:);
        P3=triangle.x(triangle.vertices(i,3),:);
        
        [Ue,Un,Uu]=unicycle.greens.computeDisplacementNikkhoo15(x(:,1),x(:,2),x(:,3), ...
            P1,P2,P3,strikeSlip(j),dipSlip(j),opening(j),nu);
        
        switch(vecsize)
            case 2
                u=[Ue,Un]';
            case 3
                u=[Ue,Un,Uu]';
            otherwise
                error('greens:computeDisplacementKernelsNikkhoo15','invalid degrees of freedom (%d)\n',...
                    unicycle.greens.model.dgf);
        end
        
        G(:,i)=u(:);
        
    end
    
    KO{j}=G;
    
end

[Gs,Gd]=deal(KO{:});

textprogressbar(100);
textprogressbar('');

end % end G (Green's function)

