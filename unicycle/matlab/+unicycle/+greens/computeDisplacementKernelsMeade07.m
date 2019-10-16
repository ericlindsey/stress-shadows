function G=computeDisplacementKernelsMeade07(triangle,nu,x,vecsize)
% COMPUTEDISPLACEMENTKERNELSMEADE07 computes inversion kernels
% (Green's functions) to connect fault slip to surface
% deformation using the analytic solution of Meade (2007) for
% triangular dislocations in a half space.
%
% INTERFACE:
%
%   G = unicycle.greens.computeDisplacementKernelsMeade07(triangle,nu,x,vecsize)
%
% INPUT:
%
% triangle geometry.triangle object
% nu       Poisson's ratio
% x        coordinates of observation point (east, north, up)
%          *** but point should be at the surface ***
% vecsize number of component for displacement vector
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

assert(unicycle.greens.model.dgf==2,'greens.meade07.G assume 2=dgf');

% number of GPS stations
D=size(x,1);

% Green's function matrix
G=zeros(vecsize*D,unicycle.greens.model.dgf*triangle.N);

strikeSlip=[-1,0,0];
dipSlip=[0,-1,0];
opening=[0,0,1];

textprogressbar('# triangle displacement kernels: ');

for i=1:triangle.N
    
    % interlace east, north and up forward models
    for j=1:unicycle.greens.model.dgf
        
        U=unicycle.greens.computeDisplacementMeade07(x(:,1),x(:,2),-x(:,3), ...
            triangle.x(triangle.vertices(i,:),1),triangle.x(triangle.vertices(i,:),2),-triangle.x(triangle.vertices(i,:),3), ...
            nu,strikeSlip(j),opening(j),dipSlip(j));
        
        switch(vecsize)
            case 2
                u=[U.x,U.y];
            case 3
                u=[U.x,U.y,-U.z];
            otherwise
                error('greens:computeDisplacementKernelsMeade07','invalid degrees of freedom (%d)\n',...
                    unicycle.greens.model.dgf);
        end
        for l=1:vecsize
            % strike slip and dip slip are not interleaved
            G(l:vecsize:end,i+(j-1)*triangle.N)=u(:,l);
        end
    end
end

textprogressbar(100);
textprogressbar('');

end % end G (Green's function)

