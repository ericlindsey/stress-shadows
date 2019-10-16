function [Ks11,Kd11, ...
          Ks12,Kd12, ...
          Ks13,Kd13, ...
          Ks22,Kd22, ...
          Ks23,Kd23, ...
          Ks33,Kd33]=computeStressKernelsNikkhoo15(src,rcv,G,nu)
% computeStressKernelsNikkhoo15 computes the stress on a receiver faults due to
% the motion of triangular dislocations.
%
%   [Ks11,Kd11, ...
%    Ks12,Kd12, ...
%    Ks13,Kd13, ...
%    Ks22,Kd22, ...
%    Ks23,Kd23, ...
%    Ks33,Kd33]=COMPUTESTRESSKERNELSNIKKHOO15(src,rcv,G,nu)
%
% where nu is Poisson's ratio, G is the shear modulus. src is a structure
% with fault geometry and fault slip:
%
%   src.x         : position of vertices
%   src.vertices  : list of triangle vertices
%   src.sv        : patch strike unit vector (:,3)
%   src.dv        : patch dip unit vector (:,3)
%   src.nv        : patch normal unit vector (:,3)
%
% rcv is a structure with fault position and geometry:
%
%   rcv.xc : patch center position (:,3)
%
% OUTPUT:
%
%   Ksij is the stress components ij due to slip in the strike direction. 
%   Kdij is the stress components ij due to slip in the dip direction.
%
% Sylvain Barbot, 03/21/2015
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle

import unicycle.greens.*
import unicycle.utils.*

N=length(src.slip);
M=size(rcv.xc,1);

lambda=G*2*nu/(1-2*nu);

% initialize stress kernels
Ks11=zeros(M,N);
Ks12=zeros(M,N);
Ks13=zeros(M,N);
Ks22=zeros(M,N);
Ks23=zeros(M,N);
Ks33=zeros(M,N);

Kd11=zeros(M,N);
Kd12=zeros(M,N);
Kd13=zeros(M,N);
Kd22=zeros(M,N);
Kd23=zeros(M,N);
Kd33=zeros(M,N);

textprogressbar('# triangle  stress kernels: ');

for k=1:N
    
    if 0==mod(k-1,2)
        textprogressbar(k/N*100);
    end
    
    % strike component
    [S,~]=computeStressNikkhoo15(rcv.xc(:,1),rcv.xc(:,2),rcv.xc(:,3), ...
        src.x(src.vertices(k,1),:),src.x(src.vertices(k,2),:),src.x(src.vertices(k,3),:), ...
        1,0,0,G,lambda);
    
    % rotate stress to shear zone (primed) system of coordinates
    s11p= ( rcv.sv(:,1).*S(:,1)+rcv.sv(:,2).*S(:,4)+rcv.sv(:,3).*S(:,5)).*rcv.sv(:,1)+( rcv.sv(:,1).*S(:,4)+rcv.sv(:,2).*S(:,2)+rcv.sv(:,3).*S(:,6)).*rcv.sv(:,2)+( rcv.sv(:,1).*S(:,5)+rcv.sv(:,2).*S(:,6)+rcv.sv(:,3).*S(:,3)).*rcv.sv(:,3);
    s12p= ( rcv.sv(:,1).*S(:,1)+rcv.sv(:,2).*S(:,4)+rcv.sv(:,3).*S(:,5)).*rcv.nv(:,1)+( rcv.sv(:,1).*S(:,4)+rcv.sv(:,2).*S(:,2)+rcv.sv(:,3).*S(:,6)).*rcv.nv(:,2)+( rcv.sv(:,1).*S(:,5)+rcv.sv(:,2).*S(:,6)+rcv.sv(:,3).*S(:,3)).*rcv.nv(:,3);
    s13p=-( rcv.sv(:,1).*S(:,1)+rcv.sv(:,2).*S(:,4)+rcv.sv(:,3).*S(:,5)).*rcv.dv(:,1)-( rcv.sv(:,1).*S(:,4)+rcv.sv(:,2).*S(:,2)+rcv.sv(:,3).*S(:,6)).*rcv.dv(:,2)-( rcv.sv(:,1).*S(:,5)+rcv.sv(:,2).*S(:,6)+rcv.sv(:,3).*S(:,3)).*rcv.dv(:,3);
    s22p= ( rcv.nv(:,1).*S(:,1)+rcv.nv(:,2).*S(:,4)+rcv.nv(:,3).*S(:,5)).*rcv.nv(:,1)+( rcv.nv(:,1).*S(:,4)+rcv.nv(:,2).*S(:,2)+rcv.nv(:,3).*S(:,6)).*rcv.nv(:,2)+( rcv.nv(:,1).*S(:,5)+rcv.nv(:,2).*S(:,6)+rcv.nv(:,3).*S(:,3)).*rcv.nv(:,3);
    s23p=-( rcv.nv(:,1).*S(:,1)+rcv.nv(:,2).*S(:,4)+rcv.nv(:,3).*S(:,5)).*rcv.dv(:,1)-( rcv.nv(:,1).*S(:,4)+rcv.nv(:,2).*S(:,2)+rcv.nv(:,3).*S(:,6)).*rcv.dv(:,2)-( rcv.nv(:,1).*S(:,5)+rcv.nv(:,2).*S(:,6)+rcv.nv(:,3).*S(:,3)).*rcv.dv(:,3);
    s33p=-(-rcv.dv(:,1).*S(:,1)-rcv.dv(:,2).*S(:,4)-rcv.dv(:,3).*S(:,5)).*rcv.dv(:,1)-(-rcv.dv(:,1).*S(:,4)-rcv.dv(:,2).*S(:,2)-rcv.dv(:,3).*S(:,6)).*rcv.dv(:,2)-(-rcv.dv(:,1).*S(:,5)-rcv.dv(:,2).*S(:,6)-rcv.dv(:,3).*S(:,3)).*rcv.dv(:,3);
    
    % stress kernels in shear zone (primed) system of coordinates
    Ks11(:,k)=s11p;
    Ks12(:,k)=s12p;
    Ks13(:,k)=s13p;
    Ks22(:,k)=s22p;
    Ks23(:,k)=s23p;
    Ks33(:,k)=s33p;
    
    % dip component
    [S,~]=computeStressNikkhoo15(rcv.xc(:,1),rcv.xc(:,2),rcv.xc(:,3), ...
        src.x(src.vertices(k,1),:),src.x(src.vertices(k,2),:),src.x(src.vertices(k,3),:), ...
        0,-1,0,G,lambda);
    
    % rotate stress to shear zone (primed) system of coordinates
    s11p= ( rcv.sv(:,1).*S(:,1)+rcv.sv(:,2).*S(:,4)+rcv.sv(:,3).*S(:,5)).*rcv.sv(:,1)+( rcv.sv(:,1).*S(:,4)+rcv.sv(:,2).*S(:,2)+rcv.sv(:,3).*S(:,6)).*rcv.sv(:,2)+( rcv.sv(:,1).*S(:,5)+rcv.sv(:,2).*S(:,6)+rcv.sv(:,3).*S(:,3)).*rcv.sv(:,3);
    s12p= ( rcv.sv(:,1).*S(:,1)+rcv.sv(:,2).*S(:,4)+rcv.sv(:,3).*S(:,5)).*rcv.nv(:,1)+( rcv.sv(:,1).*S(:,4)+rcv.sv(:,2).*S(:,2)+rcv.sv(:,3).*S(:,6)).*rcv.nv(:,2)+( rcv.sv(:,1).*S(:,5)+rcv.sv(:,2).*S(:,6)+rcv.sv(:,3).*S(:,3)).*rcv.nv(:,3);
    s13p=-( rcv.sv(:,1).*S(:,1)+rcv.sv(:,2).*S(:,4)+rcv.sv(:,3).*S(:,5)).*rcv.dv(:,1)-( rcv.sv(:,1).*S(:,4)+rcv.sv(:,2).*S(:,2)+rcv.sv(:,3).*S(:,6)).*rcv.dv(:,2)-( rcv.sv(:,1).*S(:,5)+rcv.sv(:,2).*S(:,6)+rcv.sv(:,3).*S(:,3)).*rcv.dv(:,3);
    s22p= ( rcv.nv(:,1).*S(:,1)+rcv.nv(:,2).*S(:,4)+rcv.nv(:,3).*S(:,5)).*rcv.nv(:,1)+( rcv.nv(:,1).*S(:,4)+rcv.nv(:,2).*S(:,2)+rcv.nv(:,3).*S(:,6)).*rcv.nv(:,2)+( rcv.nv(:,1).*S(:,5)+rcv.nv(:,2).*S(:,6)+rcv.nv(:,3).*S(:,3)).*rcv.nv(:,3);
    s23p=-( rcv.nv(:,1).*S(:,1)+rcv.nv(:,2).*S(:,4)+rcv.nv(:,3).*S(:,5)).*rcv.dv(:,1)-( rcv.nv(:,1).*S(:,4)+rcv.nv(:,2).*S(:,2)+rcv.nv(:,3).*S(:,6)).*rcv.dv(:,2)-( rcv.nv(:,1).*S(:,5)+rcv.nv(:,2).*S(:,6)+rcv.nv(:,3).*S(:,3)).*rcv.dv(:,3);
    s33p=-(-rcv.dv(:,1).*S(:,1)-rcv.dv(:,2).*S(:,4)-rcv.dv(:,3).*S(:,5)).*rcv.dv(:,1)-(-rcv.dv(:,1).*S(:,4)-rcv.dv(:,2).*S(:,2)-rcv.dv(:,3).*S(:,6)).*rcv.dv(:,2)-(-rcv.dv(:,1).*S(:,5)-rcv.dv(:,2).*S(:,6)-rcv.dv(:,3).*S(:,3)).*rcv.dv(:,3);
    
    Kd11(:,k)=s11p;
    Kd12(:,k)=s12p;
    Kd13(:,k)=s13p;
    Kd22(:,k)=s22p;
    Kd23(:,k)=s23p;
    Kd33(:,k)=s33p;
end


textprogressbar(100);
textprogressbar('');

end
