function [Ks11,Kd11, ...
          Ks12,Kd12, ...
          Ks13,Kd13, ...
          Ks22,Kd22, ...
          Ks23,Kd23, ...
          Ks33,Kd33]=computeStressKernelsOkada92(src,rcv,G,nu)
% COMPUTESTRESSKERNELSOKADA92 computes the stress on a receiver faults due to
% motion of rectangular dislocations (Okada, 1992).
%
% SYNTAX:
%   [Ks11,Ks12,Ks13,Ks22,Ks23,Ks33, ...
%    Kd11,Kd12,Kd13,Kd22,Kd23,Kd33]=COMPUTESTRESSKERNELSOKADA92(src,rcv,G,nu)
%
% where nu is Poisson's ratio, G is the shear modulus. src is a structure
% with fault geometry and fault slip:
%
%   src.slip   : slip distribution (column vector)
%   src.x      : patch upper left corner position
%   src.L      : patch length
%   src.W      : patch width
%   src.strike : patch strike
%   src.dip    : patch dip
%   src.rake   : patch rake
%
% rcv is a structure with at least:
%
%   rcv.xc     : patch center position (:,3)
%   rcv.sv     : patch strike unit vector (:,3)
%   rcv.dv     : patch dip unit vector (:,3)
%   rcv.nv     : patch normal unit vector (:,3)
%
% OUTPUT:
%
%   [Ks11,Ks12,Ks13,Ks22,Ks23,Ks33, ...
%    Kd11,Kd12,Kd13,Kd22,Kd23,Kd33]=COMPUTESTRESSKERNELSOKADA92(src,rcv,G,nu)
%
% returns the shear stress kernels. Ksij is the stress component Sij due to
% strike slip. Kdij is the stress component Sij due to strike slip. Sij is
% in the shear-zone centric system of coordinates.
%
% Sylvain Barbot, 07/18/2016
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle

import unicycle.greens.*
import unicycle.utils.*

N=length(src.slip);
M=size(rcv.xc,1);

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

textprogressbar('# rectangle stress kernels: ');

for k=1:N
    
    if 0==mod(k-1,2)
        textprogressbar(k/N*100);
    end
    
    % distance from the fault upper left corner
    xd=rcv.xc(:,1)-src.x(k,1);
    yd=rcv.xc(:,2)-src.x(k,2);
    
    % strike component
    [~,~,~,S]=computeOkada92(1,xd(:),yd(:),rcv.xc(:,3),G,nu, ...
        -src.x(k,3),src.dip(k)/180*pi,src.L(k),src.W(k),'s',0,src.strike(k)/180*pi);
    S=reshape(S,M,3,3);
    
    % rotate tensor to shear-zone centric (primed) system of coordinates
    s11p= ( rcv.sv(:,1).*S(:,1,1)+rcv.sv(:,2).*S(:,1,2)+rcv.sv(:,3).*S(:,1,3)).*rcv.sv(:,1)+( rcv.sv(:,1).*S(:,1,2)+rcv.sv(:,2).*S(:,2,2)+rcv.sv(:,3).*S(:,2,3)).*rcv.sv(:,2)+( rcv.sv(:,1).*S(:,1,3)+rcv.sv(:,2).*S(:,2,3)+rcv.sv(:,3).*S(:,3,3)).*rcv.sv(:,3);
    s12p= ( rcv.sv(:,1).*S(:,1,1)+rcv.sv(:,2).*S(:,1,2)+rcv.sv(:,3).*S(:,1,3)).*rcv.nv(:,1)+( rcv.sv(:,1).*S(:,1,2)+rcv.sv(:,2).*S(:,2,2)+rcv.sv(:,3).*S(:,2,3)).*rcv.nv(:,2)+( rcv.sv(:,1).*S(:,1,3)+rcv.sv(:,2).*S(:,2,3)+rcv.sv(:,3).*S(:,3,3)).*rcv.nv(:,3);
    s13p=-( rcv.sv(:,1).*S(:,1,1)+rcv.sv(:,2).*S(:,1,2)+rcv.sv(:,3).*S(:,1,3)).*rcv.dv(:,1)-( rcv.sv(:,1).*S(:,1,2)+rcv.sv(:,2).*S(:,2,2)+rcv.sv(:,3).*S(:,2,3)).*rcv.dv(:,2)-( rcv.sv(:,1).*S(:,1,3)+rcv.sv(:,2).*S(:,2,3)+rcv.sv(:,3).*S(:,3,3)).*rcv.dv(:,3);
    s22p= ( rcv.nv(:,1).*S(:,1,1)+rcv.nv(:,2).*S(:,1,2)+rcv.nv(:,3).*S(:,1,3)).*rcv.nv(:,1)+( rcv.nv(:,1).*S(:,1,2)+rcv.nv(:,2).*S(:,2,2)+rcv.nv(:,3).*S(:,2,3)).*rcv.nv(:,2)+( rcv.nv(:,1).*S(:,1,3)+rcv.nv(:,2).*S(:,2,3)+rcv.nv(:,3).*S(:,3,3)).*rcv.nv(:,3);
    s23p=-( rcv.nv(:,1).*S(:,1,1)+rcv.nv(:,2).*S(:,1,2)+rcv.nv(:,3).*S(:,1,3)).*rcv.dv(:,1)-( rcv.nv(:,1).*S(:,1,2)+rcv.nv(:,2).*S(:,2,2)+rcv.nv(:,3).*S(:,2,3)).*rcv.dv(:,2)-( rcv.nv(:,1).*S(:,1,3)+rcv.nv(:,2).*S(:,2,3)+rcv.nv(:,3).*S(:,3,3)).*rcv.dv(:,3);
    s33p=-(-rcv.dv(:,1).*S(:,1,1)-rcv.dv(:,2).*S(:,1,2)-rcv.dv(:,3).*S(:,1,3)).*rcv.dv(:,1)-(-rcv.dv(:,1).*S(:,1,2)-rcv.dv(:,2).*S(:,2,2)-rcv.dv(:,3).*S(:,2,3)).*rcv.dv(:,2)-(-rcv.dv(:,1).*S(:,1,3)-rcv.dv(:,2).*S(:,2,3)-rcv.dv(:,3).*S(:,3,3)).*rcv.dv(:,3);
    
    % stress kernels due to strike slip
    Ks11(:,k)=s11p;
    Ks12(:,k)=s12p;
    Ks13(:,k)=s13p;
    Ks22(:,k)=s22p;
    Ks23(:,k)=s23p;
    Ks33(:,k)=s33p;
    
    % dip component
    [~,~,~,S]=computeOkada92(1,xd(:),yd(:),rcv.xc(:,3),G,nu, ...
        -src.x(k,3),src.dip(k)/180*pi,src.L(k),src.W(k),'d',0,src.strike(k)/180*pi);
    S=reshape(S,M,3,3);
    
    % rotate tensor to shear-zone centric (primed) system of coordinates
    s11p= ( rcv.sv(:,1).*S(:,1,1)+rcv.sv(:,2).*S(:,1,2)+rcv.sv(:,3).*S(:,1,3)).*rcv.sv(:,1)+( rcv.sv(:,1).*S(:,1,2)+rcv.sv(:,2).*S(:,2,2)+rcv.sv(:,3).*S(:,2,3)).*rcv.sv(:,2)+( rcv.sv(:,1).*S(:,1,3)+rcv.sv(:,2).*S(:,2,3)+rcv.sv(:,3).*S(:,3,3)).*rcv.sv(:,3);
    s12p= ( rcv.sv(:,1).*S(:,1,1)+rcv.sv(:,2).*S(:,1,2)+rcv.sv(:,3).*S(:,1,3)).*rcv.nv(:,1)+( rcv.sv(:,1).*S(:,1,2)+rcv.sv(:,2).*S(:,2,2)+rcv.sv(:,3).*S(:,2,3)).*rcv.nv(:,2)+( rcv.sv(:,1).*S(:,1,3)+rcv.sv(:,2).*S(:,2,3)+rcv.sv(:,3).*S(:,3,3)).*rcv.nv(:,3);
    s13p=-( rcv.sv(:,1).*S(:,1,1)+rcv.sv(:,2).*S(:,1,2)+rcv.sv(:,3).*S(:,1,3)).*rcv.dv(:,1)-( rcv.sv(:,1).*S(:,1,2)+rcv.sv(:,2).*S(:,2,2)+rcv.sv(:,3).*S(:,2,3)).*rcv.dv(:,2)-( rcv.sv(:,1).*S(:,1,3)+rcv.sv(:,2).*S(:,2,3)+rcv.sv(:,3).*S(:,3,3)).*rcv.dv(:,3);
    s22p= ( rcv.nv(:,1).*S(:,1,1)+rcv.nv(:,2).*S(:,1,2)+rcv.nv(:,3).*S(:,1,3)).*rcv.nv(:,1)+( rcv.nv(:,1).*S(:,1,2)+rcv.nv(:,2).*S(:,2,2)+rcv.nv(:,3).*S(:,2,3)).*rcv.nv(:,2)+( rcv.nv(:,1).*S(:,1,3)+rcv.nv(:,2).*S(:,2,3)+rcv.nv(:,3).*S(:,3,3)).*rcv.nv(:,3);
    s23p=-( rcv.nv(:,1).*S(:,1,1)+rcv.nv(:,2).*S(:,1,2)+rcv.nv(:,3).*S(:,1,3)).*rcv.dv(:,1)-( rcv.nv(:,1).*S(:,1,2)+rcv.nv(:,2).*S(:,2,2)+rcv.nv(:,3).*S(:,2,3)).*rcv.dv(:,2)-( rcv.nv(:,1).*S(:,1,3)+rcv.nv(:,2).*S(:,2,3)+rcv.nv(:,3).*S(:,3,3)).*rcv.dv(:,3);
    s33p=-(-rcv.dv(:,1).*S(:,1,1)-rcv.dv(:,2).*S(:,1,2)-rcv.dv(:,3).*S(:,1,3)).*rcv.dv(:,1)-(-rcv.dv(:,1).*S(:,1,2)-rcv.dv(:,2).*S(:,2,2)-rcv.dv(:,3).*S(:,2,3)).*rcv.dv(:,2)-(-rcv.dv(:,1).*S(:,1,3)-rcv.dv(:,2).*S(:,2,3)-rcv.dv(:,3).*S(:,3,3)).*rcv.dv(:,3);
    
    % stress kernels due to dip slip
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
