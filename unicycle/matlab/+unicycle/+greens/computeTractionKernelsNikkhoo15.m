function [Kss,Kds,Ksd,Kdd,Ksn,Kdn]=computeTractionKernelsNikkhoo15(src,rcv,G,nu)
% computeTractionKernelsNikkhoo15 computes the stress on a receiver faults due to
% the motion of triangular dislocations.
%
%   [Kss,Ksd,Ksn,Kds,Kdd,Kdn]=COMPUTETRACTIONKERNELSNIKKHOO15(src,rcv,G,nu)
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
%   rcv.sv : patch strike unit vector (:,3)
%   rcv.dv : patch dip unit vector (:,3)
%   rcv.nv : patch normal unit vector (:,3)
%
% OUTPUT:
%
%   [Kss,Ksd,Ksn,Kds,Kdd,Kdn]=TRIANGLESTRESSKERNELS(src,rcv,G,nu)
%
% returns the shear stress kernels. Kss is the traction in the strike
% direction due to strike slip. Ksd is the traction in the dip direction
% due to strike slip. Kds is the tarction in the strike direction due to
% dip slip and Kdd is the traction in the dip direction due to dip slip.
%
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
Kss=zeros(M,N);
Ksd=zeros(M,N);
Ksn=zeros(M,N);
Kds=zeros(M,N);
Kdd=zeros(M,N);
Kdn=zeros(M,N);

textprogressbar('# triangle traction kernels: ');

for k=1:N
    
    if 0==mod(k-1,2)
        textprogressbar(k/N*100);
    end
    
    % strike component
    [S,~]=computeStressNikkhoo15(rcv.xc(:,1),rcv.xc(:,2),rcv.xc(:,3), ...
        src.x(src.vertices(k,1),:),src.x(src.vertices(k,2),:),src.x(src.vertices(k,3),:), ...
        1,0,0,G,lambda);
    
    % full traction vector on receiver faults
    t=[...
        S(:,1).*rcv.nv(:,1)+S(:,4).*rcv.nv(:,2)+S(:,5).*rcv.nv(:,3), ...
        S(:,4).*rcv.nv(:,1)+S(:,2).*rcv.nv(:,2)+S(:,6).*rcv.nv(:,3), ...
        S(:,5).*rcv.nv(:,1)+S(:,6).*rcv.nv(:,2)+S(:,3).*rcv.nv(:,3)];
    
    % shear stress in strike direction
    Kss(:,k)=sum(t.*rcv.sv,2);
    
    % shear stress in dip direction
    Ksd(:,k)=sum(t.*rcv.dv,2);
    
    % stress in normal direction
    Ksn(:,k)=sum(t.*rcv.nv,2);
    
    % dip component
    [S,~]=computeStressNikkhoo15(rcv.xc(:,1),rcv.xc(:,2),rcv.xc(:,3), ...
        src.x(src.vertices(k,1),:),src.x(src.vertices(k,2),:),src.x(src.vertices(k,3),:), ...
        0,-1,0,G,lambda);
    
    % full traction vector on receiver faults
    t=[...
        S(:,1).*rcv.nv(:,1)+S(:,4).*rcv.nv(:,2)+S(:,5).*rcv.nv(:,3), ...
        S(:,4).*rcv.nv(:,1)+S(:,2).*rcv.nv(:,2)+S(:,6).*rcv.nv(:,3), ...
        S(:,5).*rcv.nv(:,1)+S(:,6).*rcv.nv(:,2)+S(:,3).*rcv.nv(:,3)];
    
    % shear stress in strike direction due to dip slip
    Kds(:,k)=sum(t.*rcv.sv,2);
    
    % shear stress in dip direction due to dip slip
    Kdd(:,k)=sum(t.*rcv.dv,2);
    
    % stress in normal direction due to dip slip
    Kdn(:,k)=sum(t.*rcv.nv,2);
    
end

textprogressbar(100);
textprogressbar('');

end
