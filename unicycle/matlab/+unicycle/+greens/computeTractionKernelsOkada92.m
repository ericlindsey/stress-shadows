function [Kss,Kds,Ksd,Kdd,Ksn,Kdn]=computeTractionKernelsOkada92(src,rcv,G,nu)
% COMPUTETRACTIONKERNELSOKADA92 computes the stress on a receiver faults due to
% motion of rectangular dislocations (Okada, 1992).
%
% SYNTAX:
%   [Kss,Ksd,Ksn,Kds,Kdd,Kdn]=COMPUTESTRESSKERNELSOKADA92(src,rcv,G,nu)
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
%   [Kss,Ksd,Ksn,Kds,Kdd,Kdn]=COMPUTESTRESSKERNELSOKADA92(src,rcv,G,nu)
%
% returns the shear stress kernels. Kss is the traction in the strike
% direction due to strike slip. Ksd is the traction in the dip direction
% due to strike slip. Ksn is the traction is the normal direction due to 
% strike slip. Kds is the traction in the strike direction due to
% dip slip and Kdd is the traction in the dip direction due to dip slip.
% Kdn is the traction in the normal direction due to dip slip.
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
Kss=zeros(M,N);
Ksd=zeros(M,N);
Ksn=zeros(M,N);
Kds=zeros(M,N);
Kdd=zeros(M,N);
Kdn=zeros(M,N);

textprogressbar('# rectangle traction kernels: ');

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
    
    % full traction vector on receiver faults
    t=[...
        S(:,1,1).*rcv.nv(:,1)+S(:,2,1).*rcv.nv(:,2)+S(:,3,1).*rcv.nv(:,3), ...
        S(:,1,2).*rcv.nv(:,1)+S(:,2,2).*rcv.nv(:,2)+S(:,3,2).*rcv.nv(:,3), ...
        S(:,1,3).*rcv.nv(:,1)+S(:,2,3).*rcv.nv(:,2)+S(:,3,3).*rcv.nv(:,3)];
    
    % normal traction vector on receiver faults
    tn=sum(t.*rcv.nv,2);
    
    % shear traction vector
    ts=t-repmat(tn,1,3).*rcv.nv;
    
    % shear stress in strike direction
    Kss(:,k)=sum(ts.*rcv.sv,2);
    
    % shear stress in dip direction
    Ksd(:,k)=sum(ts.*rcv.dv,2);
    
    % normal stress due to strike slip
    Ksn(:,k)=sum(ts.*rcv.nv,2);
    
    % dip component
    [~,~,~,S]=computeOkada92(1,xd(:),yd(:),rcv.xc(:,3),G,nu, ...
        -src.x(k,3),src.dip(k)/180*pi,src.L(k),src.W(k),'d',0,src.strike(k)/180*pi);
    S=reshape(S,M,3,3);
    
    % full traction vector on receiver faults
    t=[...
        S(:,1,1).*rcv.nv(:,1)+S(:,2,1).*rcv.nv(:,2)+S(:,3,1).*rcv.nv(:,3), ...
        S(:,1,2).*rcv.nv(:,1)+S(:,2,2).*rcv.nv(:,2)+S(:,3,2).*rcv.nv(:,3), ...
        S(:,1,3).*rcv.nv(:,1)+S(:,2,3).*rcv.nv(:,2)+S(:,3,3).*rcv.nv(:,3)];
    
    % normal traction vector on receiver faults
    tn=sum(t.*rcv.nv,2);
    
    % shear traction vector
    ts=t-repmat(tn,1,3).*rcv.nv;
    
    % shear stress in strike direction due to dip slip
    Kds(:,k)=sum(ts.*rcv.sv,2);
    
    % shear stress in dip direction due to dip slip
    Kdd(:,k)=sum(ts.*rcv.dv,2);
    
    % normal stress due to dip slip
    Kdn(:,k)=sum(ts.*rcv.nv,2);
    
end

textprogressbar(100);
textprogressbar('');

end
