function s=stress(src,rcv,G,nu)
% STRESS computes the stress on a receiver faults due to motion of
% dislocations:
%
%   s=STRESS(src,rcv,G,nu); 
%
% INPUT:
% 
% src      unicycle.geometry.source object
% rcv      unicycle.geometry.receiver object
% G        shear modulus
% nu       Poisson's ratio
%
% OUTPUT:
%
% s is a structure with components:
%
%   s.tau : the norm of the shear stress (column vector)
%   s.t   : the traction vector (:,3)
%   s.tn  : the normal stress (:,3)
%   s.ts  : the shear traction vector (:,3)
%   s.tss : the shear traction in the strike direction
%   s.tsd : the shear traction in the dip direction
%
%
%   src.slip   : slip distribution (column vector)
%   src.x      : patch upper left corner position
%   src.L      : patch length
%   src.W      : patch width
%   src.strike : patch strike
%   src.dip    : patch dip
%   src.rake   : patch rake
%   src.sv     : patch strike unit vector (:,3)
%   src.dv     : patch dip unit vector (:,3)
%   src.nv     : patch normal unit vector (:,3)
%
% rcv is a structure with fault geometry:
%
%   rcv.x      : patch upper left corner position (:,3)
%   rcv.xc     : patch center position (:,3)
%   rcv.L      : patch length (column vector)
%   rcv.W      : patch width
%   rcv.strike : patch strike
%   rcv.dip    : patch dip
%   rcv.sv     : patch strike unit vector (:,3)
%   rcv.dv     : patch dip unit vector (:,3)
%   rcv.nv     : patch normal unit vector (:,3)
%

import unicycle.greens.*
import unicycle.utils.*

% initialize stress
s=struct('t',0,'ts',0,'tn',0,'tss',0,'tsd',0,'tsr',0,'tau',0);

N=length(src.slip);
M=length(rcv.L);

textprogressbar('stress: ');

for k=1:N
    
    textprogressbar(k/N*100);
    
    % distance from the fault upper left corner
    xd=rcv.xc(:,1)-src.x(k,1);
    yd=rcv.xc(:,2)-src.x(k,2);
    
    % strike component
    [~,~,~,Ss]=computeOkada92(+src.slip(k),xd(:),yd(:),rcv.xc(:,3),G,nu, ...
        -src.x(k,3),src.dip(k)/180*pi,src.L(k),src.W(k),'s',0,src.strike(k)/180*pi);
    
    % dip component
    [~,~,~,Sd]=computeOkada92(+src.slip(k),xd(:),yd(:),rcv.xc(:,3),G,nu, ...
        -src.x(k,3),src.dip(k)/180*pi,src.L(k),src.W(k),'d',0,src.strike(k)/180*pi);
    
    % stress tensor
    S=reshape(cosd(src.rake(k))*Ss+sind(src.rake(k))*Sd,M,3,3);
    
    % full traction vector on receiver faults
    s.t=s.t+[...
        S(:,1,1).*rcv.nv(:,1)+S(:,2,1).*rcv.nv(:,2)+S(:,3,1).*rcv.nv(:,3), ...
        S(:,1,2).*rcv.nv(:,1)+S(:,2,2).*rcv.nv(:,2)+S(:,3,2).*rcv.nv(:,3), ...
        S(:,1,3).*rcv.nv(:,1)+S(:,2,3).*rcv.nv(:,2)+S(:,3,3).*rcv.nv(:,3)];
    s.t
end

% normal traction vector on receiver faults
s.tn=sum(s.t.*rcv.nv,2);

% shear traction vector
s.ts=s.t-repmat(s.tn,1,3).*rcv.nv;

% shear stress in strike direction
s.tss=sum(s.ts.*rcv.sv,2);

% shear stress in dip direction
s.tsd=sum(s.ts.*rcv.dv,2);

% shear stress in rake direction
s.tsr=cosd(rcv.rake).*s.tss+sind(rcv.rake).*s.tsd;

% norm of shear stress
s.tau=sqrt(sum(s.ts.*s.ts,2));

textprogressbar('');

end