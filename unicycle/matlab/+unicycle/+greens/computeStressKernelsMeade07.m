function [varargout]=computeStressKernelsMeade07(src,rcv,G,nu)
% computeStressKernelsMeade07 computes the stress on a receiver faults due to
% the motion of triangular dislocations.
%
%   [varargout]=COMPUTESTRESSKERNELSMEADE07(src,rcv,G,nu)
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
%   [Kss,Ksd,Kds,Kdd]=TRIANGLESTRESSKERNELS(src,rcv,G,nu)
%
% returns the shear stress kernels. Kss is the traction in the strike
% direction due to strike slip. Ksd is the traction in the dip direction
% due to strike slip. Kds is the tarction in the strike direction due to
% dip slip and Kdd is the traction in the dip direction due to dip slip.
%
%   [Kss,Ksd,Kds,Kdd,SS]=TRIANGLESTRESSKERNELS(src,rcv,G,nu)
%
%   [Kss,Ksd,Kds,Kdd,SS,DS]=TRIANGLESTRESSKERNELS(src,rcv,G,nu)
%
% provides all the stress tensor components Sxx, Sxy, Sxz, Syy, Syz, Szz
% and the normal tractions tn. SS contains the stress changes due to strike
% slip and DS contains the stress changes due to dip slip. The two
% variables are structures with SS.xx, SS.xy, SS.xz, SS.yy, SS.yz, SS.zz
% and SS.tn.
%
% Sylvain Barbot, 03/21/2015
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle

import unicycle.greens.*
import unicycle.utils.*

N=length(src.slip);
M=size(rcv.xc,1);

% initialize stress kernels
if (0<nargout), Kss=zeros(M,N); end;
if (1<nargout), Ksd=zeros(M,N); end;
if (2<nargout), Ksn=zeros(M,N); end;
if (3<nargout), Kds=zeros(M,N); end;
if (4<nargout), Kdd=zeros(M,N); end;
if (5<nargout), Kdn=zeros(M,N); end;

% stress components due to strike slip
if (6<nargout)
    SS.Kxx=zeros(M,N);
    SS.Kxy=zeros(M,N);
    SS.Kxz=zeros(M,N);
    SS.Kyy=zeros(M,N);
    SS.Kyz=zeros(M,N);
    SS.Kzz=zeros(M,N);
    SS.Ktn=zeros(M,N);
end

% stress components due to dip slip
if (7<nargout)
    DS.Kxx=zeros(M,N);
    DS.Kxy=zeros(M,N);
    DS.Kxz=zeros(M,N);
    DS.Kyy=zeros(M,N);
    DS.Kyz=zeros(M,N);
    DS.Kzz=zeros(M,N);
    DS.Ktn=zeros(M,N);
end

textprogressbar('# triangle  stress kernels: ');

for k=1:N
    
    if 0==mod(k-1,2)
        textprogressbar(k/N*100);
    end
    
    % strike component
    E=computeStrainMeade07(rcv.xc(:,1),rcv.xc(:,2),-rcv.xc(:,3), ...
        src.x(src.vertices(k,:),1),src.x(src.vertices(k,:),2),-src.x(src.vertices(k,:),3), ...
        nu,-1,0,0);
    S=hookesLaw(E,G,nu);
    
    % full traction vector on receiver faults
    t=[...
        S.xx.*rcv.nv(:,1)+S.xy.*rcv.nv(:,2)+S.xz.*rcv.nv(:,3), ...
        S.xy.*rcv.nv(:,1)+S.yy.*rcv.nv(:,2)+S.yz.*rcv.nv(:,3), ...
        S.xz.*rcv.nv(:,1)+S.yz.*rcv.nv(:,2)+S.zz.*rcv.nv(:,3)];
    
    % normal traction vector on receiver faults
    tn=sum(t.*rcv.nv,2);
    
    % shear stress in strike direction
    if (0<nargout), Kss(:,k)=sum(t.*rcv.sv,2); end
    
    % shear stress in dip direction
    if (1<nargout), Ksd(:,k)=sum(t.*rcv.dv,2); end
    
    % normal stress due to strike slip
    if (2<nargout), Ksn(:,k)=sum(t.*rcv.nv,2); end
    
    if (6<nargout)
        SS.Kxx(:,k)=S.xx;
        SS.Kxy(:,k)=S.xy;
        SS.Kxz(:,k)=S.xz;
        SS.Kyy(:,k)=S.yy;
        SS.Kyz(:,k)=S.yz;
        SS.Kzz(:,k)=S.zz;
        SS.Ktn(:,k)=tn;
    end
    
    % dip component
    E=computeStrainMeade07(rcv.xc(:,1),rcv.xc(:,2),-rcv.xc(:,3), ...
        src.x(src.vertices(k,:),1),src.x(src.vertices(k,:),2),-src.x(src.vertices(k,:),3), ...
        nu,0,0,1);
    S=hookesLaw(E,G,nu);
    
    % full traction vector on receiver faults
    t=[...
        S.xx.*rcv.nv(:,1)+S.xy.*rcv.nv(:,2)+S.xz.*rcv.nv(:,3), ...
        S.xy.*rcv.nv(:,1)+S.yy.*rcv.nv(:,2)+S.yz.*rcv.nv(:,3), ...
        S.xz.*rcv.nv(:,1)+S.yz.*rcv.nv(:,2)+S.zz.*rcv.nv(:,3)];
    
    % normal traction vector on receiver faults
    tn=sum(t.*rcv.nv,2);
    
    % shear traction vector
    ts=t-repmat(tn,1,3).*rcv.nv;
    
    % shear stress in strike direction due to dip slip
    if (3<nargout), Kds(:,k)=sum(t.*rcv.sv,2); end
    
    % shear stress in dip direction due to dip slip
    if (4<nargout), Kdd(:,k)=sum(t.*rcv.dv,2); end
    
    % normal stress due to dip slip
    if (5<nargout), Kdn(:,k)=sum(t.*rcv.nv,2); end
    
    if (7<nargout)
        DS.Kxx(:,k)=S.xx;
        DS.Kxy(:,k)=S.xy;
        DS.Kxz(:,k)=S.xz;
        DS.Kyy(:,k)=S.yy;
        DS.Kyz(:,k)=S.yz;
        DS.Kzz(:,k)=S.zz;
        DS.Ktn(:,k)=tn;
    end
end

if (0<nargout), varargout{1}=Kss; end
if (1<nargout), varargout{2}=Ksd; end
if (2<nargout), varargout{3}=Ksn; end
if (3<nargout), varargout{4}=Kds; end
if (4<nargout), varargout{5}=Kdd; end
if (5<nargout), varargout{6}=Kdn; end
if (6<nargout), varargout{7}=SS; end
if (7<nargout), varargout{8}=DS; end

textprogressbar(100);
textprogressbar('');

end

function S = hookesLaw(E,G,nu)
lambda=G*2*nu/(1-2*nu);
S.xx=2*G*E.xx+lambda*(E.xx+E.yy+E.zz);
S.yy=2*G*E.yy+lambda*(E.xx+E.yy+E.zz);
S.zz=2*G*E.zz+lambda*(E.xx+E.yy+E.zz);
S.xy=2*G*E.xy;
S.xz=2*G*E.xz;
S.yz=2*G*E.yz;
end
