function [K1111,K1211,K1311,K2211,K2311,K3311, ...
          K1112,K1212,K1312,K2212,K2312,K3312, ...
          K1113,K1213,K1313,K2213,K2313,K3313, ...
          K1122,K1222,K1322,K2222,K2322,K3322, ...
          K1123,K1223,K1323,K2223,K2323,K3323, ...
          K1133,K1233,K1333,K2233,K2333,K3333]=computeStressKernelsVerticalShearZone(shz,rcv,G,nu)
% computeStressKernelsVerticalShearZone computes the stress on a shear zone due to
% strain on a shear zone.
%
%   [K1111,K1112,K1113,K1122,K1123,K1133, ...
%    K1211,K1212,K1213,K1222,K1223,K1233, ...
%    K1311,K1312,K1313,K1322,K1323,K1333, ...
%    K2211,K2212,K2213,K2222,K2223,K2233, ...
%    K2311,K2312,K2313,K2322,K2323,K2333, ...
%    K3311,K3312,K3313,K3322,K3323,K3333]=COMPUTESTRESSKERNELSVERTICALSHEARZONE(shz,rcv,G,nu)
%
% where nu is Poisson's ratio, G is the shear modulus. shz is a structure
% with shear zone geometry:
%
%   shz.x         : position of vertices
%   shz.L         : array of length of each shear zone
%   shz.W         : array of width of each shear zone
%   shz.T         : array of thickness of each shear zone
%   shz.strike    : array of strike of each shear zone
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
% Kijkl is the stress component Skl due to shear Eij, both in the shear 
% zone centric (primed) system of coordinates.
%
% Sylvain Barbot and James D. P. Moore, 07/16/2016
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle

import unicycle.greens.*
import unicycle.utils.*

N=length(shz.L);
M=size(rcv.xc,1);

lambda=G*2*nu/(1-2*nu);

out=cell(6,6);

textprogressbar('# shear zone stress kernels: ');

e=eye(6);

% loop of causes (strain on shear zones)
for j=1:6
    
    % initialize stress kernels
    K11=zeros(M,N);
    K12=zeros(M,N);
    K13=zeros(M,N);
    K22=zeros(M,N);
    K23=zeros(M,N);
    K33=zeros(M,N);
    
    for k=1:N
        
        if 0==mod(k-1,2)
            textprogressbar(((j-1)*N+k)/(6*N)*100);
        end
        
        % stress components in Matlab x,y,z system of coordinates
        [s22,s12,s23,s11,s13,s33]=computeStressVerticalShearZone(rcv.xc(:,2),rcv.xc(:,1),-rcv.xc(:,3), ...
            shz.x(k,2),shz.x(k,1),-shz.x(k,3),shz.L(k),shz.T(k),shz.W(k),shz.strike(k), ...
            e(1,j),e(2,j),e(3,j),e(4,j),e(5,j),e(6,j), ...
            G,nu);
        s13=-s13;
        s23=-s23;
        
        % rotate the stress in the receiver direction
        % e_i' = A_ij e_j
        % sigma_ij' e_i' e_j' = sigma_ij' Aik e_k Ajl e_l
        % sigma_kl = Aik sigma_ij' Ajl
        % sigma = A' sigma' A
        % sigma' = A sigma A'
        % A = [shz.sv(1,:);shz.nv(1,:);-shz.dv(1,:);]
        
        s11p= ( rcv.sv(:,1).*s11+rcv.sv(:,2).*s12+rcv.sv(:,3).*s13).*rcv.sv(:,1)+( rcv.sv(:,1).*s12+rcv.sv(:,2).*s22+rcv.sv(:,3).*s23).*rcv.sv(:,2)+( rcv.sv(:,1).*s13+rcv.sv(:,2).*s23+rcv.sv(:,3).*s33).*rcv.sv(:,3);
        s12p= ( rcv.sv(:,1).*s11+rcv.sv(:,2).*s12+rcv.sv(:,3).*s13).*rcv.nv(:,1)+( rcv.sv(:,1).*s12+rcv.sv(:,2).*s22+rcv.sv(:,3).*s23).*rcv.nv(:,2)+( rcv.sv(:,1).*s13+rcv.sv(:,2).*s23+rcv.sv(:,3).*s33).*rcv.nv(:,3);
        s13p=-( rcv.sv(:,1).*s11+rcv.sv(:,2).*s12+rcv.sv(:,3).*s13).*rcv.dv(:,1)-( rcv.sv(:,1).*s12+rcv.sv(:,2).*s22+rcv.sv(:,3).*s23).*rcv.dv(:,2)-( rcv.sv(:,1).*s13+rcv.sv(:,2).*s23+rcv.sv(:,3).*s33).*rcv.dv(:,3);
        s22p= ( rcv.nv(:,1).*s11+rcv.nv(:,2).*s12+rcv.nv(:,3).*s13).*rcv.nv(:,1)+( rcv.nv(:,1).*s12+rcv.nv(:,2).*s22+rcv.nv(:,3).*s23).*rcv.nv(:,2)+( rcv.nv(:,1).*s13+rcv.nv(:,2).*s23+rcv.nv(:,3).*s33).*rcv.nv(:,3);
        s23p=-( rcv.nv(:,1).*s11+rcv.nv(:,2).*s12+rcv.nv(:,3).*s13).*rcv.dv(:,1)-( rcv.nv(:,1).*s12+rcv.nv(:,2).*s22+rcv.nv(:,3).*s23).*rcv.dv(:,2)-( rcv.nv(:,1).*s13+rcv.nv(:,2).*s23+rcv.nv(:,3).*s33).*rcv.dv(:,3);
        s33p=-(-rcv.dv(:,1).*s11-rcv.dv(:,2).*s12-rcv.dv(:,3).*s13).*rcv.dv(:,1)-(-rcv.dv(:,1).*s12-rcv.dv(:,2).*s22-rcv.dv(:,3).*s23).*rcv.dv(:,2)-(-rcv.dv(:,1).*s13-rcv.dv(:,2).*s23-rcv.dv(:,3).*s33).*rcv.dv(:,3);
        
        % stress kernels in shear zone (primed) system of coordinates
        K11(:,k)=s11p(:);
        K12(:,k)=s12p(:);
        K13(:,k)=s13p(:);
        K22(:,k)=s22p(:);
        K23(:,k)=s23p(:);
        K33(:,k)=s33p(:);
        
    end
    
    out{j,1}=K11;
    out{j,2}=K12;
    out{j,3}=K13;
    out{j,4}=K22;
    out{j,5}=K23;
    out{j,6}=K33;
    
end

[K1111,K1211,K1311,K2211,K2311,K3311, ...
 K1112,K1212,K1312,K2212,K2312,K3312, ...
 K1113,K1213,K1313,K2213,K2313,K3313, ...
 K1122,K1222,K1322,K2222,K2322,K3322, ...
 K1123,K1223,K1323,K2223,K2323,K3323, ...
 K1133,K1233,K1333,K2233,K2333,K3333]=deal(out{:});

textprogressbar(100);
textprogressbar('');

end
