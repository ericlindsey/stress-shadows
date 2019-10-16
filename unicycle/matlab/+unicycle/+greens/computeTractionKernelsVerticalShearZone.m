function [K11s,K12s,K13s,K22s,K23s,K33s, ...
          K11d,K12d,K13d,K22d,K23d,K33d, ...
          K11n,K12n,K13n,K22n,K23n,K33n]=computeTractionKernelsVerticalShearZone(shz,rcv,G,nu)
% computeTractionKernelsVerticalShearZone computes the stress on a receiver faults due to
% strain on a shear zone.
%
%   [K11s,K11d,K11n, ...
%    K12s,K12d,K12n, ...
%    K13s,K13d,K13n, ...
%    K22s,K22d,K22n, ...
%    K23s,K23d,K23n, ...
%    K33s,K33d,K33n]=COMPUTETRACTIONKERNELSVERTICALSHEARZONE(shz,rcv,G,nu)
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
%   [Kss,Ksd,Ksn,Kds,Kdd,Kdn]=TRIANGLESTRESSKERNELS(shz,rcv,G,nu)
%
% returns the shear stress kernels. Kss is the traction in the strike
% direction due to strike slip. Ksd is the traction in the dip direction
% due to strike slip. Kds is the tarction in the strike direction due to
% dip slip and Kdd is the traction in the dip direction due to dip slip.
%
%
% Sylvain Barbot, 07/16/2016
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle

import unicycle.greens.*
import unicycle.utils.*

N=length(shz.L);
M=size(rcv.xc,1);

lambda=G*2*nu/(1-2*nu);

out=cell(6,3);

textprogressbar('# shear zone traction kernels: ');

e=eye(6);

% loop of causes (strain on shear zones)
for j=1:6
    
    % initialize stress kernels
    Kes=zeros(M,N);
    Ked=zeros(M,N);
    Ken=zeros(M,N);
    
    for k=1:N
        
        if 0==mod(k-1,2)
            textprogressbar(((j-1)*N+k)/(6*N)*100);
        end
        
        % stress components
        [s22,s12,s23,s11,s13,s33]=computeStressVerticalShearZone(rcv.xc(:,2),rcv.xc(:,1),-rcv.xc(:,3), ...
            shz.x(k,2),shz.x(k,1),-shz.x(k,3),shz.L(k),shz.T(k),shz.W(k),shz.strike(k), ...
            e(1,j),e(2,j),e(3,j),e(4,j),e(5,j),e(6,j), ...
            G,nu);
        s13=-s13;
        s23=-s23;
        
        % full traction vector on receiver faults
        t=[...
            s11.*rcv.nv(:,1)+s12.*rcv.nv(:,2)+s13.*rcv.nv(:,3), ...
            s12.*rcv.nv(:,1)+s22.*rcv.nv(:,2)+s23.*rcv.nv(:,3), ...
            s13.*rcv.nv(:,1)+s23.*rcv.nv(:,2)+s33.*rcv.nv(:,3)];
        
        % shear stress in strike direction
        Kes(:,k)=sum(t.*rcv.sv,2);
        
        % shear stress in dip direction
        Ked(:,k)=sum(t.*rcv.dv,2);
        
        % stress in normal direction
        Ken(:,k)=sum(t.*rcv.nv,2);
        
    end
    
    out{j,1}=Kes;
    out{j,2}=Ked;
    out{j,3}=Ken;
    
end

[K11s,K12s,K13s,K22s,K23s,K33s, ...
 K11d,K12d,K13d,K22d,K23d,K33d, ...
 K11n,K12n,K13n,K22n,K23n,K33n]=deal(out{:});

textprogressbar(100);
textprogressbar('');

end
