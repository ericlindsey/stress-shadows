function [shz,wn,lt]=shz2shz(xo,L,W,T,strike,dip,lo,wo,to,alphal,alphaw,alphat)
% function SHZ2SHZ subsamples a shear zone into smaller shear zones
% of length, width, thickness starting from lo, wo, and to, respectively
% increasing geometrically with down-dip distance with increment 
% alphal, alphaw, and alphat (alphal>1 for increase).
%
%   [shz]=shz2shz(xo,L,W,T,strike,dip,lo,wo,to,alphal,alphaw,alphat)
%
% input:
%   xo     origin position vector [x1 (north);x2 (east);x3 (down)]
%   L      total length
%   W      total width
%   T      total thickness
%   strike strike angle in degrees
%   dip    dip angle in degrees
%   lo     approximative initial length of output shear zones
%   wo     approximative initial width of output shear zones
%   to     approximative initial thickness of output shear zones
%   alpha1 geometric factor for length increase
%   alphaw geometric factor for width increase
%   alphat geometric factor for thickness increase
%
% output:
%   shz    list of shear zones in the format
%          x1,x2,x3,length,width,thickness,strike,dip
%
% to write an ascii output compatible with unicycle input:
%
%   fprintf('%3i %f %f %f %f %f %f %f %f\n',[[1:length(shz)]',shz]');

% create wi
Wc=W;
k=0;
w=0;
while Wc>0
    Wt=wo*alphaw^k;
    if Wt > Wc/2
        Wt = Wc;
    end
    wn=min([Wt,Wc]);
    w=[w; wn];
    k=k+1;
    Wc=Wc-wn;
end
Nw=k;

% strike and dip direction normal vectors
Sv=[ cosd(strike); sind(strike); 0];
Dv=[-cosd(dip)*sind(strike); cosd(dip)*cosd(strike); sind(dip)];
Nv=[-sind(strike)*sind(dip);+cosd(strike)*sind(dip);-cosd(dip)];

shz=[];

% loop in dip direction
for k=1:Nw
    lt=lo*alphal^(k-1);
    Nl=ceil(L/lt);
    lt=L/Nl;
    
    tt=to*alphat^(k-1);
    Nt=ceil(T/tt);
    tt=T/Nt;

    % look in normal direction
    for j=1:Nt
        % loop in strike direction
        for i=1:Nl
            x=xo+sum(w(1:k))*Dv+((j-0.5)*tt-T/2)*Nv+(i-1)*lt*Sv;
            shz=[shz; [x',lt,w(k+1),tt,strike,dip]]; 
        end
    end
end




        