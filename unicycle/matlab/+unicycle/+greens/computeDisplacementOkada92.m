function [ux,uy,uz]=computeDisplacementOkada92(X,Y,Z,G,nu,depth,angle,L,W,type,strike)
% COMPUTESTRESSOKADA92 calculates the displacement due to shear or tensile
% faults in a half space using the analytical solution of Okada (1992).
%
% INPUT:
%   x, y, z  - coordinate of observation point. units: m. column vector; (z negative down)
%   G        - rigidity
%   nu       - Poisson's ratio
%   depth    - depth of reference point. units: m
%              same point as Relax and EDCMP software (depth positive down)
%   angle    - dip angle (radians)
%   L        - fault length along strike. units: m
%   W        - fault width along dip. units: m
%   type     - fault type ('S':Strike slip; 'D':Dip slip; 'T':Tensile)
%
% OUTPUT:
% [ux,uy,uz]=computeStressOkada92(x,y,z,depth,angle,L,W,type)
% returns the displacement components at the observing points. 
%
% SEE ALSO: unicycle

assert(0<=depth,'unicycle.greens.computeDisplacementOkada92: depth must be positive.');
assert(0>=max(Z),'unicycle.greens.computeDisplacementOkada92: observation z-component must be negative.');

L1=-L/2;
L2=+L/2;
W1=-W/2;
W2=+W/2;

% description of fault patch is relative to the UPPER-LEFT corner of the
% fault, following the convention of Wang or the Relax series.
% Internal computations use the fault center coordinates.
% s=[sin(strike) cos(strike) 0];
% d=[cos(strike)*cos(angle) -sin(strike)*cos(angle) sin(angle)];
depth=depth+W*sin(angle)/2; % from top to center of the fault (depth is negative).
X=X-W/2*cos(angle)*cos(strike)-sin(strike)*L/2;
Y=Y+W/2*cos(angle)*sin(strike)-cos(strike)*L/2;

global critical lama lamu lamequot;

% lama and lamu are the Lamé and the shear modulus, respectively. units: Pa.
% critical: shortest distance to the points with singularities, this value is prescribed to
% avoid singularities in the calculation

WARNING('OFF');

% shear modulus
lamu=G;
% lamé parameter
lama=lamu*2*nu/(1-2*nu);
% Poisson's ratio
lamequot=(lama+lamu)./(lama+2.*lamu);

strike=-strike+pi/2;
x =  X*cos(strike)+Y*sin(strike);
y = -X*sin(strike)+Y*cos(strike);
z = Z;

critical=1e-6;

depth=-depth;
p0=y.*cos(angle)+(-z-depth).*sin(angle);        % Vector // to fault surface
p1=y.*cos(angle)+( z-depth).*sin(angle);        % Vector // to IMAGE fault surface

% Calculate the displacement
UA = FA(x,y,z,depth,angle,x-L1,p0-W1,type)-FA(x,y,z,depth,angle,x-L1,p0-W2,type)...
    -FA(x,y,z,depth,angle,x-L2,p0-W1,type)+FA(x,y,z,depth,angle,x-L2,p0-W2,type);
UAIMAGE = FA(x,y,-z,depth,angle,x-L1,p1-W1,type)-FA(x,y,-z,depth,angle,x-L1,p1-W2,type)...
    -FA(x,y,-z,depth,angle,x-L2,p1-W1,type)+FA(x,y,-z,depth,angle,x-L2,p1-W2,type);
UB = FB(x,y,z,depth,angle,x-L1,p0-W1,type)-FB(x,y,z,depth,angle,x-L1,p0-W2,type)...
    -FB(x,y,z,depth,angle,x-L2,p0-W1,type)+FB(x,y,z,depth,angle,x-L2,p0-W2,type);
UC = FC(x,y,z,depth,angle,x-L1,p0-W1,type)-FC(x,y,z,depth,angle,x-L1,p0-W2,type)...
    -FC(x,y,z,depth,angle,x-L2,p0-W1,type)+FC(x,y,z,depth,angle,x-L2,p0-W2,type);

Ux = 1./(2.*pi).*( UA{1}- UAIMAGE{1}+UB{1}+ z.*UC{1});
Uy = 1./(2.*pi).*((UA{2}-UAIMAGE{2}+UB{2}+z.*UC{2}).*cos(angle) ...
    -(UA{3}-UAIMAGE{3}+UB{3}+z.*UC{3}).*sin(angle));
uz = 1./(2.*pi).*((UA{2}-UAIMAGE{2}+UB{2}-z.*UC{2}).*sin(angle) ...
    +(UA{3}-UAIMAGE{3}+UB{3}-z.*UC{3}).*cos(angle));

% rotate in the original reference system
ux=-Uy*sin(strike)+Ux*cos(strike);
uy= Ux*sin(strike)+Uy*cos(strike);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions called for calculating the displacement and its derivative for a rectangular source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displacement caused by infinite medium terms
function FA=FA(x,y,z,depth,angle,epson,inta,type)
global critical lamequot;
WARNING('OFF');

[~,~,q,R,~,~,~,X11,~,~,Y11,~,~,~,~,~,~,~] ...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[theta,~,~,~,~,~]=DISPPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FA{1}=theta./2+lamequot./2.*epson.*q.*Y11;
        FA{2}=lamequot./2.*q./R;
        FA{3}=(1-lamequot)./2.*log(R+inta)-lamequot./2.*q.^2.*Y11;
        %erase singularities
        m=find(abs(R+inta)<critical);
        FA(2.*length(R(:))+m)=(1-lamequot)./2.*(-log(R(m)-inta(m)))-lamequot./2.*q(m).^2.*Y11(m);
    case 'D'
        FA{1}=lamequot./2.*q./R;
        FA{2}=theta./2+lamequot./2.*inta.*q.*X11;
        FA{3}=(1-lamequot)./2.*log(R+epson)-lamequot./2.*q.^2.*X11;
        % erase singularity
        m=find(abs(R+epson)<critical);
        FA(2.*length(R(:))+m)=(1-lamequot)./2.*(-log(R(m)-epson(m)))-lamequot./2.*q(m).^2.*X11(m);
    case 'T'
        FA{1}=-(1-lamequot)./2.*log(R+inta) -lamequot./2.*q.^2.*Y11;
        FA{2}=-(1-lamequot)./2.*log(R+epson)-lamequot./2.*q.^2.*X11;
        FA{3}=theta./2-lamequot./2.*q.*(inta.*X11+epson.*Y11);
        %erase singularity
        m=find(abs(R+inta)<critical);
        FA(m)=-(1-lamequot)./2.*(-log(R(m)-inta(m))) -lamequot./2.*q(m).^2.*Y11(m);
        m=find(abs(R+epson)<critical);
        FA(length(R(:))+m)=-(1-lamequot)./2.*(-log(R(m)-epson(m)))-lamequot./2.*q(m).^2.*X11(m);
    otherwise
        error('no such type');
end

end

% Displacement caused by surface deformation related term
function FB=FB(x,y,z,depth,angle,epson,inta,type)
global lamequot;
WARNING('OFF');
[~,~,q,R,yg,dg,~,X11,~,~,Y11,~,~,~,~,~,~,~]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[theta,~,I4,I3,I2,I1]=DISPPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FB{1}=-epson.*q.*Y11-theta-(1-lamequot)./lamequot.*I1.*sin(angle);
        FB{2}=-q./R+(1-lamequot)./lamequot.*yg./(R+dg).*sin(angle);
        FB{3}=q.^2.*Y11-(1-lamequot)./lamequot.*I2.*sin(angle);
    case 'D'
        FB{1}=-q./R+(1-lamequot)./lamequot.*I3.*sin(angle).*cos(angle);
        FB{2}=-inta.*q.*X11-theta-(1-lamequot)./lamequot.*epson./(R+dg).*sin(angle).*cos(angle);
        FB{3}= q.^2.*X11+(1-lamequot)./lamequot.*I4.*sin(angle).*cos(angle);
    case 'T'
        FB{1}=q.^2.*Y11-(1-lamequot)./lamequot.*I3.*(sin(angle)).^2;
        FB{2}=q.^2.*X11+(1-lamequot)./lamequot.*epson./(R+dg).*(sin(angle)).^2;
        FB{3}=q.*(inta.*X11+epson.*Y11)-theta-(1-lamequot)./lamequot.*I4.*(sin(angle)).^2;
    otherwise
        error('no such type');
end

end

% Displacement caused by depth dependent term
function FC=FC(x,y,z,depth,angle,epson,inta,type)
global lamequot;
WARNING('OFF');
[~,~,q,R,yg,dg,cg,X11,X32,~,Y11,~,~,~,Z32,~,~,~]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FC{1}=(1-lamequot).*epson.*Y11.*cos(angle)-lamequot.*epson.*q.*Z32;
        FC{2}=(1-lamequot).*(cos(angle)./R+2.*q.*Y11.*sin(angle))-lamequot.*cg.*q./R.^3;
        FC{3}=(1-lamequot).*q.*Y11.*cos(angle)-lamequot.*(cg.*inta./R.^3-z.*Y11+epson.^2.*Z32);
    case 'D'
        FC{1}=(1-lamequot).*cos(angle)./R-q.*Y11.*sin(angle)-lamequot.*cg.*q./R.^3;
        FC{2}=(1-lamequot).*yg.*X11-lamequot.*cg.*inta.*q.*X32;
        FC{3}=-dg.*X11-epson.*Y11.*sin(angle)-lamequot.*cg.*(X11-q.^2.*X32);
    case 'T'
        FC{1}=-(1-lamequot).*(sin(angle)./R+q.*Y11.*cos(angle))-lamequot.*(z.*Y11-q.^2.*Z32);
        FC{2}=(1-lamequot).*2.*epson.*Y11.*sin(angle)+dg.*X11-lamequot.*cg.*(X11-q.^2.*X32);
        FC{3}=(1-lamequot).*(yg.*X11+epson.*Y11.*cos(angle))+lamequot.*q.*(cg.*inta.*X32+epson.*Z32);
    otherwise
        error('no such type');
end

end


%----------------------------------------------------------------------------------------------------
% Functions for calculating internal deformation (see Okada, 1992, Table 2-9);
% The output variables are theta,X,I4,I3,I2,I1;
%----------------------------------------------------------------------------------------------------
function [theta,X,I4,I3,I2,I1]=DISPPARA(x,y,z,depth,angle,epson,inta)
global critical;
WARNING('OFF');

[~,~,q,R,yg,dg,~,~,~,~,~,~,~,~,~,~,~,~]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
theta=atan(epson.*inta./(q.*R));
% erase theta's singularity occurring on the planes that include the fault surface and its image;
m=find(abs(q)<critical);
theta(m)=0;
X=(epson.^2+q.^2).^0.5;
if abs(cos(angle))<critical
    I3=0.5.*(inta./(R+dg)+yg.*q./(R+dg).^2-log(R+inta));
    I4=0.5.*epson.*yg./(R+dg).^2;
    % erase I3's singularity occurring along the lines extending the edges that are
    % perpendicular to the fault strike and p<0;
    m=find(abs(R+inta)<critical);
    I3(m)=0.5.*(inta(m)./(R(m)+dg(m))+yg(m).*q(m)./(R(m)+dg(m)).^2+log(R(m)-inta(m)));
else
    I4=tan(angle).*epson./(R+dg)+2./(cos(angle)).^2.*atan((inta.*(X+q.*cos(angle))...
        +X.*(R+X).*sin(angle))./(epson.*(R+X).*cos(angle)));
    I3=1./cos(angle).*yg./(R+dg)-1./(cos(angle)).^2.*(log(R+inta)-sin(angle).*log(R+dg));
    % erase I3's singularity occurring along the lines extending the edges that are
    % perpendicular to the fault strike and p<0;
    m=find(abs(R+inta)<critical);
    I3(m)=1./cos(angle).*yg(m)./(R(m)+dg(m))-1./(cos(angle)).^2.*(-log(R(m)-inta(m))...
        -sin(angle).*log(R(m)+dg(m)));
end
% erase I4's singularity occurring on the vertical planes that include the edges that are
% perpendicular to the fault strike;
m=find(abs(epson)<critical);
I4(m)=0;
I1=-epson./(R+dg).*cos(angle)-I4.*sin(angle);
I2=log(R+dg)+I3.*sin(angle);
end

%----------------------------------------------------------------------------------------------------
% common variables;
% The output variables are
% d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0;
%----------------------------------------------------------------------------------------------------
function [d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta)
global critical;
WARNING('OFF');
d=-z-depth;
p=y.*cos(angle)+d.*sin(angle);
q=y.*sin(angle)-d.*cos(angle);
R=(epson.^2+inta.^2+q.^2).^0.5;
% erase R's singularity
m=find(abs(R)<critical);
R(m)=NaN;
yg=inta.*cos(angle)+q.*sin(angle);
dg=inta.*sin(angle)-q.*cos(angle);
cg=dg+z;

X11=1./(R.*(R+epson));
X32=(2.*R+epson)./(R.^3.*(R+epson).^2);
X53=(8.*R.^2+9.*R.*epson+3.*epson.^2)./(R.^5.*(R+epson).^3);
%erase X11,X53 and X32's singularities
m=find(abs(R+epson)<critical);
X11(m)=0;
X32(m)=0;
X53(m)=0;

Y11=1./(R.*(R+inta));
Y32=(2.*R+inta)./(R.^3.*(R+inta).^2);
Y53=(8.*R.^2+9.*R.*inta+3.*inta.^2)./(R.^5.*(R+inta).^3);
%erase Y11 singularities
m=find(abs(R+inta)<critical);
Y11(m)=0;
Y32(m)=0;
Y53(m)=0;

h=q.*cos(angle)-z;
Z32=sin(angle)./R.^3-h.*Y32;
Z53=3.*sin(angle)./R.^5-h.*Y53;
Y0=Y11-epson.^2.*Y32;
Z0=Z32-epson.^2.*Z53;
end

