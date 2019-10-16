function [varargout]=computeOkada92(U,X,Y,Z,G,nu,depth,angle,L,W,type,gridflag,strike)
% COMPUTEOKADA92 calculates internal deformation due to shear or tensile
% faults in a half space using the analytical solution of Okada (1992).
%
% INPUT:
%
%   U        - fault dislocation (slip). units: km;
%   x, y, z  - coordinate of observation point. units: km. column vector; (z negative down)
%   nu       - Poisson's ratio
%   depth    - depth of reference point. units: km; 
%              same point as Relax and EDCMP software (depth positive down!)
%   angle    - dip angle (radians);
%   L        - fault length along strike. units: km;
%   W        - fault width along dip. units: km;
%   type     - fault type ('S':Strike slip; 'D':Dip slip; 'T':Tensile)
%
% R. Kanda added the last input argument so that this function can be run multiple times without command line input.
% gridflag: if this flag is 0, the observation points are NOT GRIDIZED, otherwise, MESHGRID is used to compute the velocity
% field in a rectangular volume encompassed by [(xmin, ymin, zmin), (xmax, ymax, zmax)].  This takes a lot of computational effort!
%
% OUTPUT:
%
% [displacement(km), derivative, strain, stress(bar)].
%
% M=computeOkada92(U, x, y, z, depth, angle, L, W, type, gridflag) 
% returns displacements at the observing point. M is a 4D array with 
% length LENGTH(x)* LENGTH(y)* LENGTH(z)*3.
% M(:, :, :, i) denotes the corresponding displacement component.
% If gridflag = 1, use meshgrid to "grid" the observation points.
%
% [M, N]=computeOkada92(U, x, y, z, depth, angle, L, W, type, gridflag) 
% returns displacements and spatial derivatives at the observation point. 
% M is the same as above. N is a LENGTH(x)*LENGTH(y)*LENGTH(z)*3*3 5D array. 
% N (:, :, :, i, j) denotes the spatial derivative of the i-th displacement 
% component with respect to the j-th axis-direction at the observing points.
%
% [M, N, P]=computeOkada92(U, x, y, z, depth, angle, L, W, type, gridflag) 
% return displacements and spatial derivatives and strain at the observing 
% point. M and N are the same as above. 
% P(:, :, :, i, j) denotes the strain tensor at the observing point;
%
% [M, N, P, Q]=computeOkada92(U, x, y, z, depth, angle, L, W, type, gridflag)
% return displacements and spatial derivatives and strains at the observing
% point. M, N and P are the same as above. 
% Q(:, :, :, i, j) denotes the stress tensor at the observing point;
%
% Modified by Ravi Kanda: October, 2006.
% Modified by S. Barbot: October, 2010: fix a missing code line.
%
% SEE ALSO: unicycle

assert(0<=depth,'unicycle.greens.computeOkada92: depth must be positive.');
assert(0>=max(Z),'unicycle.greens.computeOkada92: observation z-component must be negative.');

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
warning off;

% shear modulus
lamu=G;
% lamé parameter
lama=lamu*2*nu/(1-2*nu);
% Poisson's ratio
lamequot=(lama+lamu)./(lama+2.*lamu);

strike=-strike+pi/2;
rotx =  X*cos(strike)+Y*sin(strike);
roty = -X*sin(strike)+Y*cos(strike);
X=rotx;
Y=roty;

critical=1e-6;
% check input parameters
if ndims(X)<3
    if gridflag==1
        [x,y,z]=meshgrid(X,Y,Z);
        disp(['In Elastipatch: Shape of input x,y,z vectors = (',int2str(size(X)),'), (',int2str(size(Y)),'), (',int2str(size(Z)),')'])
        disp(['In Elastipatch: Shape of Grid X,Y,Z coordinates = (',int2str(size(x)),'), (',int2str(size(y)),'), (',int2str(size(z)),')'])
    else
        if (length(X)~=length(Y))|(length(X)~=length(Z))
            error('*** If you do not gridize, vectors x, y, & z MUST have the same size!');
        end
        x(1,1,:)=X;
        y(1,1,:)=Y;
        z(1,1,:)=Z;
    end
end

depth=-depth;
p0=y.*cos(angle)+(-z-depth).*sin(angle);        % Vector // to fault surface
p1=y.*cos(angle)+(z-depth).*sin(angle);         % Vector // to IMAGE fault surface

% Calculate the displacement/varargout{1};
UA = FA(x,y,z,depth,angle,x-L1,p0-W1,type)-FA(x,y,z,depth,angle,x-L1,p0-W2,type)...
    -FA(x,y,z,depth,angle,x-L2,p0-W1,type)+FA(x,y,z,depth,angle,x-L2,p0-W2,type);
UAIMAGE = FA(x,y,-z,depth,angle,x-L1,p1-W1,type)-FA(x,y,-z,depth,angle,x-L1,p1-W2,type)...
    -FA(x,y,-z,depth,angle,x-L2,p1-W1,type)+FA(x,y,-z,depth,angle,x-L2,p1-W2,type);
UB = FB(x,y,z,depth,angle,x-L1,p0-W1,type)-FB(x,y,z,depth,angle,x-L1,p0-W2,type)...
    -FB(x,y,z,depth,angle,x-L2,p0-W1,type)+FB(x,y,z,depth,angle,x-L2,p0-W2,type);
UC = FC(x,y,z,depth,angle,x-L1,p0-W1,type)-FC(x,y,z,depth,angle,x-L1,p0-W2,type)...
    -FC(x,y,z,depth,angle,x-L2,p0-W1,type)+FC(x,y,z,depth,angle,x-L2,p0-W2,type);

varargout{1}(:,:,:,1) = U./(2.*pi).*( UA(:,:,:,1)- UAIMAGE(:,:,:,1)+UB(:,:,:,1)+ z.*UC(:,:,:,1));
varargout{1}(:,:,:,2) = U./(2.*pi).*((UA(:,:,:,2)-UAIMAGE(:,:,:,2)+UB(:,:,:,2)+z.*UC(:,:,:,2)).*cos(angle) ...
    -(UA(:,:,:,3)-UAIMAGE(:,:,:,3)+UB(:,:,:,3)+z.*UC(:,:,:,3)).*sin(angle));
varargout{1}(:,:,:,3) = U./(2.*pi).*((UA(:,:,:,2)-UAIMAGE(:,:,:,2)+UB(:,:,:,2)-z.*UC(:,:,:,2)).*sin(angle) ...
    +(UA(:,:,:,3)-UAIMAGE(:,:,:,3)+UB(:,:,:,3)-z.*UC(:,:,:,3)).*cos(angle));

% rotate in the original reference system
ux=-varargout{1}(:,:,:,2)*sin(strike)+varargout{1}(:,:,:,1)*cos(strike);
uy=+varargout{1}(:,:,:,1)*sin(strike)+varargout{1}(:,:,:,2)*cos(strike);
varargout{1}(:,:,:,1)=ux;
varargout{1}(:,:,:,2)=uy;

if nargout>=2
    % Calculate the derivative/varargout{2};
    UAD = FAD(x,y,z,depth,angle,x-L1,p0-W1,type)-FAD(x,y,z,depth,angle,x-L1,p0-W2,type)...
        -FAD(x,y,z,depth,angle,x-L2,p0-W1,type)+FAD(x,y,z,depth,angle,x-L2,p0-W2,type);
    UADIMAGE = FAD(x,y,-z,depth,angle,x-L1,p1-W1,type)-FAD(x,y,-z,depth,angle,x-L1,p1-W2,type)...
        -FAD(x,y,-z,depth,angle,x-L2,p1-W1,type)+FAD(x,y,-z,depth,angle,x-L2,p1-W2,type);
    UBD = FBD(x,y,z,depth,angle,x-L1,p0-W1,type)-FBD(x,y,z,depth,angle,x-L1,p0-W2,type)...
        -FBD(x,y,z,depth,angle,x-L2,p0-W1,type)+FBD(x,y,z,depth,angle,x-L2,p0-W2,type);
    UCD = FCD(x,y,z,depth,angle,x-L1,p0-W1,type)-FCD(x,y,z,depth,angle,x-L1,p0-W2,type)...
        -FCD(x,y,z,depth,angle,x-L2,p0-W1,type)+FCD(x,y,z,depth,angle,x-L2,p0-W2,type);
    for i=1:2,
        varargout{2}(:,:,:,1,i) = U./(2.*pi).*( UAD(:,:,:,1,i)-UADIMAGE(:,:,:,1,i)+UBD(:,:,:,1,i)+z.*UCD(:,:,:,1,i));
        varargout{2}(:,:,:,2,i) = ...
            U./(2.*pi).*((UAD(:,:,:,2,i)-UADIMAGE(:,:,:,2,i)+UBD(:,:,:,2,i)+z.*UCD(:,:,:,2,i)).*cos(angle)...
            -(UAD(:,:,:,3,i)-UADIMAGE(:,:,:,3,i)+UBD(:,:,:,3,i)+z.*UCD(:,:,:,3,i)).*sin(angle));
        varargout{2}(:,:,:,3,i) = ...
            U./(2.*pi).*((UAD(:,:,:,2,i)-UADIMAGE(:,:,:,2,i)+UBD(:,:,:,2,i)-z.*UCD(:,:,:,2,i)).*sin(angle)...
            +(UAD(:,:,:,3,i)-UADIMAGE(:,:,:,3,i)+UBD(:,:,:,3,i)-z.*UCD(:,:,:,3,i)).*cos(angle));
    end
    varargout{2}(:,:,:,1,3) = U./(2.*pi).*( UAD(:,:,:,1,3)+UADIMAGE(:,:,:,1,3)+UBD(:,:,:,1,3)+UC(:,:,:,1)+z.*UCD(:,:,:,1,3));
    varargout{2}(:,:,:,2,3) = U./(2.*pi).*((UAD(:,:,:,2,3)+UADIMAGE(:,:,:,2,3)+UBD(:,:,:,2,3)+UC(:,:,:,2)+z.*UCD(:,:,:,2,3)).*cos(angle)-(UAD(:,:,:,3,3)+UADIMAGE(:,:,:,3,3)+UBD(:,:,:,3,3)+UC(:,:,:,3)+z.*UCD(:,:,:,3,3)).*sin(angle));
    varargout{2}(:,:,:,3,3) = U./(2.*pi).*((UAD(:,:,:,2,3)+UADIMAGE(:,:,:,2,3)+UBD(:,:,:,2,3)-UC(:,:,:,2)...
        -z.*UCD(:,:,:,2,3)).*sin(angle)+(UAD(:,:,:,3,3)+UADIMAGE(:,:,:,3,3)+UBD(:,:,:,3,3)-UC(:,:,:,3)-z.*UCD(:,:,:,3,3)).*cos(angle));
    
    % rotate in the original reference system
    %strike=-pi/2+strike;
    %strike=-strike;
    UXX=(cos(strike)*varargout{2}(:,:,:,1,1)-sin(strike)*varargout{2}(:,:,:,2,1))*cos(strike)...
        -(cos(strike)*varargout{2}(:,:,:,1,2)-sin(strike)*varargout{2}(:,:,:,2,2))*sin(strike);
    UXY=(cos(strike)*varargout{2}(:,:,:,1,1)-sin(strike)*varargout{2}(:,:,:,2,1))*sin(strike)...
        +(cos(strike)*varargout{2}(:,:,:,1,2)-sin(strike)*varargout{2}(:,:,:,2,2))*cos(strike);
    UXZ=cos(strike)*varargout{2}(:,:,:,1,3)...
        -sin(strike)*varargout{2}(:,:,:,2,3);
    UYX=(sin(strike)*varargout{2}(:,:,:,1,1)+cos(strike)*varargout{2}(:,:,:,2,1))*cos(strike)...
        -(sin(strike)*varargout{2}(:,:,:,1,2)+cos(strike)*varargout{2}(:,:,:,2,2))*sin(strike);
    UYY=(sin(strike)*varargout{2}(:,:,:,1,1)+cos(strike)*varargout{2}(:,:,:,2,1))*sin(strike)...
        +(sin(strike)*varargout{2}(:,:,:,1,2)+cos(strike)*varargout{2}(:,:,:,2,2))*cos(strike);
    UYZ=sin(strike)*varargout{2}(:,:,:,1,3)+cos(strike)*varargout{2}(:,:,:,2,3);
    UZX=varargout{2}(:,:,:,3,1)*cos(strike)-varargout{2}(:,:,:,3,2)*sin(strike);
    UZY=varargout{2}(:,:,:,3,1)*sin(strike)+varargout{2}(:,:,:,3,2)*cos(strike);
    
    varargout{2}(:,:,:,1,1)=UXX;
    varargout{2}(:,:,:,1,2)=UXY;
    varargout{2}(:,:,:,1,3)=UXZ;
    varargout{2}(:,:,:,2,1)=UYX;
    varargout{2}(:,:,:,2,2)=UYY;
    varargout{2}(:,:,:,2,3)=UYZ;
    varargout{2}(:,:,:,3,1)=UZX;
    varargout{2}(:,:,:,3,2)=UZY;
    
    
end
%disp(['In Elastipatch: Shape of D = ',int2str(size(varargout{2}))])

if nargout>=3
    % Calculate the strain/varargout{3}
    for i=1:3,
        for j=1:i,
            varargout{3}(:,:,:,i,j)=0.5.*(varargout{2}(:,:,:,i,j)+varargout{2}(:,:,:,j,i));
            if i~=j
                varargout{3}(:,:,:,j,i)=varargout{3}(:,:,:,i,j);
            end
        end
    end
end
%disp(['In Elastipatch: Shape of E = ',int2str(size(varargout{3}))])

if nargout==4
    % Calculate the stress/varargout{4}
    I=[1 0 0;0 1 0;0 0 1];
    for i=1:3,
        for j=1:i,
            varargout{4}(:,:,:,i,j)=lama.*(varargout{3}(:,:,:,1,1)+varargout{3}(:,:,:,2,2)+varargout{3}(:,:,:,3,3)).*I(i,j)...
                +2.*lamu.*varargout{3}(:,:,:,i,j);
            if i~=j
                varargout{4}(:,:,:,j,i)=varargout{4}(:,:,:,i,j);
            end
        end
    end
end
%disp(['In Elastipatch: Shape of S = ',int2str(size(varargout{4}))])

end     % END MAIN PROGRAM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions called for calculating the displacement and its derivative for a rectangular source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displacement caused by infinite medium terms
function FA=FA(x,y,z,depth,angle,epson,inta,type)
global critical lamequot;
warning off;

[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0] ...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[theta,X,I4,I3,I2,I1]=DISPPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FA(:,:,:,1)=theta./2+lamequot./2.*epson.*q.*Y11;
        FA(:,:,:,2)=lamequot./2.*q./R;
        FA(:,:,:,3)=(1-lamequot)./2.*log(R+inta)-lamequot./2.*q.^2.*Y11;
        %erase singularities
        m=find(abs(R+inta)<critical);
        FA(2.*length(R(:))+m)=(1-lamequot)./2.*(-log(R(m)-inta(m)))-lamequot./2.*q(m).^2.*Y11(m);
    case 'D'
        FA(:,:,:,1)=lamequot./2.*q./R;
        FA(:,:,:,2)=theta./2+lamequot./2.*inta.*q.*X11;
        FA(:,:,:,3)=(1-lamequot)./2.*log(R+epson)-lamequot./2.*q.^2.*X11;
        % erase singularity
        m=find(abs(R+epson)<critical);
        FA(2.*length(R(:))+m)=(1-lamequot)./2.*(-log(R(m)-epson(m)))-lamequot./2.*q(m).^2.*X11(m);
    case 'T'
        FA(:,:,:,1)=-(1-lamequot)./2.*log(R+inta) -lamequot./2.*q.^2.*Y11;
        FA(:,:,:,2)=-(1-lamequot)./2.*log(R+epson)-lamequot./2.*q.^2.*X11;
        FA(:,:,:,3)=theta./2-lamequot./2.*q.*(inta.*X11+epson.*Y11);
        %erase singularity
        m=find(abs(R+inta)<critical);
        FA(m)=-(1-lamequot)./2.*(-log(R(m)-inta(m))) -lamequot./2.*q(m).^2.*Y11(m);
        m=find(abs(R+epson)<critical);
        FA(length(R(:))+m)=-(1-lamequot)./2.*(-log(R(m)-epson(m)))-lamequot./2.*q(m).^2.*X11(m);
    otherwise
        error('no such type');
end
%        fprintf('....');
end

% derivative caused by infinite medium terms
function FAD=FAD(x,y,z,depth,angle,epson,inta,type)
global critical lama lamu lamequot;
warning off;

[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[D11,K1,K2,K3,K4,J1,J2,J3,J4,J5,J6,E,F,G,H,P,Q,E1,F1,G1,H1,P1,Q1] ...
    =DERIPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FAD(:,:,:,1,1)=-(1-lamequot)./2.*q.*Y11-lamequot./2.*epson.^2.*q.*Y32;
        FAD(:,:,:,2,1)=-lamequot./2.*epson.*q./R.^3;
        FAD(:,:,:,3,1)= (1-lamequot)./2.*epson.*Y11+lamequot./2.*epson.*q.^2.*Y32;
        FAD(:,:,:,1,2)=(1-lamequot)./2.*epson.*Y11.*sin(angle)+dg./2.*X11 +lamequot./2.*epson.*F;
        FAD(:,:,:,2,2)= lamequot./2.*E;
        FAD(:,:,:,3,2)= (1-lamequot)./2.*(cos(angle)./R+q.*Y11.*sin(angle))-lamequot./2.*q.*F;
        FAD(:,:,:,1,3)=(1-lamequot)./2.*epson.*Y11.*cos(angle)+yg./2.*X11 +lamequot./2.*epson.*F1;
        FAD(:,:,:,2,3)= lamequot./2.*E1;
        FAD(:,:,:,3,3)=-(1-lamequot)./2.*(sin(angle)./R-q.*Y11.*cos(angle))-lamequot./2.*q.*F1;
    case 'D'
        FAD(:,:,:,1,1)=-lamequot./2.*epson.*q./R.^3;
        FAD(:,:,:,2,1)=-q./2.*Y11-lamequot./2.*inta.*q./R.^3;
        FAD(:,:,:,3,1)= (1-lamequot)./2./R+lamequot./2.*q.^2./R.^3;
        FAD(:,:,:,1,2)= lamequot./2.*E;
        FAD(:,:,:,2,2)=(1-lamequot)./2.*dg.*X11+epson./2.*Y11.*sin(angle)+lamequot./2.*inta.*G;
        FAD(:,:,:,3,2)=(1-lamequot)./2.*yg.*X11-lamequot./2.*q.*G;
        FAD(:,:,:,1,3)= lamequot./2.*E1;
        FAD(:,:,:,2,3)=(1-lamequot)./2.*yg.*X11+epson./2.*Y11.*cos(angle)+lamequot./2.*inta.*G1;
        FAD(:,:,:,3,3)=-(1-lamequot)./2.*dg.*X11-lamequot./2.*q.*G1;
    case 'T'
        FAD(:,:,:,1,1)=-(1-lamequot)./2.*epson.*Y11+lamequot./2.*epson.*q.^2.*Y32;
        FAD(:,:,:,2,1)=-(1-lamequot)./2./R+lamequot./2.*q.^2./R.^3;
        FAD(:,:,:,3,1)=-(1-lamequot)./2.*q.*Y11-lamequot./2.*q.^3.*Y32;
        FAD(:,:,:,1,2)=-(1-lamequot)./2.*(cos(angle)./R+q.*Y11.*sin(angle))-lamequot./2.*q.*F;
        FAD(:,:,:,2,2)=-(1-lamequot)./2.*yg.*X11-lamequot./2.*q.*G;
        FAD(:,:,:,3,2)=(1-lamequot)./2.*(dg.*X11+epson.*Y11.*sin(angle))+lamequot./2.*q.*H;
        FAD(:,:,:,1,3)=(1-lamequot)./2.*(sin(angle)./R-q.*Y11.*cos(angle))-lamequot./2.*q.*F1;
        FAD(:,:,:,2,3)=(1-lamequot)./2.*dg.*X11-lamequot./2.*q.*G1;
        FAD(:,:,:,3,3)=(1-lamequot)./2.*(yg.*X11+epson.*Y11.*cos(angle))+lamequot./2.*q.*H1;
    otherwise
        error('no such type');
end
%        fprintf('....');
end

% Displacement caused by surface deformation related term
function FB=FB(x,y,z,depth,angle,epson,inta,type)
global critical lama lamu lamequot;
warning off;
[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[theta,X,I4,I3,I2,I1]=DISPPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FB(:,:,:,1)=-epson.*q.*Y11-theta-(1-lamequot)./lamequot.*I1.*sin(angle);
        FB(:,:,:,2)=-q./R+(1-lamequot)./lamequot.*yg./(R+dg).*sin(angle);
        FB(:,:,:,3)=q.^2.*Y11-(1-lamequot)./lamequot.*I2.*sin(angle);
    case 'D'
        FB(:,:,:,1)=-q./R+(1-lamequot)./lamequot.*I3.*sin(angle).*cos(angle);
        FB(:,:,:,2)=-inta.*q.*X11-theta-(1-lamequot)./lamequot.*epson./(R+dg).*sin(angle).*cos(angle);
        FB(:,:,:,3)= q.^2.*X11+(1-lamequot)./lamequot.*I4.*sin(angle).*cos(angle);
    case 'T'
        FB(:,:,:,1)=q.^2.*Y11-(1-lamequot)./lamequot.*I3.*(sin(angle)).^2;
        FB(:,:,:,2)=q.^2.*X11+(1-lamequot)./lamequot.*epson./(R+dg).*(sin(angle)).^2;
        FB(:,:,:,3)=q.*(inta.*X11+epson.*Y11)-theta-(1-lamequot)./lamequot.*I4.*(sin(angle)).^2;
    otherwise
        error('no such type');
end
%        fprintf('....');
end

% derivative caused by surface deformation related term
function FBD=FBD(x,y,z,depth,angle,epson,inta,type)
global critical lamequot;
warning off;

[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[D11,K1,K2,K3,K4,J1,J2,J3,J4,J5,J6,E,F,G,H,P,Q,E1,F1,G1,H1,P1,Q1]...
    =DERIPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FBD(:,:,:,1,1)= epson.^2.*q.*Y32-(1-lamequot)./lamequot.*J1.*sin(angle);
        FBD(:,:,:,2,1)= epson.*q./R.^3-(1-lamequot)./lamequot.*J2.*sin(angle);
        FBD(:,:,:,3,1)=-epson.*q.^2.*Y32-(1-lamequot)./lamequot.*J3.*sin(angle);
        FBD(:,:,:,1,2)=-epson.*F-dg.*X11+(1-lamequot)./lamequot.*(epson.*Y11+J4).*sin(angle);
        FBD(:,:,:,2,2)=-E+(1-lamequot)./lamequot.*(1./R+J5).*sin(angle);
        FBD(:,:,:,3,2)= q.*F-(1-lamequot)./lamequot.*(q.*Y11-J6).*sin(angle);
        FBD(:,:,:,1,3)=-epson.*F1-yg.*X11+(1-lamequot)./lamequot.*K1.*sin(angle);
        FBD(:,:,:,2,3)=-E1+(1-lamequot)./lamequot.*yg.*D11.*sin(angle);
        FBD(:,:,:,3,3)= q.*F1+(1-lamequot)./lamequot.*K2.*sin(angle);
    case 'D'
        FBD(:,:,:,1,1)=epson.*q./R.^3+(1-lamequot)./lamequot.*J4.*sin(angle).*cos(angle);
        FBD(:,:,:,2,1)=inta.*q./R.^3+q.*Y11+(1-lamequot)./lamequot.*J5.*sin(angle).*cos(angle);
        FBD(:,:,:,3,1)=-q.^2./R.^3 +(1-lamequot)./lamequot.*J6.*sin(angle).*cos(angle);
        FBD(:,:,:,1,2)=-E+(1-lamequot)./lamequot.*J1.*sin(angle).*cos(angle);
        FBD(:,:,:,2,2)=-inta.*G-epson.*Y11.*sin(angle)+(1-lamequot)./lamequot.*J2.*sin(angle).*cos(angle);
        FBD(:,:,:,3,2)=q.*G+(1-lamequot)./lamequot.*J3.*sin(angle).*cos(angle);
        FBD(:,:,:,1,3)=-E1-(1-lamequot)./lamequot.*K3.*sin(angle).*cos(angle);
        FBD(:,:,:,2,3)=-inta.*G1-epson.*Y11.*cos(angle)-(1-lamequot)./lamequot.*epson.*D11.*sin(angle).*cos(angle);
        FBD(:,:,:,3,3)=q.*G1 -(1-lamequot)./lamequot.*K4.*sin(angle).*cos(angle);
    case 'T'
        FBD(:,:,:,1,1)=-epson.*q.^2.*Y32-(1-lamequot)./lamequot.*J4.*(sin(angle)).^2;
        FBD(:,:,:,2,1)=-q.^2./R.^3-(1-lamequot)./lamequot.*J5.*(sin(angle)).^2;
        FBD(:,:,:,3,1)= q.^3.*Y32-(1-lamequot)./lamequot.*J6.*(sin(angle)).^2;
        FBD(:,:,:,1,2)= q.*F-(1-lamequot)./lamequot.*J1.*(sin(angle)).^2;
        FBD(:,:,:,2,2)= q.*G-(1-lamequot)./lamequot.*J2.*(sin(angle)).^2;
        FBD(:,:,:,3,2)=-q.*H-(1-lamequot)./lamequot.*J3.*(sin(angle)).^2;
        FBD(:,:,:,1,3)= q.*F1+(1-lamequot)./lamequot.*K3.*(sin(angle)).^2;
        FBD(:,:,:,2,3)= q.*G1+(1-lamequot)./lamequot.*epson.*D11.*(sin(angle)).^2;
        FBD(:,:,:,3,3)=-q.*H1+(1-lamequot)./lamequot.*K4.*(sin(angle)).^2;
    otherwise
        error('no such type');
end
%        fprintf('....');
end

% Displacement caused by depth dependent term
function FC=FC(x,y,z,depth,angle,epson,inta,type)
global critical lamequot;
warning off;
[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[theta,X,I4,I3,I2,I1]=DISPPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FC(:,:,:,1)=(1-lamequot).*epson.*Y11.*cos(angle)-lamequot.*epson.*q.*Z32;
        FC(:,:,:,2)=(1-lamequot).*(cos(angle)./R+2.*q.*Y11.*sin(angle))-lamequot.*cg.*q./R.^3;
        FC(:,:,:,3)=(1-lamequot).*q.*Y11.*cos(angle)-lamequot.*(cg.*inta./R.^3-z.*Y11+epson.^2.*Z32);
    case 'D'
        FC(:,:,:,1)=(1-lamequot).*cos(angle)./R-q.*Y11.*sin(angle)-lamequot.*cg.*q./R.^3;
        FC(:,:,:,2)=(1-lamequot).*yg.*X11-lamequot.*cg.*inta.*q.*X32;
        FC(:,:,:,3)=-dg.*X11-epson.*Y11.*sin(angle)-lamequot.*cg.*(X11-q.^2.*X32);
    case 'T'
        FC(:,:,:,1)=-(1-lamequot).*(sin(angle)./R+q.*Y11.*cos(angle))-lamequot.*(z.*Y11-q.^2.*Z32);
        FC(:,:,:,2)=(1-lamequot).*2.*epson.*Y11.*sin(angle)+dg.*X11-lamequot.*cg.*(X11-q.^2.*X32);
        FC(:,:,:,3)=(1-lamequot).*(yg.*X11+epson.*Y11.*cos(angle))+lamequot.*q.*(cg.*inta.*X32+epson.*Z32);
    otherwise
        error('no such type');
end
%        fprintf('....');
end

% derivative caused by depth dependent term
function FCD=FCD(x,y,z,depth,angle,epson,inta,type)
global critical lamequot;
warning off;

[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);
[D11,K1,K2,K3,K4,J1,J2,J3,J4,J5,J6,E,F,G,H,P,Q,E1,F1,G1,H1,P1,Q1]...
    =DERIPARA(x,y,z,depth,angle,epson,inta);
switch upper(type)
    case 'S'
        FCD(:,:,:,1,1)=(1-lamequot).*Y0.*cos(angle)-lamequot.*q.*Z0;
        FCD(:,:,:,2,1)=-(1-lamequot).*epson.*(cos(angle)./R.^3+2.*q.*Y32.*sin(angle))...
            +lamequot.*3.*cg.*epson.*q./R.^5;
        FCD(:,:,:,3,1)=-(1-lamequot).*epson.*q.*Y32.*cos(angle)...
            +lamequot.*epson.*(3.*cg.*inta./R.^5-z.*Y32-Z32-Z0);
        FCD(:,:,:,1,2)=-(1-lamequot).*epson.*P.*cos(angle)-lamequot.*epson.*Q;
        FCD(:,:,:,2,2)=2.*(1-lamequot).*(dg./R.^3-Y0.*sin(angle)).*sin(angle)-yg./R.^3.*cos(angle)...
            -lamequot.*((cg+dg)./R.^3.*sin(angle)-inta./R.^3-3.*cg.*yg.*q./R.^5);
        FCD(:,:,:,3,2)=-(1-lamequot).*q./R.^3+(yg./R.^3-Y0.*cos(angle)).*sin(angle)...
            +lamequot.*((cg+dg)./R.^3.*cos(angle)...
            +3.*cg.*dg.*q./R.^5-(Y0.*cos(angle)+q.*Z0).*sin(angle));
        FCD(:,:,:,1,3)=(1-lamequot).*epson.*P1.*cos(angle)-lamequot.*epson.*Q1;
        FCD(:,:,:,2,3)=2.*(1-lamequot).*(yg./R.^3-Y0.*cos(angle)).*sin(angle)+dg./R.^3.*cos(angle)...
            -lamequot.*((cg+dg)./R.^3.*cos(angle)+3.*cg.*dg.*q./R.^5);
        FCD(:,:,:,3,3)=(yg./R.^3-Y0.*cos(angle)).*cos(angle)-lamequot.*((cg+dg)./R.^3.*sin(angle)...
            -3.*cg.*yg.*q./R.^5-Y0.*(sin(angle)).^2+q.*Z0.*cos(angle));
    case 'D'
        FCD(:,:,:,1,1)=-(1-lamequot).*epson./R.^3.*cos(angle)...
            +epson.*q.*Y32.*sin(angle)+lamequot.*3.*cg.*epson.*q./R.^5;
        FCD(:,:,:,2,1)=-(1-lamequot).*yg./R.^3 +lamequot.*3.*cg.*inta.*q./R.^5;
        FCD(:,:,:,3,1)=dg./R.^3-Y0.*sin(angle) +lamequot.*cg./R.^3.*(1-3.*q.^2./R.^2);
        FCD(:,:,:,1,2)=-(1-lamequot).*inta./R.^3+Y0.*(sin(angle)).^2 ...
            -lamequot.*((cg+dg)./R.^3.*sin(angle)-3.*cg.*yg.*q./R.^5);
        FCD(:,:,:,2,2)=(1-lamequot).*(X11-yg.^2.*X32)-lamequot.*cg.*((dg+2.*q.*cos(angle)).*X32...
            -yg.*inta.*q.*X53);
        FCD(:,:,:,3,2)=epson.*P.*sin(angle)+yg.*dg.*X32+lamequot.*cg.*((yg+2.*q.*sin(angle)).*X32...
            -yg.*q.^2.*X53);
        FCD(:,:,:,1,3)=-q./R.^3+Y0.*sin(angle).*cos(angle)-lamequot.*((cg+dg)./R.^3.*cos(angle)+3.*cg.*dg.*q./R.^5);
        FCD(:,:,:,2,3)=(1-lamequot).*yg.*dg.*X32-lamequot.*cg.*((yg-2.*q.*sin(angle)).*X32+dg.*inta.*q.*X53);
        FCD(:,:,:,3,3)=-epson.*P1.*sin(angle)+X11-dg.^2.*X32...
            -lamequot.*cg.*((dg-2.*q.*cos(angle)).*X32-dg.*q.^2.*X53);
    case 'T'
        FCD(:,:,:,1,1)=(1-lamequot).*epson./R.^3.*sin(angle)+epson.*q.*Y32.*cos(angle)...
            +lamequot.*epson.*(3.*cg.*inta./R.^5-2.*Z32-Z0);
        FCD(:,:,:,2,1)=(1-lamequot).*2.*Y0.*sin(angle)-dg./R.^3+lamequot.*cg./R.^3.*(1-3.*q.^2./R.^2);
        FCD(:,:,:,3,1)=-(1-lamequot).*(yg./R.^3-Y0.*cos(angle))-lamequot.*(3.*cg.*inta.*q./R.^5-q.*Z0);
        FCD(:,:,:,1,2)=(1-lamequot).*(q./R.^3+Y0.*sin(angle).*cos(angle))...
            +lamequot.*(z./R.^3.*cos(angle)+3.*cg.*dg.*q./R.^5-q.*Z0.*sin(angle));
        FCD(:,:,:,2,2)=-(1-lamequot).*2.*epson.*P.*sin(angle)-yg.*dg.*X32...
            +lamequot.*cg.*((yg+2.*q.*sin(angle)).*X32-yg.*q.^2.*X53);
        FCD(:,:,:,3,2)=-(1-lamequot).*(epson.*P.*cos(angle)-X11+yg.^2.*X32)...
            +lamequot.*cg.*((dg+2.*q.*cos(angle)).*X32-yg.*inta.*q.*X53)+lamequot.*epson.*Q;
        FCD(:,:,:,1,3)=-inta./R.^3+Y0.*(cos(angle)).^2-lamequot.*(z./R.^3.*sin(angle)...
            -3.*cg.*yg.*q./R.^5-Y0.*(sin(angle)).^2+q.*Z0.*cos(angle));
        FCD(:,:,:,2,3)=(1-lamequot).*2.*epson.*P1.*sin(angle)-X11+dg.^2.*X32...
            -lamequot.*cg.*((dg-2.*q.*cos(angle)).*X32-dg.*q.^2.*X53);
        FCD(:,:,:,3,3)=(1-lamequot).*(epson.*P1.*cos(angle)+yg.*dg.*X32)...
            +lamequot.*cg.*((yg-2.*q.*sin(angle)).*X32+dg.*inta.*q.*X53)+lamequot.*epson.*Q1;
    otherwise
        error('no such type');
end
%        fprintf('....');
end

%----------------------------------------------------------------------------------------------------
% Functions for calculating internal deformation (see Okada, 1992, Table 2-9);
% The output variables are D11,K1,K2,K3,K4,J1,J2,J3,J4,J5,J6
% E,F,G,H,P,Q,E1,F1,G1,H1,P1,Q1
%----------------------------------------------------------------------------------------------------
function [D11,K1,K2,K3,K4,J1,J2,J3,J4,J5,J6,E,F,G,H,P,Q,E1,F1,G1,H1,P1,Q1]...
    =DERIPARA(x,y,z,depth,angle,epson,inta)
global critical;
warning off;

[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
    =COMMONPARA(x,y,z,depth,angle,epson,inta);

%solving the variables of X-derivative
D11=1./(R.*(R+dg));

if abs(cos(angle))<critical
    K1=epson.*q./(R+dg).*D11;
    K3=sin(angle)./(R+dg).*(epson.^2.*D11-1);
else
    K1=epson./cos(angle).*(D11-Y11.*sin(angle));
    K3=(q.*Y11-yg.*D11)./cos(angle);
end
K2=1./R+K3.*sin(angle);
K4=epson.*Y11.*cos(angle)-K1.*sin(angle);

J2=epson.*yg./(R+dg).*D11;
if abs(cos(angle))<critical
    J3=-epson./(R+dg).^2.*(q.^2.*D11-0.5);
else
    J3=(K1-J2.*sin(angle))./cos(angle);
end
J4=-epson.*Y11-J2.*cos(angle)+J3.*sin(angle);
J5=-(dg+yg.^2./(R+dg)).*D11;
if abs(cos(angle))<critical
    J6=-yg./(R+dg).^2.*(epson.^2.*D11-0.5);
else
    J6=(K3-J5.*sin(angle))./cos(angle);
end
J1=J5.*cos(angle)-J6.*sin(angle);

%solving the variables of Y-derivative
E=sin(angle)./R-yg.*q./R.^3;
F=dg./R.^3+epson.^2.*Y32.*sin(angle);
G=2.*X11.*sin(angle)-yg.*q.*X32;
H=dg.*q.*X32+epson.*q.*Y32.*sin(angle);
P=cos(angle)./R.^3+q.*Y32.*sin(angle);
Q=3.*cg.*dg./R.^5-(z.*Y32+Z32+Z0).*sin(angle);

%solving the variables of Z-derivative
E1=cos(angle)./R+dg.*q./R.^3;
F1=yg./R.^3+epson.^2.*Y32.*cos(angle);
G1=2.*X11.*cos(angle)+dg.*q.*X32;
H1=yg.*q.*X32+epson.*q.*Y32.*cos(angle);
P1=sin(angle)./R.^3-q.*Y32.*cos(angle);
Q1=3.*cg.*yg./R.^5+q.*Y32-(z.*Y32+Z32+Z0).*cos(angle);
end

%----------------------------------------------------------------------------------------------------
% Functions for calculating internal deformation (see Okada, 1992, Table 2-9);
% The output variables are theta,X,I4,I3,I2,I1;
%----------------------------------------------------------------------------------------------------
function [theta,X,I4,I3,I2,I1]=DISPPARA(x,y,z,depth,angle,epson,inta)
global critical;
warning off;

[d,p,q,R,yg,dg,cg,X11,X32,X53,Y11,Y32,Y53,h,Z32,Z53,Y0,Z0]...
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
warning off;
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

