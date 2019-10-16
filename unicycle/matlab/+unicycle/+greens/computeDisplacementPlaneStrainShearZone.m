function [u2,u3]=computeDisplacementPlaneStrainShearZone( ...
    x2,x3,q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu)
% function COMPUTEDISPLACEMENTPLANESTRAINSHEARZONE computes the
% displacement field associated with deforming dipping shear zones
% considering the following geometry using the analytical solution.
%
%              surface
%      -------------+-------------- E (x2)
%                   |
%                   | dip /
%                   |----/  . w
%                   |   /     . i 
%                   |  /        . d           
%                   | /           . t     
%                   |/              . h   
%           q2,q3 ->@                 .
%                  /|                   . 
%                 / :                  /  
%                /  |                 /  s
%               /   :                /  s
%              /    |               /  e
%             /     :              /  n
%               .   |             /  k
%                 . :            /  c
%                   .           /  i
%                   : .        /  h
%                   |   .     /  t
%                   :     .  /  
%                   |       .    
%                   q3 (x3)
%
%
% Input:
% x2, x3             east coordinates and depth of the observation point,
% q2, q3             east and depth coordinates of the shear zone,
% T, W               thickness and width of the shear zone,
% phi (degree)       dip angle from horizontal of the shear zone,
% epsvij             source strain component 22, 23 and 33 in the shear zone
%                    in the system of reference tied to the shear zone,
% G, nu              shear modulus and Poisson's ratio in the half space.
%
% Output:
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) and James Moore (earth@jamesdpmoore.com)
%         March 19, 2016, Singapore.
%
% SEE ALSO: unicycle
%-----------------------------------------------------------------------

assert(0<=q3,'unicycle.greens.computeDisplacementPlaneStrainShearZone: depth of source must be positive.');
assert(0<=min(x3),'unicycle.greens.computeDisplacementPlaneStrainShearZone: depth of observation point must be positive.');

% Lame parameter
lambda=G*2*nu/(1-2*nu);

% convert strain to reference coordinate system R=[[sind(phi), cosd(phi)];[-cosd(phi), sind(phi)]];
epsv22=+sind(phi)*(+epsv22p*sind(phi)+epsv23p*cosd(phi))+cosd(phi)*(+epsv23p*sind(phi)+epsv33p*cosd(phi));
epsv23=+sind(phi)*(-epsv22p*cosd(phi)+epsv23p*sind(phi))+cosd(phi)*(-epsv23p*cosd(phi)+epsv33p*sind(phi));
epsv33=-cosd(phi)*(-epsv22p*cosd(phi)+epsv23p*sind(phi))+sind(phi)*(-epsv23p*cosd(phi)+epsv33p*sind(phi));

% isotropic strain
epsvkk=epsv22+epsv33;

% translation
x2=x2-q2;

% Green's functions
y2 =@(y2p,y3p)  y2p*sind(phi)+y3p*cosd(phi);
y3 =@(y2p,y3p) -y2p*cosd(phi)+y3p*sind(phi)+q3;
r1p=@(y2p,y3p)  sqrt((x3-y3(y2p,y3p)).^2+(x2-y2(y2p,y3p)).^2);
r2p=@(y2p,y3p)  sqrt((x3+y3(y2p,y3p)).^2+(x2-y2(y2p,y3p)).^2);
p2 =@(y2p,y3p)  (x3-q3)*cosd(phi)-sind(phi)*x2+y2p;
p3 =@(y2p,y3p) -(x3-q3)*sind(phi)-cosd(phi)*x2+y3p;
p2p=@(y2p,y3p) -(q3+x3)*cosd(phi)-sind(phi)*x2+y2p;
p3p=@(y2p,y3p)  (q3+x3)*sind(phi)-cosd(phi)*x2+y3p;

r12p=@(y2p,y3p) q3.^2+x2.^2+(-2).*q3.*x3+x3.^2+y2p.^2+y3p.^2+(-2).*q3.*y2p.*cosd(phi)+ ...
  2.*x3.*y2p.*cosd(phi)+(-2).*x2.*y3p.*cosd(phi)+(-2).*x2.*y2p.*sind(phi)+2.* ...
  q3.*y3p.*sind(phi)+(-2).*x3.*y3p.*sind(phi);

r22p=@(y2p,y3p) q3.^2+x2.^2+2.*q3.*x3+x3.^2+y2p.^2+y3p.^2+(-2).*q3.*y2p.*cosd(phi)+(-2) ...
  .*x3.*y2p.*cosd(phi)+(-2).*x2.*y3p.*cosd(phi)+(-2).*x2.*y2p.*sind(phi)+2.* ...
  q3.*y3p.*sind(phi)+2.*x3.*y3p.*sind(phi);

I223=@(y2p,y3p) ...
    (-1).*p2(y2p,y3p).*atan2(p3(y2p,y3p)./p2(y2p,y3p),1).*(3-4*nu+cosd(2*phi)) ...
    +(1/4).*log(r12p(y2p,y3p)).*( ...
        -6*y3p+8*y3p*nu ...
        +x2.*(5-8*nu).*cosd(phi) ...
        +x2*cosd(3*phi) ...
        -2*(q3-x3)*(4-4*nu+cosd(2*phi))*sind(phi) ...
        +2*y2p*sind(2*phi) ...
    ) ...
    +(1/2).*atan2(-p3p(y2p,y3p)./p2p(y2p,y3p),1).*( ...
        2*y2p*(5+4*nu*(-3+2*nu)) ...
        -(13*q3+11*x3-28*(q3+x3)*nu+16*(q3+x3)*nu^2)*cosd(phi) ...
        -2*y2p*(-3+4*nu)*cosd(2*phi) ...
        -(3*q3+5*x3-4*(q3+x3)*nu)*cosd(3*phi) ...
        -x2*(7+4*nu*(-5+4*nu))*sind(phi) ...
        -x2*(3-4*nu)*sind(3*phi) ...
    ) ...
    +(1/4)*log(r22p(y2p,y3p)).*( ...
        -2*y3p*(5+4*nu*(-3+2*nu)) ...
        +x2*(7+4*nu*(-5+4*nu))*cosd(phi) ...
        +x2*(3-4*nu)*cosd(3*phi) ...
        -(13*q3+11*x3-28*(q3+x3)*nu+16*(q3+x3)*nu^2)*sind(phi) ...
        +2*y2p*(3-4*nu)*sind(2*phi) ...
        +(-3*q3-5*x3+4*(q3+x3)*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        -2*q3*x2.*cosd(phi) ...
        +(3*x2*y2p+q3*y3p-x3*y3p)*cosd(2*phi) ...
        -2*y2p*y3p*cosd(3*phi) ...
        +(-x2*y2p+(q3+x3)*y3p)*cosd(4*phi) ...
        +(-q3^2+x2.^2+x3.^2)*sind(phi) ...
        +(3.*q3.*y2p+x3.*y2p-x2*y3p)*sind(2*phi) ...
        -(x2.^2+(q3+x3).^2+2*y2p^2)*sind(3*phi) ...
        +((q3+x3)*y2p+x2*y3p)*sind(4*phi) ...
    );

I222=@(y2p,y3p) ...
    p3(y2p,y3p).*atan2(p2(y2p,y3p)./p3(y2p,y3p),1)*(-3+4*nu+cosd(2*phi)) ...
    +(1/4)*log(r12p(y2p,y3p)).*( ...
        -6*y2p ...
        +8*y2p*nu ...
        -2*(q3-x3)*cosd(phi)*(-4+4*nu+cosd(2*phi)) ...
        -2*x2*(-2+4*nu+cosd(2*phi))*sind(phi) ...
        +2*y3p*sind(2*phi) ...
	) ...
	+(1/4).*log(r22p(y2p,y3p)).*( ...
        -2*y2p*(5+4*nu*((-3)+2*nu)) ...
        +(13*q3+11*x3-28*(q3+x3)*nu+16*(q3+x3)*nu^2)*cosd(phi) ...
        +(-3*q3-5*x3+4*(q3+x3)*nu)*cosd(3*phi) ...
        +x2*(7+4*nu*(-5+4*nu))*sind(phi) ...
        +2*y3p*(3-4*nu)*sind(2*phi) ...
        +x2*(-3+4*nu)*sind(3*phi) ...
    ) ...
    +(1/2).*atan2(p2p(y2p,y3p)./p3p(y2p,y3p),1).*( ...
        -2*y3p*(5+4*nu*(-3+2*nu)) ...
        -x2*(-7+4*(5-4*nu)*nu)*cosd(phi) ...
        -2*y3p*(-3+4*nu)*cosd(2*phi) ...
        -x2*(3-4*nu)*cosd(3*phi) ...
        -(13*q3+11*x3)*sind(phi) ...
        -4*(q3+x3)*nu*(-7+4*nu)*sind(phi) ...
        -(-3*q3-5*x3+4*(q3+x3)*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        -(-q3^2+x2.^2+x3.^2)*cosd(phi) ...
        +(-q3*y2p+x3*y2p+3*x2*y3p)*cosd(2*phi)+ ...
        -(x2.^2+(q3+x3).^2+2*y3p^2)*cosd(3*phi) ...
        +((q3+x3)*y2p+x2*y3p)*cosd(4*phi) ...
        -2*q3*x2*sind(phi) ...
        +(x2*y2p+(3*q3+x3)*y3p)*sind(2*phi) ...
        -2*y2p*y3p*sind(3*phi) ...
        +(x2*y2p-(q3+x3)*y3p)*sind(4*phi) ...
	);

I233=@(y2p,y3p) ...
    4*y3p*(-1+nu)*(-1+2*nu)*atan2((x2-y2(y2p,y3p))./(x3+y3(y2p,y3p)),1) ...
    +(-1/2).*p2(y2p,y3p).*cosd(2*phi).*log(r12p(y2p,y3p)) ...
    +4*((-1)+nu).*((-1)+2.*nu).*atan2(p2p(y2p,y3p).*p3p(y2p,y3p).^(-1),1).*(x2*cosd(phi)-(q3+x3)*sind(phi)) ...
    -p2(y2p,y3p).*atan2(p2(y2p,y3p).^(-1).*p3(y2p,y3p),1)*sind(2*phi) ...
    +atan2((-1).*p2p(y2p,y3p).^(-1).*p3p(y2p,y3p),1)*sind(phi).*( ...
        2*y2p.*(3+(-4).*nu)*cosd(phi) ...
        +(-3*q3-x3+4*(q3+x3)*nu)*cosd(2*phi) ...
        +(-3+4*nu)*(q3-x3+x2*sind(2*phi)) ...
    ) ...
    +0.25*log(r22p(y2p,y3p)).*( ...
        -8*y2p*(-1+nu)*(-1+2*nu) ...
        +(11*q3+x3-4*(7*q3+3*x3)*nu+16*(q3+x3)*nu^2)*cosd(phi) ...
        +2*y2p*(-3+4*nu)*cosd(2*phi) ...
        +(3.*q3+x3-4*(q3+x3)*nu)*cosd(3*phi) ...
        +x2*(5+4*nu*(-5+4*nu))*sind(phi) ...
        +x2*(3-4*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        (-q3.^2+x2.^2+x3.^2).*cosd(phi) ...
        +(3*q3*y2p+x3*y2p-x2*y3p)*cosd(2*phi) ...
        +(-1).*(x2.^2+(q3+x3).^2+2*y2p^2).*cosd(3*phi) ...
        +((q3+x3).*y2p+x2.*y3p)*cosd(4*phi) ...
        +2*q3*x2*sind(phi) ...
        +(-3*x2.*y2p-q3*y3p+x3*y3p)*sind(2*phi) ...
        +2*y2p*y3p*sind(3*phi) ...
        +(x2*y2p-(q3+x3)*y3p)*sind(4*phi) ...
	);

I232=@(y2p,y3p) ...
    4*y2p*(-1+nu)*(-1+2*nu)*atan2( ...
        (x2-y3p.*cosd(phi)-y2p*sind(phi)) ...
        ./(q3+x3-y2p*cosd(phi)+y3p*sind(phi)),1) ...
    -(1/2)*p3(y2p,y3p)*cosd(2*phi).*log(r12p(y2p,y3p)) ...
    -4*(-1+nu).*(-1+2.*nu).*atan2(p2p(y2p,y3p).^(-1).*p3p(y2p,y3p),1).*((q3+x3)*cosd(phi)+x2*sind(phi)) ...
    +p3(y2p,y3p).*atan2(p2(y2p,y3p).*p3(y2p,y3p).^(-1),1).*sind(2*phi) ...
    -atan2(p2p(y2p,y3p).*p3p(y2p,y3p).^(-1),1)*cosd(phi).*( ...
        (3*q3+x3-4*(q3+x3)*nu)*cosd(2*phi) ...
        +(-3+4*nu)*(q3-x3+2*y3p*sind(phi)-x2.*sind(2*phi)) ...
    ) ...
    +0.25*log(r22p(y2p,y3p)).*( ...
        8*y3p.*(-1+nu)*(-1+2*nu) ...
        +x2*(-5+4*(5-4*nu)*nu)*cosd(phi) ...
        +2*y3p*(-3+4*nu)*cosd(2*phi) ...
        +x2*(3-4*nu)*cosd(3*phi) ...
        +(11*q3+x3-4*(7*q3+3*x3)*nu+16*(q3+x3)*nu^2)*sind(phi) ...
        +(-3*q3-x3+4*(q3+x3)*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        -2*q3*x2*cosd(phi) ...
        +(x2.*y2p+3*q3*y3p+x3*y3p)*cosd(2*phi) ...
        -2*y2p*y3p*cosd(3*phi) ...
        +(x2*y2p-(q3+x3)*y3p)*cosd(4*phi) ...
        +(-q3^2+x2.^2+x3.^2)*sind(phi) ...
        +(q3.*y2p-x3*y2p-3*x2*y3p)*sind(2*phi) ...
        +(x2.^2+(q3+x3).^2+2*y3p^2)*sind(3*phi) ...
        -((q3+x3)*y2p+x2*y3p)*sind(4*phi) ...
    );

I323=@(y2p,y3p) ...
    -4*y3p*(-1+nu)*(-1+2*nu).*atan2( ...
        (x2-y3p*cosd(phi)-y2p*sind(phi)) ...
        ./(q3+x3+(-1).*y2p.*cosd(phi)+y3p.*sind(phi)),1) ...
    -(1/2)*p2(y2p,y3p).*cosd(2*phi).*log(r12p(y2p,y3p)) ...
    -4*((-1)+nu).*((-1)+2.*nu).*atan2(p2p(y2p,y3p).*p3p(y2p,y3p).^(-1),1).*(x2*cosd(phi)-(q3+x3)*sind(phi)) ...
    -p2(y2p,y3p).*atan2(p2(y2p,y3p).^(-1).*p3(y2p,y3p),1)*sind(2*phi) ...
    +atan2((-1).*p2p(y2p,y3p).^(-1).*p3p(y2p,y3p),1)*sind(phi).*( ...
        2.*y2p.*(3+(-4).*nu).*cosd(phi) ...
        +((-3).*q3+(-5).*x3 ...
        +4.*(q3+x3).*nu).*cosd(2*phi) ...
        +((-3)+4.*nu).*(q3+(-1).*x3+x2.*sind(2*phi)) ...
    ) ...
    +0.25*log(r22p(y2p,y3p)).*( ...
        8*y2p*(-1+nu)*(-1+2*nu) ...
        -(x3*(19+4*nu*(-9+4*nu)) ...
        +q3*(5+4*nu*(-5+4*nu)))*cosd(phi) ...
        +2*y2p*(-3+4*nu)*cosd(2*phi) ...
        +(3*q3+5*x3-4*(q3+x3)*nu)*cosd(3*phi) ...
        +x2*(-11+4*(7-4*nu)*nu)*sind(phi) ...
        +x2*(3-4*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        (q3^2-x2.^2-x3.^2)*cosd(phi) ...
        +(-(3*q3+x3)*y2p+x2*y3p)*cosd(2*phi) ...
        +(x2.^2+(q3+x3).^2+2*y2p^2)*cosd(3*phi) ...
        -((q3+x3)*y2p+x2*y3p)*cosd(4*phi) ...
        -2*q3*x2*sind(phi) ...
        +(3*x2*y2p+q3*y3p-x3*y3p)*sind(2*phi) ...
        -2*y2p*y3p*sind(3*phi) ...
        +(-x2*y2p+(q3+x3)*y3p)*sind(4*phi) ...
    );

I322=@(y2p,y3p) ...
    -4*y2p.*((-1)+nu).*((-1)+2.*nu).*atan2( ...
        (x2+(-1).*y3p.*cosd(phi)+(-1).*y2p.*sind(phi)) ...
        .*(q3+x3+(-1).*y2p.*cosd(phi)+y3p.*sind(phi)).^(-1),1) ...
    -(0.5)*p3(y2p,y3p).*cosd(2*phi).*log(r12p(y2p,y3p)) ...
    +4*((-1)+nu).*((-1)+2.*nu).*atan2( ...
        p2p(y2p,y3p).^(-1).*p3p(y2p,y3p),1).*((q3+x3).*cosd(phi)+x2.*sind(phi)) ...
    +atan2(p2p(y2p,y3p).*p3p(y2p,y3p).^(-1),1)*cosd(phi).*( ...
        (-3*q3-5*x3+4*(q3+x3)*nu)*cosd(2*phi) ...
        -(-3+4*nu)*(q3-x3+2*(y3p-x2*cosd(phi))*sind(phi)) ...
    ) ...
    +p3(y2p,y3p).*atan2(p2(y2p,y3p)./p3(y2p,y3p),1)*sind(2*phi) ...
    +0.25*log(r22p(y2p,y3p)).*( ...
        -8*y3p*(-1+nu)*(-1+2*nu) ...
        +x2*(11+4*nu*(-7+4*nu))*cosd(phi) ...
        +2*y3p*(-3+4*nu)*cosd(2*phi) ...
        +x2*(3-4*nu)*cosd(3*phi) ...
        -(5*q3+19*x3-4*(5*q3+9*x3)*nu+16*(q3+x3)*nu^2)*sind(phi) ...
        +(-3*q3-5*x3+4*(q3+x3)*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        2.*q3.*x2.*cosd(phi) ...
        -(x2*y2p+3*q3*y3p+x3*y3p)*cosd(2*phi) ...
        +2*y2p*y3p*cosd(3*phi) ...
        +(-x2*y2p+(q3+x3)*y3p)*cosd(4*phi) ...
        -(-q3^2+x2.^2+x3.^2)*sind(phi) ...
        +(-q3*y2p+x3*y2p+3*x2*y3p)*sind(2.*phi) ...
        -(x2.^2+(q3+x3).^2+2*y3p^2)*sind(3*phi) ...
        +((q3+x3)*y2p+x2*y3p)*sind(4*phi) ...
    );

I333=@(y2p,y3p) ...
    p2(y2p,y3p).*atan2(p3(y2p,y3p)./p2(y2p,y3p),1).*((-3)+4.*nu+cosd(2*phi)) ...
    +(1/4)*log(r12p(y2p,y3p)).*( ...
        -6*y3p+8*y3p*nu ...
        +x2*(7-8*nu)*cosd(phi) ...
        -x2*cosd(3*phi) ...
        +2*(q3-x3)*(-2+4*nu+cosd(2*phi))*sind(phi) ...
        -2*y2p*sind(2*phi) ...
    ) ...
    +0.5*atan2(-p3p(y2p,y3p)./p2p(y2p,y3p),1).*( ...
        2*y2p*(5+4*nu*(-3+2*nu)) ...
        -(7*q3+5*x3-20*(q3+x3)*nu+16*(q3+x3)*nu^2)*cosd(phi) ...
        +2.*y2p*(-3+4*nu)*cosd(2*phi) ...
        +(3*q3+x3-4*(q3+x3)*nu)*cosd(3*phi) ...
        +x2*(-13+4*(7-4*nu)*nu)*sind(phi) ...
        +x2*(3-4*nu)*sind(3*phi) ...
    ) ...
    +0.25*log(r22p(y2p,y3p)).*( ...
        -2*y3p*(5+4*nu*(-3+2*nu)) ...
        +x2*(13+4*nu*(-7+4*nu))*cosd(phi) ...
        +x2*(-3+4*nu)*cosd(3*phi) ...
        -(7*q3+5*x3-20*(q3+x3)*nu+16*(q3+x3)*nu^2)*sind(phi) ...
        +2*y2p*(-3+4*nu)*sind(2*phi) ...
        +(3*q3+x3-4*(q3+x3)*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        -2*q3*x2*cosd(phi) ...
        +(3*x2*y2p+q3*y3p-x3*y3p)*cosd(2*phi) ...
        -2*y2p.*y3p.*cosd(3*phi) ...
        +(-x2*y2p+(q3+x3)*y3p)*cosd(4*phi)+(-q3.^2+x2.^2+x3.^2)*sind(phi) ...
        +(3*q3*y2p+x3*y2p-x2*y3p)*sind(2*phi) ...
        -(x2.^2+(q3+x3).^2+2*y2p^2)*sind(3*phi) ...
        +((q3+x3)*y2p+x2*y3p)*sind(4*phi) ...
    );

I332=@(y2p,y3p) ...
    -p3(y2p,y3p).*atan2(p2(y2p,y3p)./p3(y2p,y3p),1)*(3-4*nu+cosd(2*phi)) ...
    +(1/4).*log(r12p(y2p,y3p)).*( ...
        -6*y2p ...
        +8*y2p*nu ...
        +2*(q3-x3)*cosd(phi)*(2-4*nu+cosd(2*phi)) ...
        +2*x2*(4-4*nu+cosd(2*phi))*sind(phi) ...
        -2*y3p*sind(2*phi) ...
    ) ...
    +0.25*log(r22p(y2p,y3p)).*( ...
        (-2).*y2p.*(5+4.*nu.*((-3)+2.*nu)) ...
        +(7*q3+5.*x3+(-20).*(q3+x3).*nu+16.*(q3+x3).*nu.^2).*cosd(phi) ...
        +(3.*q3+x3+(-4).*(q3+x3).*nu).*cosd(3*phi) ...
        +x2.*(13+4.*nu.*((-7)+4.*nu)).*sind(phi) ...
        +2.*y3p.*((-3)+4.*nu).*sind(2*phi) ...
        +x2.*(3+(-4).*nu).*sind(3*phi) ...
    ) ...
    +0.5*atan2(p2p(y2p,y3p)./p3p(y2p,y3p),1).*( ...
        -2*y3p*(5+4*nu*(-3+2*nu)) ...
        -x2*(-13+4*(7-4*nu)*nu)*cosd(phi) ...
        -2*y3p*(3-4*nu)*cosd(2*phi) ...
        -x2.*((-3)+4.*nu)*cosd(3*phi) ...
        -(7.*q3+5.*x3)*sind(phi) ...
        -4*(q3+x3).*nu.*((-5)+4.*nu).*sind(phi) ...
        -(3.*q3+x3+(-4).*(q3+x3)*nu)*sind(3*phi) ...
    ) ...
    +r2p(y2p,y3p).^(-2).*x3.*( ...
        (-1).*((-1).*q3.^2+x2.^2+x3.^2)*cosd(phi) ...
        +((-1).*q3.*y2p+x3.*y2p+3.*x2.*y3p)*cosd(2*phi) ...
        -(x2.^2+(q3+x3).^2+2.*y3p.^2)*cosd(3*phi) ...
        +((q3+x3).*y2p+x2.*y3p)*cosd(4*phi) ...
        +(-2).*q3.*x2.*sind(phi) ...
        +(x2.*y2p+(3.*q3+x3).*y3p)*sind(2*phi) ...
        -2*y2p*y3p*sind(3*phi) ...
        +(x2.*y2p+(-1).*(q3+x3).*y3p)*sind(4*phi) ...
	);

IU2=@(y2p,y3p) ...
    1/(8*pi*G*(1-nu))*( ...
     sind(phi)*((lambda*epsvkk+2*G*epsv22)*I223(y2p,y3p) ...
                              +2*G*epsv23*(I222(y2p,y3p)+I323(y2p,y3p)) ...
               +(lambda*epsvkk+2*G*epsv33)*I322(y2p,y3p) ) ...
    +cosd(phi)*((lambda*epsvkk+2*G*epsv22)*I222(y2p,y3p) ...
                              +2*G*epsv23*(I322(y2p,y3p)-I223(y2p,y3p)) ...
               -(lambda*epsvkk+2*G*epsv33)*I323(y2p,y3p) ) ...
        );

IU3=@(y2p,y3p) ...
    1/(8*pi*G*(1-nu))*( ...
     sind(phi)*((lambda*epsvkk+2*G*epsv22)*I233(y2p,y3p) ...
                              +2*G*epsv23*(I232(y2p,y3p)+I333(y2p,y3p)) ...
               +(lambda*epsvkk+2*G*epsv33)*I332(y2p,y3p) ) ...
    +cosd(phi)*((lambda*epsvkk+2*G*epsv22)*I232(y2p,y3p) ...
                              +2*G*epsv23*(I332(y2p,y3p)-I233(y2p,y3p)) ...
               -(lambda*epsvkk+2*G*epsv33)*I333(y2p,y3p) ) ...
        );

u2=IU2(T/2,W)-IU2(-T/2,W)+IU2(-T/2,0)-IU2(T/2,0);
u3=IU3(T/2,W)-IU3(-T/2,W)+IU3(-T/2,0)-IU3(T/2,0);

end







