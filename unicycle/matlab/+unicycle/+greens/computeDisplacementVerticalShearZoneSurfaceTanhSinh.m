function [u1,u2,u3]=computeDisplacementVerticalShearZoneSurfaceTanhSinh( ...
    x1,x2,q1,q2,q3,L,T,W,theta, ...
    epsv11,epsv12,epsv13,epsv22,epsv23,epsv33,G,nu)
% function COMPUTEDISPLACEMENTVERTICALSHEARZONESURFACETANHSINH computes the
% displacement field associated with deforming vertical shear zones
% considering the following geometry.
%
%                N (x1)
%   observation /
%     point    /
%    x1,x2 -> +----------- E (x2)
%
%                        N (x1)
%                       /
%                      /\
%                     /  \ strike     .-------+
%         source     /    \  .--------        |    E (x2)
%        q1,q2,q3 ->@-------------------------+----   
%                   |                         |   + s
%                w  :                         |  / s
%                i  |                         | / e
%                d  :                         |/ n
%                t  |                 .-------+ k
%                h  :        .--------       / c
%                   +---------              / i
%                   :       l e n g t h    / h
%                   |                     + t
%                   :                      
%                   |                 
%                   D (x3)
%
% Input:
% x1, x2             north and east coordinates of the observation point
% q1, q2, q3         north, east and depth coordinates of the shear zone
% L, T, W            length, thickness and width of the shear zone
% theta (degree)     strike angle from north (from x1) of the shear zone
% epsvij             source strain component ij in the shear zone in the 
%                    reference system tied to the shear zone.
% G, nu              shear modulus and Poisson's ratio in the half space.
%
% Output:
% u1                 displacement component in the north direction,
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - Feb 22, 2016, Singapore.
%
% SEE ALSO: unicycle

% Lame parameter
lambda=G*2*nu/(1-2*nu);

% isotropic strain
epsvkk=epsv11+epsv22+epsv33;

% array size
s=size(x1);

% rotate to the shear-zone-centric system of coordinates
t1= (x1-q1)*cosd(theta)+(x2-q2)*sind(theta);
x2=-(x1-q1)*sind(theta)+(x2-q2)*cosd(theta);
x1=t1;

% Green's functions
r=@(y1,y2,y3) sqrt((x1-y1).^2+(x2-y2).^2+y3.^2);

% numerical solution
h=0.01;
n=fix(1/h*2);
u1=zeros(s);
u2=zeros(s);
u3=zeros(s);
for k=-n:n
    wk=(0.5*h*pi*cosh(k*h))./(cosh(0.5*pi*sinh(k*h))).^2;
    xk=tanh(0.5*pi*sinh(k*h));
    for j=-n:n
        wj=(0.5*h*pi*cosh(j*h))./(cosh(0.5*pi*sinh(j*h))).^2;
        xj=tanh(0.5*pi*sinh(j*h));
        u1=u1+wk*wj*IU1(xk,xj);
        u2=u2+wk*wj*IU2(xk,xj);
        u3=u3+wk*wj*IU3(xk,xj);
    end
end

% rotate displacement field to reference system of coordinates
t1=u1*cosd(theta)-u2*sind(theta);
u2=u1*sind(theta)+u2*cosd(theta);
u1=t1;


    function d = G11(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=1/(16*pi*(1-nu))*( ...
            (3-4*nu)./lr+1./lr+(x1-y1).^2./lr.^3+(3-4*nu)*(x1-y1).^2./lr.^3 ...
            +4*(1-2*nu)*(1-nu)*(lr.^2-(x1-y1).^2+lr.*y3)./(lr.*(lr+y3).^2));
    end
    function d = G12(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=(x1-y1).*(x2-y2)/(16*pi*G*(1-nu)).*( ...
            1./lr.^3+(3-4*nu)./lr.^3 ...
            -4*(1-2*nu)*(1-nu)./(lr.*(lr+y3).^2) ...
            );
    end
    function d = G13(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=(x1-y1)/(16*pi*G*(1-nu)).*( ...
            (-y3)./lr.^3+(3-4*nu)*(-y3)./lr.^3 ...
            +4*(1-2*nu)*(1-nu)./(lr.*(lr+y3)) ...
            );
    end
    function d = G21(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=(x1-y1).*(x2-y2)/(16*pi*G*(1-nu)).*( ...
            1./lr.^3+(3-4*nu)./lr.^3 ...
            -4*(1-2*nu)*(1-nu)./(lr.*(lr+y3).^2) ...
            );
    end
    function d = G22(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=1/(16*pi*(1-nu))*( ...
            (3-4*nu)./lr+1./lr+(x2-y2).^2./lr.^3 ...
            +(3-4*nu)*(x2-y2).^2./lr.^3 ...
            +4*(1-2*nu)*(1-nu)*(lr.^2-(x2-y2).^2+lr.*(y3))./(lr.*(lr+y3).^2)...
            );
    end
    function d = G23(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=1/(16*pi*G*(1-nu))*(x2-y2).*( ...
            +(4-4*nu)*(-y3)./lr.^3 ...
            +4*(1-2*nu)*(1-nu)./(lr.*(lr+y3)) ...
            );
    end
    function d = G31(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=(x1-y1)/(16*pi*G*(1-nu)).*( ...
            +(4-4*nu)*(-y3)./lr.^3 ...
            -4*(1-2*nu)*(1-nu)./(lr.*(lr+y3)) ...
            );
    end
    function d = G32(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=(x2-y2)/(16*pi*G*(1-nu)).*( ...
            (-y3)./lr.^3+(3-4*nu)*(-y3)./lr.^3 ...
            -4*(1-2*nu)*(1-nu)./(lr.*(lr+y3)) ...
            );
    end
    function d = G33(y1,y2,y3)
        lr=r(y1,y2,y3);
        d=1/(16*pi*G*(1-nu))*( ...
            (3-4*nu)./lr+(5-12*nu+8*nu^2)./lr+(-y3).^2./lr.^3 ...
            +((3-4*nu)*(+y3).^2)./lr.^3 ...
            );
    end

    function u = IU1(x,y)
        % function IU1 is the integrand for displacement component u1
        u=zeros(s);
        if epsv11 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv11)*T*W/4*(G11(L,x*T/2,(1+y)*W/2+q3)-G11(0,x*T/2,(1+y)*W/2+q3));
        end
        if epsv12 ~= 0
            u=u+2*G*epsv12*T*W/4*(G21(L,x*T/2,(1+y)*W/2+q3)-G21(0,x*T/2,(1+y)*W/2+q3)) ...
               +2*G*epsv12*L*W/4*(G11((1+x)*L/2,T/2,(1+y)*W/2+q3)-G11((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv13 ~= 0
            u=u+2*G*epsv13*T*W/4*(G31(L,x*T/2,(1+y)*W/2+q3)-G31(0,x*T/2,(1+y)*W/2+q3)) ...
               +2*G*epsv13*L*T/4*(G11((1+x)*L/2,y*T/2,W+q3)-G11((1+x)*L/2,y*T/2,q3));
        end
        if epsv22 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv22)*L*W/4*(G21((1+x)*L/2,T/2,(1+y)*W/2+q3)-G21((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv23 ~= 0
            u=u+2*G*epsv23*L*T/4*(G21((1+x)*L/2,y*T/2,W+q3)-G21((1+x)*L/2,y*T/2,q3)) ...
               +2*G*epsv23*L*W/4*(G31((1+x)*L/2,T/2,(1+y)*W/2+q3)-G31((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv33 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv33)*L*T/4*(G31((1+x)*L/2,y*T/2,W+q3)-G31((1+x)*L/2,y*T/2,q3));
        end
    end

    function u = IU2(x,y)
        % function IU2 is the integrand for displacement component u1
        u=zeros(s);
        if epsv11 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv11)*T*W/4*(G12(L,x*T/2,(1+y)*W/2+q3)-G12(0,x*T/2,(1+y)*W/2+q3));
        end
        if epsv12 ~= 0
            u=u+2*G*epsv12*T*W/4*(G22(L,x*T/2,(1+y)*W/2+q3)-G22(0,x*T/2,(1+y)*W/2+q3)) ...
               +2*G*epsv12*L*W/4*(G12((1+x)*L/2,T/2,(1+y)*W/2+q3)-G12((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv13 ~= 0
            u=u+2*G*epsv13*T*W/4*(G32(L,x*T/2,(1+y)*W/2+q3)-G32(0,x*T/2,(1+y)*W/2+q3)) ...
               +2*G*epsv13*L*T/4*(G12((1+x)*L/2,y*T/2,W+q3)-G12((1+x)*L/2,y*T/2,q3));
        end
        if epsv22 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv22)*L*W/4*(G22((1+x)*L/2,T/2,(1+y)*W/2+q3)-G22((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv23 ~= 0
            u=u+2*G*epsv23*L*T/4*(G22((1+x)*L/2,y*T/2,W+q3)-G22((1+x)*L/2,y*T/2,q3)) ...
               +2*G*epsv23*L*W/4*(G32((1+x)*L/2,T/2,(1+y)*W/2+q3)-G32((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv33 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv33)*L*T/4*(G32((1+x)*L/2,y*T/2,W+q3)-G32((1+x)*L/2,y*T/2,q3));
        end
    end

    function u = IU3(x,y)
        % function IU1 is the integrand for displacement component u1
        u=zeros(s);
        if epsv11 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv11)*T*W/4*(G13(L,x*T/2,(1+y)*W/2+q3)-G13(0,x*T/2,(1+y)*W/2+q3));
        end
        if epsv12 ~= 0
            u=u+2*G*epsv12*T*W/4*(G23(L,x*T/2,(1+y)*W/2+q3)-G23(0,x*T/2,(1+y)*W/2+q3)) ...
               +2*G*epsv12*L*W/4*(G13((1+x)*L/2,T/2,(1+y)*W/2+q3)-G13((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv13 ~= 0
            u=u+2*G*epsv13*L*T/4*(G13((1+x)*L/2,y*T/2,W+q3)-G13((1+x)*L/2,y*T/2,q3)) ...
               +2*G*epsv13*T*W/4*(G33(L,x*T/2,(1+y)*W/2+q3)-G33(0,x*T/2,(1+y)*W/2+q3));
        end
        if epsv22 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv22)*L*W/4*(G23((1+x)*L/2,T/2,(1+y)*W/2+q3)-G23((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv23 ~= 0
            u=u+2*G*epsv23*L*T/4*(G23((1+x)*L/2,y*T/2,W+q3)-G23((1+x)*L/2,y*T/2,q3)) ...
               +2*G*epsv23*L*W/4*(G33((1+x)*L/2,T/2,(1+y)*W/2+q3)-G33((1+x)*L/2,-T/2,(1+y)*W/2+q3));
        end
        if epsv33 ~= 0 || epsvkk ~= 0
            u=u+(lambda*epsvkk+2*G*epsv33)*L*T/4*(G33((1+x)*L/2,y*T/2,W+q3)-G33((1+x)*L/2,y*T/2,q3));
        end
    end

end
