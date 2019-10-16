function u1=computeDisplacementAntiplaneShearZone( ...
    x2,x3,q2,q3,T,W,phi,epsv12p,epsv13p)
% function COMPUTEDISPLACEMENTANTIPLANESHEARZONE computes the
% displacement field associated with deforming shear zones
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
%                   D (x3)
%
% Input:
% x2, x3             east coordinates and depth of the observation point,
% q2, q3             east and depth coordinates of the shear zone,
% T, W               thickness and width of the shear zone,
% epsv12, epsv13     source strain component 22, 23 and 33 in the shear zone
%                    in the system of reference tied to the shear zone
%
% Output:
% u1                 displacement component in the north (along-strike)
%                    direction.
%
% AUTHOR: Sylvain Barbot (sbarbot@ntu.edu.sg) - March 6, 2016, Singapore.
%-----------------------------------------------------------------------

epsv12= sind(phi)*epsv12p+cosd(phi)*epsv13p;
epsv13=-cosd(phi)*epsv12p+sind(phi)*epsv13p;

y2=@(y2p,y3p) +y2p*sind(phi)+y3p*cosd(phi)+q2;
y3=@(y2p,y3p) -y2p*cosd(phi)+y3p*sind(phi)+q3;

% Green's functions
r2=@(y2p,y3p) x2-y2(y2p,y3p);
r3=@(y2p,y3p) x3-y3(y2p,y3p);

r2p=@(y2p,y3p) r2(y2p,y3p)*sind(phi)-r3(y2p,y3p)*cosd(phi);
r3p=@(y2p,y3p) r2(y2p,y3p)*cosd(phi)+r3(y2p,y3p)*sind(phi);

s3=@(y2p,y3p) x3+y3(y2p,y3p);

s2p=@(y2p,y3p) r2(y2p,y3p)*sind(phi)+s3(y2p,y3p)*cosd(phi);
s3p=@(y2p,y3p) r2(y2p,y3p)*cosd(phi)-s3(y2p,y3p)*sind(phi);

J12=@(y2p,y3p) -2*r2p(y2p,y3p).*atan(r3p(y2p,y3p)./r2p(y2p,y3p))-r3p(y2p,y3p).*log(r2(y2p,y3p).^2+r3(y2p,y3p).^2);
K12=@(y2p,y3p) -2*s2p(y2p,y3p).*atan(s3p(y2p,y3p)./s2p(y2p,y3p))-s3p(y2p,y3p).*log(r2(y2p,y3p).^2+s3(y2p,y3p).^2);

J13=@(y2p,y3p) -2*r3p(y2p,y3p).*atan(r2p(y2p,y3p)./r3p(y2p,y3p))-r2p(y2p,y3p).*log(r2(y2p,y3p).^2+r3(y2p,y3p).^2);
K13=@(y2p,y3p) -2*s3p(y2p,y3p).*atan(s2p(y2p,y3p)./s3p(y2p,y3p))-s2p(y2p,y3p).*log(r2(y2p,y3p).^2+s3(y2p,y3p).^2);

IU1=@(y2p,y3p) -1/2/pi*( ...
     sind(phi)*(epsv12*(J12(y2p,y3p)+K12(y2p,y3p))+epsv13*(J13(y2p,y3p)+K13(y2p,y3p))) ...
    +cosd(phi)*(epsv12*(J13(y2p,y3p)+K13(y2p,y3p))-epsv13*(J12(y2p,y3p)+K12(y2p,y3p))) ...
    );

u1=IU1(T/2,W)-IU1(-T/2,W)-IU1(T/2,0)+IU1(-T/2,0);

end

