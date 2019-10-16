function [omega]=euler_vector(latp,lonp,degmyr)
%Convert geodetic euler pole (lat,lon,deg/myr) to a scaled ECEF rotation
%vector scaled to give units of mm/yr when multiplied by a rotation matrix based on 
%spherical earth coordinates.
px=cosd(lonp)*cosd(latp);
py=sind(lonp)*cosd(latp);
pz=sind(latp);
omega=6371*degmyr*pi/180*[px;py;pz];
end