function [xy] = polyconic(lat, diffLon, lat0)
	
%%	Polyconic Projection
%
%	Input:
%   lat 		= latitude (decimal seconds)
%	diffLon 	= Differential Longitude (decimal seconds) relative to Central Meridian
%   lat0        = latitude of Origin (decimal seconds)
%
%	Output: 
%   x =	Distance from Central Meridian
%	y = Distance from Origin to latitude
%
% From PCAIM

p1 = lat0; 
p2 = lat;
il = diffLon;

arcone = 4.8481368e-6;
esq = 6.7686580e-3;
la = 6378206.4;
a0 = 6367399.7;
a2 = 32433.888;
a4 = 34.4187;
a6 = .0454;
a8 = 6.0e-5;

ip = p2 - p1;
sinp2 = sin(p2 * arcone);
cosp2 = cos(p2 * arcone);
theta = il * sinp2;
a = sqrt(1.0 - (esq * (2. * sinp2))) / (la * arcone);
cot = cosp2 / sinp2;
x = (cot * sin(theta * arcone)) / (a * arcone);
ipr = ip * arcone;
pr = ((p2 + p1) / 2.)*arcone;
y = ((((a0*ipr) - ((a2*cos(2.*pr))*sin(ipr))) + ((a4*cos(4.*pr))*sin(2.*ipr))) - ((a6*cos(6.*pr))*sin(3.*ipr))) + ((a8*cos(8.*pr))*sin(4.*ipr));
xy = [x, y];

