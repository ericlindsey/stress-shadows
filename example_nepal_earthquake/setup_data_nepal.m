
% set up and check faults and data for Nepal inversion
clear all, close all

%%
% load, test and plot geometry
mu=30e9; nu=1/4;
fltname='faults/qiu+15_1';
geom=unicycle.geometry.triangleReceiver(fltname,unicycle.greens.nikkhoo15(mu,nu));

figure(1),clf
geom.plotPatch(geom.xc(:,3)), hold on
%geom.plotPatch() % an empty call to this function adds a wireframe
view(40,30)
daspect([1 1 1])

% fault azimuth from USGS GCMT
az = 190;
% rake file columns: [rake_angle,rake_ss,rake_ds,rakeE,rakeN]
rake = create_rake_file( geom, az, 1, 'faults/nepal_approx_rake.dat');
geom.plotSlipVectors(rake(:,2),rake(:,3),4000)


%%

% GPS data reference: Huang et al., 2017. http://dx.doi.org/10.1016/j.geog.2017.03.003
offsets = readtable('data/aria_offsets_reordered.txt');
coordinates = readtable('data/gps_projected_coordinates.txt');

%%

figure(1)

plot(coordinates.x2,coordinates.x1,'k.')

scaled_quiver(coordinates.x2,coordinates.x1,offsets.E_cm_,offsets.N_cm_,1e3,{'r'})

% format GPS data for jointinv
% For type 'dat': data is a simple ASCII table of size N x 6, containing locations and displacements:
% [lon,lat,hgt,uE,uN,uZ,sigE,sigN,sigZ,[wgt]]

% also, convert to meters
gps_out = [coordinates.x2,coordinates.x1,coordinates.x3, offsets.E_cm_/100,offsets.N_cm_/100,offsets.U_cm_/100, ...
    offsets.E_sig_/100,offsets.N_sig_/100,offsets.U_sig_/100];
save('data/nepal_aria_gps_formatted.dat','gps_out','-ASCII');

