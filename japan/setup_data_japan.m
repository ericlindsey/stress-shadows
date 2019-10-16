% Create necessary jointinv2 input files for Japan backslip model

clear, close all

run('../functions/set_jointinv_path')
unicyclepath='/Users/elindsey/Dropbox/code/geodesy/unicycle';
set_unicycle_path(unicyclepath)
import unicycle.*

% set projection center
lon0=143;
lat0=39;


%% convert Japan trench megathrust geometry to local XYZ

!cp Slab2/japan_trench_geo2.tri Slab2/japan_trench_xyz.tri

flt_geo=load('Slab2/japan_trench_geo2.ned');
[fltXkm,fltYkm]=latlon_to_xy_polyconic(flt_geo(:,2),flt_geo(:,3),lat0,lon0);

!mkdir faults

%% extract polygons and coastline

figure(4),clf
plot3(flt_geo(:,3),flt_geo(:,2),-flt_geo(:,4),'k.'),hold on
%
daspect([1 1 100])
Izero = find(flt_geo(:,4)==0 & flt_geo(:,3)> 142 & flt_geo(:,2)< 42.2);
plot(flt_geo(Izero,3),flt_geo(Izero,2),'rx')
% sort the points by latitude
[tracelat,Isort]=sort(flt_geo(Izero,2),'ascend');
tracelon=flt_geo(Izero,3);
tracelon=tracelon(Isort);
% save fault trace
flt_trace = [tracelon,tracelat];
save('faults/japan_trace.dat','flt_trace','-ASCII');
%
I100 = find(flt_geo(:,4)==100);
plot3(flt_geo(I100,3),flt_geo(I100,2),-flt_geo(I100,4),'gx')
% sort the points reversed by latitude
[deeptracelat,Isort]=sort(flt_geo(I100,2),'descend');
deeptracelon=flt_geo(I100,3);
deeptracelon=deeptracelon(Isort);
flt_polygon = [tracelon, tracelat; deeptracelon, deeptracelat; tracelon(1), tracelat(1)];
save('faults/japan_polygon.dat','flt_polygon','-ASCII');

% not exactly the polygon... but good enough since we will crop the edges
% in the plot anyway

plot(flt_polygon(:,1),flt_polygon(:,2),'-b')
view(2)


%%


fileID = fopen('Slab2/japan_trench_xyz.ned','w');
for i=1:length(fltXkm)
    fprintf(fileID,'%d %.9f %.9f %.9f\n',flt_geo(i,1), fltYkm(i)*1e3, fltXkm(i)*1e3, flt_geo(i,4)*1e3);
end
fclose(fileID);

%%


% load and plot geometry
mu=30e9; nu=1/4;
fltname='Slab2/japan_trench_xyz';
megathrust_full=unicycle.geometry.triangleReceiver(fltname,unicycle.greens.nikkhoo15(mu,nu));

fltname_geo='Slab2/japan_trench_geo2';
megathrust_geo=unicycle.geometry.triangleReceiver(fltname_geo,unicycle.greens.nikkhoo15(mu,nu));

% drop some values

crop_sz=false;

% if crop_sz
% 
%     x=linspace(-2e5,3e5,10);
%     y=3e5-1.0*x;
% 
% 
%     figure(2), clf
%     megathrust_full.plotPatch()
%     view(2)
%     daspect([1 1 1]), hold on
% 
%     plot(x,y);
% 
% 
%     Idrop=find(megathrust_full.xc(:,2)>(2.98e5-.7*megathrust_full.xc(:,1)));
% 
%     vals = megathrust_full.xc(:,1)*0;
%     vals(Idrop)=1;
%     %megathrust.plotPatch(vals)
%     keepVertices = megathrust_full.vertices;
%     keepVertices(Idrop,:) = [];
% 
%     fileID = fopen('Slab2/japan_trench_xyz_crop.tri','w');
%     for i=1:length(keepVertices)
%         fprintf(fileID,'%d %d %d %d 0\n',i,keepVertices(i,1), keepVertices(i,2), keepVertices(i,3));
%     end
%     fclose(fileID);
% 
%     !cp Slab2/japan_trench_xyz.ned Slab2/japan_trench_xyz_crop.ned
% 
%     fltname='Slab2/japan_trench_xyz_crop';
%     megathrust=unicycle.geometry.triangleReceiver(fltname,unicycle.greens.nikkhoo15(mu,nu));
% 
% else
    megathrust = megathrust_full;
% end
%%

figure(3), clf
megathrust.plotPatch(megathrust.xc(:,3))
view(2)
daspect([1 1 1]), hold on
colorbar

%%





% figure(1), clf   
% 
% quiver(megathrust.xc(:,1),megathrust.xc(:,2),megathrust.sv(:,1),megathrust.sv(:,2),'k'), hold on
%     
%%

load coastlines.mat
In=coastlat>32 & coastlat <45 & coastlon > 138 & coastlon < 148 | isnan(coastlat);
[coastXkm,coastYkm]=latlon_to_xy_polyconic(coastlat(In),coastlon(In),lat0,lon0);

%% load Japan GPS and convert to xy
% data: Kreemer et al., 2014 (?)
% projection center: (143,39).
% Use polyconic projection

figure(2), clf
megathrust.plotPatch()
view(2)
daspect([1 1 1]), hold on
%

% gps=load('gps_data/GPS_OK.dat');
% gps = deuplicate_GSRM_gps_data(gps);
% save('gps_data/GPS_OK_dedup.dat','outgps','-ASCII');
gps=load('gps_data/GPS_OK_dedup.dat');


x=linspace(135,145,100);
y=lat0 + 8.5 + 2*(x - lon0 );
[lineX,lineY]=latlon_to_xy_polyconic(y,x,lat0,lon0);
plot(lineX*1e3,lineY*1e3,'k'), hold on

if crop_sz
    In=gps(:,1)<147 & gps(:,1)>138.5 & gps(:,2)+0.2*gps(:,1)>64.3 & gps(:,2)<41.9;
else
    In=gps(:,1)<147 & gps(:,1)>138.5 & gps(:,2)+0.2*gps(:,1)>64.3 & gps(:,2)+0.7*gps(:,1)< 143.9 & gps(:,2) < lat0 + 8.5 + 2*(gps(:,1) - lon0);
    
end
gpslat=gps(In,2);
gpslon=gps(In,1);
gpsE=gps(In,3);
gpsN=gps(In,4);
gpsEsig=gps(In,5);
gpsNsig=gps(In,6);
[gpsX,gpsY]=latlon_to_xy_polyconic(gpslat,gpslon,lat0,lon0);

%quiver(gpsX*1e3,gpsY*1e3,gpsE,gpsN)
%figure(1),clf, hold on

Iseafloor=gpslon>142 & gpslat<39;

quiver(gpsX(~Iseafloor)*1e3,gpsY(~Iseafloor)*1e3,gpsE(~Iseafloor),gpsN(~Iseafloor))

plot(coastXkm*1e3,coastYkm*1e3);

if crop_sz
    fname_all='gps_data/all_gps_tohoku_backslip_crop.dat';
    fname_allcoords='gps_data/all_gps_coords_crop.dat';
    fname_land='gps_data/land_gps_tohoku_backslip_crop.dat';
    fname_landcoords = 'gps_data/land_gps_coords_crop.dat';
else
    fname_all='gps_data/all_gps_tohoku_backslip.dat';
    fname_allcoords='gps_data/all_gps_coords.dat';
    fname_land='gps_data/land_gps_tohoku_backslip.dat';
    fname_landcoords = 'gps_data/land_gps_coords.dat';
end

% [lon,lat,hgt,uE,uN,uZ,sigE,sigN,sigZ,[wgt]]
all_backslip=[gpsX*1e3, gpsY*1e3, 0*gpsY, gpsE, gpsN, 0*gpsN, gpsEsig, gpsNsig, 0*gpsNsig];
land_backslip=all_backslip(~Iseafloor,:);

save(fname_all,'all_backslip', '-ASCII');
save(fname_land,'land_backslip', '-ASCII');

ylim([-6e5,8e5])

gpslonlat_land=[gpslon(~Iseafloor),gpslat(~Iseafloor)];
save(fname_landcoords,'gpslonlat_land','-ASCII')

gpslonlat_all=[gpslon,gpslat];
save(fname_allcoords,'gpslonlat_all','-ASCII')

gps_gmt = [gpslon(~Iseafloor),gpslat(~Iseafloor),gpsE(~Iseafloor),gpsN(~Iseafloor),gpsEsig(~Iseafloor),gpsNsig(~Iseafloor),0*gpsNsig(~Iseafloor)];
save('gps_data/all_gps_tohoku_backslip.gmt','gps_gmt','-ASCII')

%% create the plate motion model "rake file"


% get imposed rake
poles=readtable('gps_data/kreemer_pole_neg.txt');
%
omega_block = euler_vector(poles.Latp(1),poles.Lonp(1),poles.Omega(1)); % PA-OK

% flip sign
omega = - omega_block;

% get fault lat/lon
[fltlat,fltlon] = xy_to_latlon_polyconic(megathrust.xc(:,1)/1e3, megathrust.xc(:,2)/1e3, lon0, lat0);
Vrake_en = rot_matrix(fltlat,fltlon)*omega; % interlaced by station
VrakeE = Vrake_en(1:2:end);
VrakeN = Vrake_en(2:2:end);


%
figure(7),clf
plot3(fltlat,fltlon,sqrt(VrakeE.^2+VrakeN.^2),'k.'), hold on
%xlim([-8.5,5])
set(gca,'xdir','rev')
view([0,-1,0])
%

Vrake_ss = VrakeE .* megathrust.sv(:,1) + VrakeN .* megathrust.sv(:,2);
Vrake_ds = -VrakeE .* megathrust.sv(:,2) + VrakeN .* megathrust.sv(:,1);

rake_angle=rad2deg(atan2(Vrake_ds,Vrake_ss));

figure(5), clf
megathrust.plotPatch(sqrt(VrakeE.^2+VrakeN.^2)), hold on
megathrust.plotSlipVectors(Vrake_ss,Vrake_ds,400,'r')
colorbar
view(2)
title('Predicted rate (PA-OK Euler Pole, Kreemer et al. 2014)')
daspect([1 1 1])

if crop_sz
    fname_rake = 'Slab2/rake_PA_OK_kreemer_crop.dat';
else
    fname_rake = 'Slab2/rake_PA_OK_kreemer.dat';
end

rake_out=[rake_angle,Vrake_ss,Vrake_ds,VrakeE,VrakeN];
save(fname_rake,'rake_out','-ASCII')

Iplot=1;
mindist=50e3;
for i=1:length(VrakeE)
    dist_i = get_min_distance(megathrust.xc(i,1),megathrust.xc(i,2),megathrust.xc(Iplot,1),megathrust.xc(Iplot,2));
    if dist_i > mindist
        Iplot = [Iplot, i];
    end
end


slipvectors_out=[megathrust_geo.xc(Iplot,1), megathrust_geo.xc(Iplot,2), VrakeE(Iplot), VrakeN(Iplot) ];
save('Slab2/japan_trench_kreemer_slipvectors.dat','slipvectors_out','-ASCII')

% save depth for plotting: positive downward for ease of GMT
save_geom_for_GMT('Slab2/japan_trench_depths_polygons_geo.dat', megathrust, -megathrust.xc(:,3)/1e3, lat0, lon0) 
    
% save plate rate for plotting: interp/no interp
save_geom_for_interp_GMT('Slab2/japan_trench_rates_geo.dat', megathrust, sqrt(Vrake_ss.^2 + Vrake_ds.^2), lat0, lon0) 
save_geom_for_GMT('Slab2/japan_trench_rates_polygons_geo.dat', megathrust, sqrt(Vrake_ss.^2 + Vrake_ds.^2), lat0, lon0) 





