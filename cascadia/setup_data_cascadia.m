% Create necessary jointinv2 input files for Japan backslip model

clear, close all

run('../functions/set_jointinv_path')
unicyclepath='/Users/elindsey/Dropbox/code/geodesy/unicycle';
set_unicycle_path(unicyclepath)
import unicycle.*

%%
% set projection center
lon0=237;
lat0=45;

% convert Cacadia trench megathrust geometry to local XYZ
%!cp faults/cascadia_100km_geo.tri faults/cascadia_100km_xyz.tri

flt_geo=load('faults/cascadia_100km_geo.ned');

%% extract polygons and coastline

figure(4),clf
plot3(flt_geo(:,3),flt_geo(:,2),-flt_geo(:,4),'k.'),hold on
daspect([1 1 100])
Izero = find(flt_geo(:,4)==0);
plot(flt_geo(Izero,3),flt_geo(Izero,2),'rx')
% sort the points by latitude
[tracelat,Isort]=sort(flt_geo(Izero,2),'ascend');
tracelon=flt_geo(Izero,3);
tracelon=tracelon(Isort);
% save fault trace
flt_trace = [tracelon,tracelat];
save('faults/cascadia_trace.dat','flt_trace','-ASCII');

I100 = find(flt_geo(:,4)==100);
plot3(flt_geo(I100,3),flt_geo(I100,2),-flt_geo(I100,4),'gx')
% sort the points reversed by latitude
[deeptracelat,Isort]=sort(flt_geo(I100,2),'descend');
deeptracelon=flt_geo(I100,3);
deeptracelon=deeptracelon(Isort);
flt_polygon = [tracelon, tracelat; deeptracelon, deeptracelat; tracelon(1), tracelat(1)];
save('faults/cascadia_100km_polygon.dat','flt_polygon','-ASCII');

plot(flt_polygon(:,1),flt_polygon(:,2),'-b')
view(2)

%%

[fltXkm,fltYkm]=latlon_to_xy_polyconic(flt_geo(:,2),flt_geo(:,3),lat0,lon0);
%
fileID = fopen('faults/cascadia_100km_xyz.ned','w');
for i=1:length(fltXkm)
    fprintf(fileID,'%d %.9f %.9f %.9f\n',flt_geo(i,1), fltYkm(i)*1e3, fltXkm(i)*1e3, flt_geo(i,4)*1e3);
end
fclose(fileID);

% load, test and plot geometry
mu=30e9; nu=1/4;
fltname='faults/cascadia_100km_xyz';
megathrust=unicycle.geometry.triangleReceiver(fltname,unicycle.greens.nikkhoo15(mu,nu));

fltname_geo='faults/cascadia_100km_geo';
megathrust_geo=unicycle.geometry.triangleReceiver(fltname_geo,unicycle.greens.nikkhoo15(mu,nu));


% drop the deep values
if false
    Idrop=find(megathrust.xc(:,3) < -80e3);
    vals = megathrust.xc(:,1)*0;
    vals(Idrop)=1;
    %megathrust.plotPatch(vals)
    keepVertices = megathrust.vertices;
    keepVertices(Idrop,:) = [];

    fileID = fopen('faults/cascadia_80km_xyz.tri','w');
    for i=1:length(keepVertices)
        fprintf(fileID,'%d %d %d %d 0\n',i,keepVertices(i,1), keepVertices(i,2), keepVertices(i,3));
    end
    fclose(fileID);
    !cp faults/cascadia_100km_xyz.ned faults/cascadia_80km_xyz.ned
    fltname='faults/cascadia_80km_xyz';
    megathrust=unicycle.geometry.triangleReceiver(fltname,unicycle.greens.nikkhoo15(mu,nu));
end



% IMPORTANT: check sign consistency of the megathrust triangles
check_triangleReceiver_signs(megathrust)

%
figure(1), clf   
megathrust.plotPatch(megathrust.xc(:,3)), hold on
megathrust.plotPatch(), hold on

%
%scaled_quiver(megathrust.xc(:,1),megathrust.xc(:,2),megathrust.dv(:,1),megathrust.dv(:,2),2e4,{'k'}), hold on
view(2)
%
daspect([1 1 1])
set(gcf, 'position',[1800 1500 600 1200])
%
load coastlines.mat
In=coastlat>39 & coastlat <51 & coastlon > -130 & coastlon < -115 | isnan(coastlat);
[coastXkm,coastYkm]=latlon_to_xy_polyconic(coastlat(In),coastlon(In),lat0,lon0-360);


plot(coastXkm*1e3,coastYkm*1e3,'g'), hold on

% load Cascadia GPS, select region, and convert to xy
% Use polyconic projection

% data: Kreemer et al., 2014

% gps=load('gps/GPS_NA.dat');
% In=gps(:,1)<240 & gps(:,1)>230 & gps(:,2)<49.5 & gps(:,2)>41; %gps(:,2)+0.2*gps(:,1)>64.6 & gps(:,2)+ gps(:,1)<187.5 & gps(:,2)<44.8;
% gpslat=gps(In,2);
% gpslon=gps(In,1);
% gpsE_orig=gps(In,3);
% gpsN_orig=gps(In,4);
% gpsEsig=gps(In,5);
% gpsNsig=gps(In,6);
% [gpsX,gpsY]=latlon_to_xy_polyconic(gpslat,gpslon,lat0,lon0);
% 
% scaled_quiver(gpsX*1e3,gpsY*1e3,gpsE_orig,gpsN_orig,5e3,{'b'})

%% data: Michel et al., 2018

gps=readtable('gps/Secular_Velocities_for_Eric_TXT.txt');
In=gps.Var1<240 & gps.Var1>220 & gps.Var2<50.8 & gps.Var2>40.7 & gps.Var3<19 & gps.Var6<3 & gps.Var7<3; %gps(:,2)+0.2*gps(:,1)>64.6 & gps(:,2)+ gps(:,1)<187.5 & gps(:,2)<44.8;
gpslat=gps.Var2(In);
gpslon=gps.Var1(In);
gpsE=gps.Var3(In);
gpsN=gps.Var4(In);
gpsEsig=gps.Var6(In);
gpsNsig=gps.Var7(In);
[gpsX,gpsY]=latlon_to_xy_polyconic(gpslat,gpslon,lat0,lon0);
figure(1)
%scaled_quiver(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2),scenario.dataVector(1:2:end),scenario.dataVector(2:2:end),5e3,{'k'})
hold on
scaled_quiver(gpsX*1e3,gpsY*1e3,gpsE,gpsN,5e3,{'r'})
daspect([ 1 1 1])
% [X(m),Y(m),hgt(m),uE,uN,uZ,sigE,sigN,sigZ,[wgt]]
michel_backslip=[gpsX*1e3, gpsY*1e3, 0*gpsY, gpsE, gpsN, 0*gpsN, gpsEsig, gpsNsig, 0*gpsNsig];
save('gps/michel_gps.dat','michel_backslip', '-ASCII');

gpslonlat_michel=[gpslon,gpslat];
save('gps/michel_gps_coords.dat','gpslonlat_michel','-ASCII')

%% create the matching plate rake file from Schmalzle et al. 2014, table TS09

schmalzle_rates = readtable('gps/schmalzle_ts09_plate_rates_noheader.txt');
sch_lon = schmalzle_rates.Var3;
sch_lat = schmalzle_rates.Var4;
sch_depth = schmalzle_rates.Var5;
[schX,schY]=latlon_to_xy_polyconic(sch_lat,sch_lon,lat0,lon0);
sch_ve = schmalzle_rates.Var7;
sch_vn = schmalzle_rates.Var8;

%megathrust_ve = interp2(schX,schY,sch_ve, megathrust.xc(:,1), megathrust.xc(:,2))

quiv_sc = 300;

figure(3),clf
subplot(1,3,1)
plot(schX,schY,'.'),hold on
scatter3(schX*1e3,schY*1e3,-sch_depth*1e3,30,sch_ve,'.')
daspect([1 1 1])
colorbar
%scaled_quiver(schX,schY,sch_ve,sch_vn,quiv_sc,{'r'})
%plot(megathrust.xc(:,1)/1e3,megathrust.xc(:,2)/1e3,'kx')

megathrust.plotPatch(megathrust_ve)

%
ve_interp = scatteredInterpolant([schX,schY],sch_ve);
ve_interp.Method = 'natural';
ve_interp.ExtrapolationMethod = 'nearest';
megathrust_ve = ve_interp([megathrust.xc(:,1)/1e3,megathrust.xc(:,2)/1e3]);


subplot(1,3,2)
megathrust.plotPatch(megathrust_ve)
colorbar
view(2)


vn_interp = scatteredInterpolant([schX,schY],sch_vn);
vn_interp.Method = 'natural';
vn_interp.ExtrapolationMethod = 'nearest';
megathrust_vn = vn_interp([megathrust.xc(:,1)/1e3,megathrust.xc(:,2)/1e3]);

subplot(1,3,3)
megathrust.plotPatch(megathrust_vn)
colorbar
view(2)

subplot(1,3,1)
scaled_quiver(megathrust.xc(:,1),megathrust.xc(:,2),megathrust_ve,megathrust_vn,quiv_sc,{'k'})
scaled_quiver(schX*1e3,schY*1e3,sch_ve,sch_vn,quiv_sc,{'r'})


% rotate to strike and dip coordinates, using the megathrust strike vector
megathrust_vs = megathrust_ve .* megathrust.sv(:,1) + megathrust_vn .* megathrust.sv(:,2);
megathrust_vd = -megathrust_ve .* megathrust.sv(:,2) + megathrust_vn .* megathrust.sv(:,1);

% compute plate motion rake angle for each patch
sch_rake_angle=rad2deg(atan2(megathrust_vd,megathrust_vs));

% plot plate motion magnitude and direction in strike, dip coordinates
figure(5), clf
megathrust.plotPatch(sqrt(megathrust_ve.^2+megathrust_vn.^2)), hold on
megathrust.plotSlipVectors(megathrust_vs,megathrust_vd,400,'r')
colorbar
view(2)

figure(6), clf
scatter3(megathrust.xc(:,1),megathrust.xc(:,2),sqrt(megathrust_ve.^2+megathrust_vn.^2)), hold on


% save rake file output
sch_rake_out=[sch_rake_angle,megathrust_vs,megathrust_vd,megathrust_ve,megathrust_vn];
save('faults/cascadia_100km_rake_schmalzle.dat','sch_rake_out','-ASCII')




%%

figure(1),clf, hold on
% data: Mccaffrey et al., 2013
gps=readtable('gps/mccaffrey_2013_suppl_tableS1.txt');
In=gps.Var2<240 & gps.Var2>220 & gps.Var3<50 & gps.Var3>41.1 & gps.Var8<30 & gps.Var9<30; %gps(:,2)+0.2*gps(:,1)>64.6 & gps(:,2)+ gps(:,1)<187.5 & gps(:,2)<44.8;
gpslat=gps.Var3(In);
gpslon=gps.Var2(In);
gpsE_orig=gps.Var5(In);
gpsN_orig=gps.Var6(In);
gpsEsig=gps.Var8(In);
gpsNsig=gps.Var9(In);
[gpsX,gpsY]=latlon_to_xy_polyconic(gpslat,gpslon,lat0,lon0);
scaled_quiver(gpsX*1e3,gpsY*1e3,gpsE_orig,gpsN_orig,5e3,{'k'})


%
% %% add a block correction following the WWW model of Li et al.

% add variable euler pole correction for non-rigid motion
pole_li=readtable('gps/li_euler_pole.txt');
omega_li = euler_vector(pole_li.Latp(1),pole_li.Lonp(1),pole_li.Omega(1)); % pole reported by Li et al.
Vpl_en = rot_matrix(gpslat,gpslon)*omega_li; % interlaced by station
VplE = Vpl_en(1:2:end);
VplN = Vpl_en(2:2:end);

% apply full correction for sites South of 46N, and linearly decrease the
% correction to zero at 49N
corr_mag = ((49 - gpslat)/3).*(gpslat<=49 & gpslat>=46) + 1*(gpslat<46);
VplE_corr = corr_mag.*VplE;
VplN_corr = corr_mag.*VplN;

% corrected GPS
gpsE = gpsE_orig - VplE_corr;
gpsN = gpsN_orig - VplN_corr;

%plot corrected vectors
figure(1)
scaled_quiver(gpsX*1e3,gpsY*1e3,gpsE,gpsN,5e3,{'r'})
title('Fault mesh,depth, and uncorrected/corrected GPS vectors')

figure(10)
subplot(2,1,1)
plot(gpsY,gpsE_orig,'r.')
subplot(2,1,2)
plot(gpsY,gpsN_orig,'r.')

% save GPS data output

% [X(m),Y(m),hgt(m),uE,uN,uZ,sigE,sigN,sigZ,[wgt]]
all_backslip=[gpsX*1e3, gpsY*1e3, 0*gpsY, gpsE, gpsN, 0*gpsN, gpsEsig, gpsNsig, 0*gpsNsig];
save('gps/mccaffrey_gps_www_corrected.dat','all_backslip', '-ASCII');

gpslonlat_all=[gpslon,gpslat];
save('gps/mccaffrey_gps_coords.dat','gpslonlat_all','-ASCII')

all_gmt=[gpslon,gpslat, gpsE, gpsN, gpsEsig, gpsNsig, 0*gpsNsig];
save('gps/mccaffrey_gps_www_corrected_geo.gmt','all_gmt', '-ASCII');


% %% create the plate motion model "rake file"

% get imposed rake
poles_kreemer=readtable('gps/kreemer_euler_poles.txt'); %pa-na (1); pa-jf (2)
%
omega_na = euler_vector(poles_kreemer.Latp(1),poles_kreemer.Lonp(1),poles_kreemer.Omega(1)); % PA-NA
omega_jf = euler_vector(poles_kreemer.Latp(2),poles_kreemer.Lonp(2),poles_kreemer.Omega(2)); % PA-JF

% flip sign
%omega = - omega_block;

%  subtract vectors
omega_jf_na = omega_na - omega_jf;

% get fault lat/lon
[fltlat,fltlon] = xy_to_latlon_polyconic(megathrust.xc(:,1)/1e3, megathrust.xc(:,2)/1e3, lon0, lat0);
Vrake_en = rot_matrix(fltlat,fltlon)*omega_jf_na; % interlaced by station
VrakeE_kreemer = Vrake_en(1:2:end);
VrakeN_kreemer = Vrake_en(2:2:end);

% add the li et al. correction term
Vrake_en_li = rot_matrix(fltlat,fltlon)*omega_li; % interlaced by station
VrakeE_li = Vrake_en_li(1:2:end);
VrakeN_li = Vrake_en_li(2:2:end);

% apply full correction for patches South of 46N, and linearly decrease the
% correction to zero at 49N
corr_mag = ((49 - fltlat)/3).*(fltlat<=49 & fltlat>=46) + 1*(fltlat<46);
%note, correction now has the opposite sign compared to the GPS data
VrakeE = VrakeE_kreemer + corr_mag.*VrakeE_li;
VrakeN = VrakeN_kreemer + corr_mag.*VrakeN_li;

%plot plate motion in E, N coordinates
figure(2),clf
megathrust.plotPatch(), hold on
scaled_quiver(megathrust.xc(:,1),megathrust.xc(:,2),VrakeE,VrakeN,5e2,{'r'}), hold on
view(2)
daspect([1 1 1])
set(gcf, 'position',[1200 1500 600 1200])
title('Predicted Vpl in E,N coordinates')

% rotate to strike and dip coordinates, using the megathrust strike vector
Vrake_ss = VrakeE .* megathrust.sv(:,1) + VrakeN .* megathrust.sv(:,2);
Vrake_ds = -VrakeE .* megathrust.sv(:,2) + VrakeN .* megathrust.sv(:,1);

% compute plate motion rake angle for each patch
rake_angle=rad2deg(atan2(Vrake_ds,Vrake_ss));

% plot plate motion magnitude and direction in strike, dip coordinates
figure(5), clf
megathrust.plotPatch(sqrt(VrakeE.^2+VrakeN.^2)), hold on
megathrust.plotSlipVectors(Vrake_ss,Vrake_ds,400,'r')
colorbar
view(2)

title('Predicted Vpl in strike, dip coordinates')
daspect([1 1 1])
set(gcf, 'position',[600 1500 600 1200])

% save rake file output
rake_out=[rake_angle,Vrake_ss,Vrake_ds,VrakeE,VrakeN];
save('faults/cascadia_100km_rake_kreemer_li_corrected.dat','rake_out','-ASCII')

Iplot=1;
mindist=50e3;
for i=1:length(VrakeE)
    dist_i = get_min_distance(megathrust.xc(i,1),megathrust.xc(i,2),megathrust.xc(Iplot,1),megathrust.xc(Iplot,2));
    if dist_i > mindist
        Iplot = [Iplot, i];
    end
end


slipvectors_out=[megathrust_geo.xc(Iplot,1), megathrust_geo.xc(Iplot,2), VrakeE(Iplot), VrakeN(Iplot) ];
save('faults/cascadia_100km_slipvectors.dat','slipvectors_out','-ASCII')

% save depth for plotting: positive downward for ease of GMT
save_geom_for_GMT('faults/cascadia_100km_depths_polygons_geo.dat', megathrust, -megathrust.xc(:,3)/1e3, lat0, lon0) 
    
% save plate rate for plotting: interp/no interp
save_geom_for_interp_GMT('faults/cascadia_100km_rates_geo.dat', megathrust, sqrt(Vrake_ss.^2 + Vrake_ds.^2), lat0, lon0) 
save_geom_for_GMT('faults/cascadia_100km_rates_polygons_geo.dat', megathrust, sqrt(Vrake_ss.^2 + Vrake_ds.^2), lat0, lon0) 
