% create synthetic 3D forward model 
clear all, close all

% set jointinv and unicycle paths automatically
% run this cell from inside the 'example_nepal_earthquake' folder.
here = pwd;
cd ../functions
set_jointinv_path()
cd(here)

%%
% write Unicycle patch file
xc=0;
yc=0;
dip = 20;
strike=0;
fault_width = 300*1e3; % meters, along-dip width
fault_length = 600*1e3; % meters, along-dip width
fault_depth = 0; % top of the fault edge
nW = 20;
nL = 30;
patchfname='faults/ramp_3d.seg';
write_3d_ramp(patchfname, xc, yc, strike, dip, fault_width, fault_length, fault_depth, nW, nL)
%% create fault object
% use shear modulus of 30 GPa, poisson's ratio of 0.25
earthModel = unicycle.greens.okada92(30e3,0.25);
geom = unicycle.geometry.source(patchfname, earthModel);

% plot fault object
figure(1),clf, hold on
geom.plotPatch(geom.xc(:,3));
colorbar
title('Fault segments colored by depth in meters')
box on
view(37.5,37.5)
daspect([1,1,1])

% generate synthetic GPS site locations
Ngps = 60;
gpsx=linspace(-150e3,450e3,Ngps)';
gpsy=repmat(2e3*[-125,125,-75,75,-25,25]',10,1);
gpscoords=[gpsx, gpsy 0*gpsx];

plot(gpscoords(:,1),gpscoords(:,2),'r^')

%% create a simple slip distribution

slipcenterx = 50e3;
slipcentery = 100e3;
slipwidth = 75e3;
sliplength = 175e3;
slipmaximum = 3;

% equation for top half of an ellipsoid
slipmag = real( slipmaximum * sqrt( 1 - (geom.xc(:,1) - slipcenterx).^2/slipwidth^2 - (geom.xc(:,2) - slipcentery).^2/sliplength^2 ) );

% save true model
save('./faults/synthetic_slip_model','slipmag');


%% get the greens function matrix and predict surface displacements

slipComponents = 2; % (keep dip component of G only)
kernelFolder = 'faults'; % place to save/load the displacement kernel file after calculating it
G = unicycle_displacement_kernel(geom,gpscoords,slipComponents,kernelFolder);

figure(3),clf
pcolor(G),shading flat
set(gca,'ydir','rev')
caxis([-0.1 0.1])
colormap(bluewhitered)
colorbar
title('Greens function matrix for synthetic example')
%daspect([1 1 1])
xlabel('model columns')
ylabel('data rows')
box on


% add synthetic noise
sigma = 0.01;

d_obs=G*slipmag + sigma*randn(Ngps*3,1);
obs_E = d_obs(1:3:end);
obs_N = d_obs(2:3:end);
obs_U = d_obs(3:3:end);

figure(2),clf
geom.plotPatch(slipmag); shading flat, hold on
geom.plotPatch()
scatter(gpscoords(:,1),gpscoords(:,2),200,obs_U,'filled','MarkerEdgeColor',[0 0 0 ])
scaled_quiver(gpscoords(:,1),gpscoords(:,2),obs_E,obs_N,2e5,{'k','linewidth',2})
view(2)
caxis([-3,3])
colormap(bluewhitered)
colorbar
title('input model')

view(37.5,37.5)
daspect([1,1,1])

% save synthetic GPS data file: format [E N U vE vN vU sigE sigN sigU]
gps_out = [gpscoords, obs_E, obs_N, obs_U, sigma.*ones(Ngps,3)];
save('data/synthetic_gps.dat','gps_out','-ASCII');



%% compute a naiive, unsmoothed model
naiivemodel = G \ d_obs;

d_pred  = G*naiivemodel;
pred_E = d_pred(1:3:end);
pred_N = d_pred(2:3:end);
pred_U = d_pred(3:3:end);

figure(5),clf
geom.plotPatch(naiivemodel); shading flat, hold on
geom.plotPatch()
scatter(gpscoords(:,1),gpscoords(:,2),300,obs_U,'filled','MarkerEdgeColor',[0 0 0 ])
scatter(gpscoords(:,1),gpscoords(:,2),100,pred_U,'filled','MarkerEdgeColor',[1 0 0 ])
scaled_quiver(gpscoords(:,1),gpscoords(:,2),obs_E,obs_N,2e5,{'k','linewidth',2})
scaled_quiver(gpscoords(:,1),gpscoords(:,2),pred_E,pred_N,2e5,{'r','linewidth',1})
view(2)
%caxis([-3,3])
colormap(bluewhitered)
colorbar
title('naiive model')

view(37.5,37.5)
daspect([1,1,1])


%% create a smoothing matrix

L = 1e8*compute_laplacian(geom.xc(:,1),geom.xc(:,2),geom.xc(:,3),4);

figure(6),clf
imagesc(L),shading flat
set(gca,'ydir','rev')
%caxis([-2 2])
colormap(bluewhitered)
colorbar
title('Laplacian matrix for synthetic example')
daspect([1 1 1]) 

%% compute a better model

lambda=0.3;
G_aug = [G; lambda*L];
d_aug = [d_obs; zeros(geom.N,1)];

bettermodel = lsqnonneg(G_aug,d_aug);


d_pred  = G*bettermodel;
pred_E = d_pred(1:3:end);
pred_N = d_pred(2:3:end);
pred_U = d_pred(3:3:end);

figure(7),clf
geom.plotPatch(bettermodel); shading flat, hold on
geom.plotPatch()
scatter(gpscoords(:,1),gpscoords(:,2),300,obs_U,'filled','MarkerEdgeColor',[0 0 0 ])
scatter(gpscoords(:,1),gpscoords(:,2),100,pred_U,'filled','MarkerEdgeColor',[1 0 0 ])
scaled_quiver(gpscoords(:,1),gpscoords(:,2),obs_E,obs_N,2e5,{'k','linewidth',2})
scaled_quiver(gpscoords(:,1),gpscoords(:,2),pred_E,pred_N,2e5,{'r','linewidth',1})
view(2)
caxis([-3,3])
colormap(bluewhitered)
colorbar
title('better model')

view(37.5,37.5)
daspect([1,1,1])

