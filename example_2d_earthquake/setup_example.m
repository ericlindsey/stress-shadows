% create synthetic 2D forward model - 
clear all, close all

%%
% write Unicycle patch file
dip = 10;
v_plate = 1;
fault_width = 300*1e3; % meters, along-dip width
npatch = 50;
patchfname='faults/ramp_2d.seg';
write_2d_ramp(patchfname, dip, fault_width, npatch)

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
daspect([1,1e3,1])

%% generate slip model and synthetic GPS site locations
Ngps = 50;
gpsx=linspace(-150e3,450e3,Ngps)';
gpscoords=[gpsx, repmat(0*gpsx,1,2)];


% create a simple slip distribution

slipcenter = 49e3;
slipwidth = 75e3;
slipmaximum = 3;

% equation for top half of an ellipse
slipmag = real( (slipmaximum/slipwidth)*sqrt( slipwidth.^2 - (geom.xc(:,1) - slipcenter).^2 ) );

figure(2),clf,hold on
geom.plotPatch(slipmag); shading flat
colorbar
title('Synthetic slip model')
box on
view(37.5,37.5)
daspect([1,2e3,1])

plot(gpscoords(:,1),gpscoords(:,2),'k^')


%
% % save true model
% save('faults/synthetic_slip_model','slipmag');


% get the greens function matrix and predict surface displacements

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
daspect([1 1 1])
xlabel('model')
ylabel('data')

% add synthetic noise
sigma = 0.1;

d_obs=G*slipmag + sigma*randn(Ngps*3,1);
obs_E = d_obs(1:3:end);
obs_N = d_obs(2:3:end);
obs_U = d_obs(3:3:end);

figure(4),clf, hold on
plot(gpscoords(:,1),obs_E,'-bx')
plot(gpscoords(:,1),obs_N,'-rx')
plot(gpscoords(:,1),obs_U,'-mx')
legend('East','North','Up')

% save synthetic GPS data file: format [E N U vE vN vU sigE sigN sigU]
gps_out = [gpscoords, obs_E, obs_N, obs_U, ones(Ngps,3)];
save('data/synthetic_gps.dat','gps_out','-ASCII');

%%

% compute a naiive, unsmoothed model
naiivemodel = lsqlin(G,d_obs);
% compute a naiive, unsmoothed model
%naiivemodel = G \ d_obs;
d_pred  = G*naiivemodel;
pred_E = d_pred(1:3:end);
pred_N = d_pred(2:3:end);
pred_U = d_pred(3:3:end);


plot(gpscoords(:,1),pred_E,'-c','linewidth',2)
plot(gpscoords(:,1),pred_N,'-m','linewidth',2)
plot(gpscoords(:,1),pred_U,'-k','linewidth',2)
legend('East','North','Up','fitE','fitN','fitU')

figure(5),clf, hold on
plot(geom.xc(:,1),slipmag,'-kx')
plot(geom.xc(:,1),naiivemodel,'r')
legend('truth','naiive model')

%ylim([-1 4])
%% create a smoothing matrix

L = compute_laplacian_1d(geom.N);

figure(6),clf
imagesc(L),shading flat
set(gca,'ydir','rev')
caxis([-2 2])
colormap(bluewhitered)
colorbar
title('Laplacian matrix for synthetic example')
daspect([1 1 1]) 

%% compute a better model

lambda=1;
G_aug = [G; lambda*L];
d_aug = [d_obs; zeros(geom.N,1)];

bettermodel = lsqlin(G_aug, d_aug);

d_pred  = G*bettermodel;
pred_E = d_pred(1:3:end);
pred_N = d_pred(2:3:end);
pred_U = d_pred(3:3:end);


figure(4),clf, hold on
plot(gpscoords(:,1),obs_E,'-b.')
plot(gpscoords(:,1),obs_N,'-r.')
plot(gpscoords(:,1),obs_U,'-m.')
plot(gpscoords(:,1),pred_E,'-c','linewidth',2)
plot(gpscoords(:,1),pred_N,'-m','linewidth',2)
plot(gpscoords(:,1),pred_U,'-k','linewidth',2)
legend('East','North','Up','predicted East','predicted North','predicted Up')


figure(5),clf, hold on
plot(geom.xc(:,1),slipmag,'-kx')
%plot(geom.xc(:,1),naiivemodel,'-r')
plot(geom.xc(:,1),bettermodel,'-b','linewidth',2)
legend('true slip','smoothed model')


%% compute an even better model, without negative slip.

lambda=1;
G_aug = [G; lambda*L];
d_aug = [d_obs; zeros(geom.N,1)];

bestmodel = lsqnonneg(G_aug,d_aug);

d_pred  = G*bestmodel;
pred_E = d_pred(1:3:end);
pred_N = d_pred(2:3:end);
pred_U = d_pred(3:3:end);


figure(6),clf, hold on
plot(gpscoords(:,1),obs_E,'-b.')
plot(gpscoords(:,1),obs_N,'-r.')
plot(gpscoords(:,1),obs_U,'-m.')
plot(gpscoords(:,1),pred_E,'-c','linewidth',2)
plot(gpscoords(:,1),pred_N,'-m','linewidth',2)
plot(gpscoords(:,1),pred_U,'-k','linewidth',2)
legend('East','North','Up','predicted East','predicted North','predicted Up')


figure(7),clf, hold on
plot(geom.xc(:,1),slipmag,'-kx')
%plot(geom.xc(:,1),naiivemodel,'-r')
plot(geom.xc(:,1),bettermodel,'-b','linewidth',2)
plot(geom.xc(:,1),bestmodel,'-r','linewidth',2)
legend('true slip','smoothed with non-negativity')

