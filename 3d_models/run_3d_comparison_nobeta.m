% run 3D alpha/beta optimization with the full dataset
% Eric, May 2019

set_jointinv_path()
clear all, close all
%


% plot optimum model results

% Determined with cross-validation and a grid search over values
alpha_best=0.02;
beta_best=0;

% load data

expNumber=4;
scenario = Jointinv(expNumber);
scenario.run_setup();


patchx = scenario.sources{1}.geom.xc(:,1);
patchy = scenario.sources{1}.geom.xc(:,2);

% compute the driving (long-term) slip rate on the fault patches
Vpl = [100; 0; 0]; % include a zero z component for ease of dot products
Vds = - scenario.sources{1}.geom.dv * Vpl ./ cosd(scenario.sources{1}.geom.dip).^2;
Vss = scenario.sources{1}.geom.sv * Vpl;
Vtransl=[Vss; Vds];
%coupling = ( (scenario0_train.modelVector(1:length(patchx)).^2 + scenario0_train.modelVector(length(patchx)+1:2*length(patchx)).^2).^0.5) ./ Vds;

true_model=load('test_data_3d/test_true_modelVector.dat');
true_coupling =(Vds-true_model(:,4))./Vds;

% make plot

 % assign weight values
scenario.userParams.smoothingWeights{1}{1}=alpha_best;
scenario.userParams.smoothingWeights{1}{2}=beta_best;

% run the inversion on the full dataset
scenario.run_inversion()

coupling = ( scenario.modelVector(length(patchx)+1:2*length(patchx))) ./ Vds;

%
figure(1),clf
fs=14;

subplot(3,3,1)

scenario.sources{1}.geom.plotPatch(coupling), hold on, colorbar
title('Smoothing only')
colormap(gca,flipud(hot))
%daspect([1 1 1])
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
caxis([0 1])
set(gca,'fontsize',fs)
shading flat

surfacerate=scenario.designMatrix * scenario.modelVector;
truegps=load('test_data_3d/gps_synthetic.dat');
quiver(truegps(:,1),truegps(:,2),-truegps(:,4),-truegps(:,5))
quiver(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2), -surfacerate(1:2:end), -surfacerate(2:2:end))

subplot(3,3,4)

scenario.sources{1}.geom.plotPatch(coupling - true_coupling), hold on, colorbar
title('Difference from synthetic')
caxis([-1 1])
colormap(gca,bluewhitered)
%daspect([1 1 1])
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
set(gca,'fontsize',fs)
shading flat


chi2 = (scenario.predVector - scenario.dataVector)'*scenario.datasets{1}.covarianceMatrix*(scenario.predVector - scenario.dataVector)/length(scenario.dataVector)

outfname = 'test_data_3d/coupling_model_smoothonly_nobeta.dat';
save_jointinv_model_trench0(scenario, coupling, outfname)

outfname = 'test_data_3d/coupling_difference_smoothonly_nobeta.dat';
save_jointinv_model_trench0(scenario, coupling - true_coupling, outfname)


% plot stresses
stress_ss_ds = scenario.sources{1}.KK*scenario.modelVector;
stress_ss = stress_ss_ds(1:length(scenario.modelVector)/2);
stress_ds = stress_ss_ds(length(scenario.modelVector)/2+1:end);

subplot(3,3,7)
scenario.sources{1}.geom.plotPatch(-stress_ds), hold on
%plot(gpsX,gpsY,'k.')
colorbar
caxis([-40, 40])
title('Tau_d')
shading flat
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
set(gca,'fontsize',fs)
colormap(gca,bluewhitered)

outfname = 'test_data_3d/stress_smoothonly_nobeta.dat';
save_jointinv_model_trench0(scenario, -stress_ds, outfname)

%

expNumber=5;
scenario = Jointinv(expNumber);
scenario.run_setup();

 % assign weight values
scenario.userParams.smoothingWeights{1}{1}=alpha_best;
scenario.userParams.smoothingWeights{1}{2}=beta_best;

% run the inversion on the full dataset
scenario.run_inversion()
%
coupling = ( scenario.modelVector(length(patchx)+1:2*length(patchx))) ./ Vds;



subplot(3,3,2)

scenario.sources{1}.geom.plotPatch(coupling), hold on, colorbar
title('Smoothing + constraint')
colormap(gca,flipud(hot))
%daspect([1 1 1])
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
caxis([0 1])
set(gca,'fontsize',fs)
shading flat
 
surfacerate=scenario.designMatrix * scenario.modelVector;
truegps=load('test_data_3d/gps_synthetic.dat');
quiver(truegps(:,1),truegps(:,2),-truegps(:,4),-truegps(:,5))
quiver(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2), -surfacerate(1:2:end), -surfacerate(2:2:end))
%
subplot(3,3,5)

scenario.sources{1}.geom.plotPatch(coupling - true_coupling), hold on, colorbar
title('Difference from synthetic')
caxis([-1 1])
colormap(gca,bluewhitered)
%daspect([1 1 1])
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
set(gca,'fontsize',fs)
shading flat

%

chi2 = (scenario.predVector - scenario.dataVector)'*scenario.datasets{1}.covarianceMatrix*(scenario.predVector - scenario.dataVector)/length(scenario.dataVector)

outfname = 'test_data_3d/coupling_model_constrainedsmooth_nobeta.dat';
save_jointinv_model_trench0(scenario, coupling, outfname)

outfname = 'test_data_3d/coupling_difference_constrainedsmooth_nobeta.dat';
save_jointinv_model_trench0(scenario, coupling - true_coupling, outfname)


% plot stresses
stress_ss_ds = scenario.sources{1}.KK*scenario.modelVector;
stress_ss = stress_ss_ds(1:length(scenario.modelVector)/2);
stress_ds = stress_ss_ds(length(scenario.modelVector)/2+1:end);

subplot(3,3,8)
scenario.sources{1}.geom.plotPatch(-stress_ds), hold on
%plot(gpsX,gpsY,'k.')
colorbar
caxis([-40, 40])
title('Tau_d')
shading flat
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
set(gca,'fontsize',fs)
colormap(gca,bluewhitered)

outfname = 'test_data_3d/stress_constrainedsmooth_nobeta.dat';
save_jointinv_model_trench0(scenario, -stress_ds, outfname)


%
expNumber=5;
scenario = Jointinv(expNumber);
scenario.run_setup();

 % assign weight values
scenario.userParams.smoothingWeights{1}{1}=0*alpha_best;
scenario.userParams.smoothingWeights{1}{2}=beta_best;

% run the inversion on the full dataset
scenario.run_inversion()
%
coupling = ( scenario.modelVector(length(patchx)+1:2*length(patchx))) ./ Vds;

%
subplot(3,3,3)

scenario.sources{1}.geom.plotPatch(coupling), hold on, colorbar
title('Constraints only')
colormap(gca,flipud(hot))
%daspect([1 1 1])
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
caxis([0 1])
set(gca,'fontsize',fs)
shading flat
 
surfacerate=scenario.designMatrix * scenario.modelVector;
truegps=load('test_data_3d/gps_synthetic.dat');
quiver(truegps(:,1),truegps(:,2),-truegps(:,4),-truegps(:,5))
quiver(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2), -surfacerate(1:2:end), -surfacerate(2:2:end))

subplot(3,3,6)

scenario.sources{1}.geom.plotPatch(coupling - true_coupling), hold on, colorbar
title('Difference from synthetic')
caxis([-1 1])
colormap(gca,bluewhitered)
%daspect([1 1 1])
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
set(gca,'fontsize',fs)
shading flat

chi2 = (scenario.predVector - scenario.dataVector)'*scenario.datasets{1}.covarianceMatrix*(scenario.predVector - scenario.dataVector)/length(scenario.dataVector)

outfname = 'test_data_3d/coupling_model_constrainedonly_nobeta.dat';
save_jointinv_model_trench0(scenario, coupling, outfname)

outfname = 'test_data_3d/coupling_difference_constrainedonly_nobeta.dat';
save_jointinv_model_trench0(scenario, coupling - true_coupling, outfname)


% plot stresses
stress_ss_ds = scenario.sources{1}.KK*scenario.modelVector;
stress_ss = stress_ss_ds(1:length(scenario.modelVector)/2);
stress_ds = stress_ss_ds(length(scenario.modelVector)/2+1:end);

subplot(3,3,9)
scenario.sources{1}.geom.plotPatch(-stress_ds), hold on
%plot(gpsX,gpsY,'k.')
colorbar
caxis([-40, 40])
title('Tau_d')
shading flat
xlim([0,2.6e5])
ylim([-2.6e5,2.6e5])
set(gca,'fontsize',fs)
colormap(gca,bluewhitered)

outfname = 'test_data_3d/stress_constrainedonly_nobeta.dat';
save_jointinv_model_trench0(scenario, -stress_ds, outfname)


