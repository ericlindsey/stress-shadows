% run a 2D forward model - use the stress kernel (computed as smoothing
% matrix) to do the stress shadow computation
clear all, close all

set_jointinv_path()


% create ramp and points files


% write Unicycle patch file
dip = 10;
v_plate = 1;
fault_width  = 250*1e3; %meters - along-dip width
fault_length = 700*1e3; %meters - along-strike width
npatchw = 25;
npatchl = 35;
patch_width  = fault_width/npatchw;
patch_length = fault_length/npatchl;

% saved results: for dip 10, w 220, l 700,
% size  22x35   ( 1170 patches) took   65 sec. on Eric's desktop  File: ( 0.06 GB)
% size  22x35   (  770 patches) took   37 sec. on Eric's desktop  File: ( 0.03 GB) - without long patches
% size  22x70   ( 1540 patches) took  150 sec. on Eric's desktop  File: ( 0.11 GB) - without long patches
% size  44x120  ( 5280 patches) took 2023 sec. on Eric's desktop  File: ( 1.36 GB) - without long patches


% saved results: for dip 10, w 320, l 700,
% size  32x35   ( 1120 patches) took   67 sec. on Eric's desktop  File: ( 0.06 GB) - without long patches


% saved results: for dip 10, w 250, l 700,
% size  25x35   ( 875 patches) took   45 sec. on Eric's desktop  File: ( 0.03 GB) - without long patches
% size  50x70   (3500 patches) took  539 sec. on Eric's desktop  File: ( 0.6  GB) - without long patches


% % write Unicycle patch file
% dip = 10;
% v_plate = 1;
% fault_width  = 250*1e3; %meters - along-dip width
% fault_length = 1000*1e3; %meters - along-strike width
% npatchw = 100;
% npatchl = 100;
% patch_width  = fault_width/npatchw;
% patch_length = fault_length/npatchl;

% saved results: for dip 10, w 250, l 1000, plus extra long patches,
% size  25x20   (  900 patches) took   47 sec. on Eric's desktop  File: ( 0.04 GB)
% size  25x40   ( 1400 patches) took   89 sec. on Eric's desktop  File: ( 0.09 GB)
% size  25x100  (2,900 patches) took  292 sec. on Eric's desktop  File: ( 0.40 GB - deleted)
% size  50x50   (2,900 patches) took  292 sec. on Eric's desktop  File: ( 0.40 GB)
% size  50x100  (5,800 patches) took  936 sec. on Eric's desktop  File: ( 1.08 GB - deleted)
% size 100x100 (10,400 patches) took    x sec. on Eric's desktop  File: ( 6.5  GB - deleted)


 % write patch file
patchfname='test_data_3d/ramp_3d.seg';
fileID = fopen(patchfname,'w');
fprintf(fileID,'%s\n','#patch file generated automatically - 3D ramp model, constant patch size');
fprintf(fileID,'%s\n',...
    '# n  Vpl      x1              x2   x3  Length        Width   Strike  Dip  Rake      L0       W0        qL  qW');
fprintf(fileID,'%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %d %d\n',...
      1, v_plate, -fault_length/2, 0,    0, fault_length, fault_width, 0,   dip, 90, patch_length, patch_width, 1, 1);

% % add extra-long patches on the ends
% fprintf(fileID,'%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %d %d\n',...
%       1, v_plate, -8.5*fault_length, 0,    0, 8*fault_length, fault_width, 0,   dip, 90, 1*fault_length, patch_width, 1, 1);
% 
% fprintf(fileID,'%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %d %d\n',...
%       1, v_plate,  0.5*fault_length, 0,    0, 8*fault_length, fault_width, 0,   dip, 90, 1*fault_length, patch_width, 1, 1);

fclose(fileID);

%
% write GPS data file with locations only
gpsfname='test_data_3d/gps_synthetic.dat';

sigma=0;

nx=5;
ny=21;
xmin=150e3;
xmax=250e3;
ymin=-250e3;
ymax=250e3;
[gpsx,gpsy] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

%stagger the grid in x
xshift = ((xmax-xmin)/(nx-1))/2;
gpsx(1:2:end,:)=gpsx(1:2:end,:)+xshift/2;
gpsx(2:2:end,:)=gpsx(2:2:end,:)-xshift/2;
% 
% %stagger the grid in y
 yshift = 0;%((ymax-ymin)/(ny-1))/2;
 gpsy(:,1:2:end)=gpsy(:,1:2:end)+yshift/2;
 gpsy(:,2:2:end)=gpsy(:,2:2:end)-yshift/2;

figure(1),clf
plot(gpsx,gpsy,'b.')


daspect([1 1 1])

gpscoords=[gpsx(:), gpsy(:), repmat(0*gpsx(:),1,7)];
save(gpsfname,'gpscoords','-ASCII');

% write polygon file

poly=[0 -1e9; 0 1e9; 1e9 1e9; 1e9 -1e9; 0 -1e9];
polyfname='test_data_3d/polygon_simple_east.dat';
save(polyfname,'poly','-ASCII');

%
%clear all`

% create the master object, and load matrices
expNumber = 40;
scenario = Jointinv(expNumber);
scenario.run_setup();


%

% now assign some locked patches and directly compute the slip rates.
% Note, this is not done inside jointinv.
close all

figure(1),clf

patchx = scenario.sources{1}.geom.xc(:,1);
patchy = scenario.sources{1}.geom.xc(:,2);

randvals=rand(length(patchx),1);


% assign locked patches based on patch (center) coordinate
% lockstart1 = 90e3;
% lockend1   = 90e3;
% locky1     = 150e3;
% cond1 = (patchx > lockstart1  &  patchx < lockend1 & patchy>locky1);
% cond1a = ( (patchx - (lockstart1+lockend1)/2).^2 ./ ((lockstart1-lockend1)/2).^2  +  (patchy - locky1).^2 ./ ...
%     ((lockstart1-lockend1)/2).^2  < 1);

% 
% lockstart2 = 70e3;
% lockend2   = 150e3;
% locky2     = -300e3;
% end_radius=(lockstart2-lockend2)/2;
% cond1 = (patchx > lockstart2  &  patchx < lockend2 & patchy<locky2 );
% cond1a = ( (patchx - (lockstart2+lockend2)/2).^2 ./ end_radius.^2  +  (patchy - locky2).^2 ./ ...
%     (4*end_radius).^2  < 1);
% 
% lockx0     = 90e3;
% lockx_r0   = 40e3;
% locky0     = 200e3;
% locky_r0   = 150e3;
% cond3 = ( (patchx - lockx0).^2 ./ lockx_r0.^2  +  (patchy - locky0).^2 ./ locky_r0.^2  < 1+0.1*randvals);
% 
% lockx0     = 110e3;
% lockx_r0   = 40e3;
% locky0     = -1400e3;
% locky_r0   = 1290e3;
% cond2 = ( (patchx - lockx0).^2 ./ lockx_r0.^2  +  (patchy - locky0).^2 ./ locky_r0.^2  < 1+0.1*randvals);
% 
% lockx0     = 100e3;
% lockx_r0   = 50e3;
% locky0     = -250e3;
% locky_r0   = 200e3;
% cond2a = ( (patchx - lockx0).^2 ./ lockx_r0.^2  +  (patchy - locky0).^2 ./ locky_r0.^2  < 1+0.1*randvals);
%
%
lockx0     = 97.5e3;
lockx_r0   = 50e3;
locky0     = 0e3;
locky_r0   = 100e3;
cond2 = ( (patchx - lockx0).^2 ./ lockx_r0.^2  +  (patchy - locky0).^2 ./ locky_r0.^2  < 1+0*randvals);

creepstart = 250e3;
cond4 = patchx > creepstart;

Ilocked =  ( cond2 | cond4 ); 
%Ilocked =  ( cond1 | cond1a| cond2 | cond2a | cond3 | cond4 ); 


figure(1),clf
scenario.sources{1}.geom.plotPatch(1*Ilocked)
daspect([1 1 1])
xlim([0,2.5e5])
ylim([-5.25e5,5.25e5])
%

% compute the driving (long-term) slip rate on the fault patches
Vpl = [100; 0; 0]; % include a zero z component for ease of dot products
Vds = - scenario.sources{1}.geom.dv * Vpl ./ cosd(scenario.sources{1}.geom.dip).^2;
Vss = scenario.sources{1}.geom.sv * Vpl;
Vtransl=[Vss; Vds];

%in: 0 or 1 if on the moving plate
%in=scenario.datasets{1}.coordinates(:,1) > 0;
poly = load('test_data_3d/polygon_simple_east.dat');
in = inpolygon(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2),poly(:,1),poly(:,2));

% compute loading stress rate
sig0 = scenario.smoothingMatrix * Vtransl;


% force the deep slip rate
Icreep=patchx>creepstart;

strikesliprate = zeros(size(patchx));
dipsliprate = zeros(size(patchx));
dipsliprate(Icreep)=Vds(Icreep);
sliprate=[strikesliprate; dipsliprate];

deepcreepstressrate = -scenario.smoothingMatrix * sliprate;


% extract the stress kernel and stress rate for only free patches
subK = get_free_patches_4x(scenario.smoothingMatrix, ~Ilocked);
Ilocked_doubled = [Ilocked; Ilocked];
substressrate=sig0(~Ilocked_doubled) + deepcreepstressrate(~Ilocked_doubled); %note: for 3D case add strike


% use this subset of patches to balance the stress rate
subsliprate=subK\substressrate;


% compute the final slip rate
sliprate(~Ilocked_doubled)=subsliprate; %note: for 3D case add strike

%
figure(2)

subplot(3,1,1)

% what we "actually observe"
% surfacerate=scenario.designMatrix*(-Vtransl + sliprate);
% flipping both signs 

surfacerate=scenario.designMatrix*(Vtransl - sliprate) + sigma * randn(length(scenario.datasets{1}.coordinates(:,1))*3,1);

plot(scenario.datasets{1}.coordinates(:,1)/1e3, surfacerate(1:3:end), 'b.'),hold on
plot(scenario.datasets{1}.coordinates(:,1)/1e3, surfacerate(2:3:end), 'k.'),hold on
plot(scenario.datasets{1}.coordinates(:,1)/1e3, surfacerate(3:3:end), 'r.'),hold on
ylabel('surface velocity (mm/yr)')

subplot(3,1,2)
plot(patchx/1e3, sliprate(1:length(patchx)),'b.'), hold on
ylabel('strike slip rate (mm/yr)')

subplot(3,1,3)
plot(patchx/1e3, sliprate(length(patchx)+1:2*length(patchx)),'b.'), hold on
ylabel('dip slip rate (mm/yr)')
xlabel('distance from fault (km)')


%
%close all

%figure('pos',[1200 800 720 550]), clf, hold on
figure(3),clf,hold on
fs=14;

subplot(1,2,1)
lock=Ilocked*1 + Icreep*1;
lock(lock==0)=nan;
scenario.sources{1}.geom.plotPatch(1-lock), hold on
plot(gpsx,gpsy,'b.')
quiver(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2), -surfacerate(1:3:end), -surfacerate(2:3:end))
%quiver(scenario.datasets{1}.coordinates(:,1),scenario.datasets{1}.coordinates(:,2), surfacerate(1:3:end),  surfacerate(2:3:end))
%title('locked patches, and vectors')
title('Locking')
daspect([1 2 1])
xlim([0,3.5e5])
ylim([-5.1e5,5.1e5])
shading flat
% set(gca,'fontsize',fs)
% xlabel('Distance East (km)')
% ylabel('Distance North (km)')
% set(gca,'XTickLabel',[0 50 100 150 200 250] )
% set(gca,'YTickLabel',[-500 -400 -300 -200 -100 0 100 200 300 400 500] )
% set(gca,'Position',[0.1   0.1200    0.32    0.8150]);
% box on
daspect([1 1 1])

Inan=Ilocked*1;
Inan(Ilocked)=nan;

% subplot(1,4,2)
% scenario.sources{1}.geom.plotPatch(sliprate(1:length(patchx)) + Inan), hold on, colorbar
% title('strike slip rate')
% %daspect([1 1 1])
% xlim([0,2.5e5])
% ylim([-5.1e5,5.1e5])
% caxis([-20,20])
% shading flat
% set(gca,'fontsize',fs)
% %
% subplot(1,4,3)
% scenario.sources{1}.geom.plotPatch(sliprate(length(patchx)+1:2*length(patchx))+ Inan), hold on, colorbar
% title('dip slip rate')
% %daspect([1 1 1])
% xlim([0,2.5e5])
% ylim([-5.1e5,5.1e5])
% shading flat
% set(gca,'fontsize',fs)

subplot(1,2,2)

coupling = (Vds - (sliprate(length(patchx)+1:2*length(patchx)).^2).^0.5) ./ Vds;

scenario.sources{1}.geom.plotPatch(coupling), hold on, colorbar
title('Coupling')
colormap(flipud(hot))
daspect([1 1 1])
xlim([0,2.5e5])
ylim([-5.1e5,5.1e5])
shading flat

% set(gca,'fontsize',fs)
% xlabel('Distance East (km)')
% ylabel('Distance North (km)')
% set(gca,'XTickLabel',[0 50 100 150 200 250] )
% set(gca,'YTickLabel',[-500 -400 -300 -200 -100 0 100 200 300 400 500] )
%  set(gca,'Position',[0.48    0.1200    0.42    0.8150]);
% box on

%saveas2('forward_model_3d.pdf','-r200')

% save output gps data
gpsfname='test_data_3d/gps_synthetic.dat';
% columns: E N U vE vN vU sigE sigN sigU
gpscoords=[scenario.datasets{1}.coordinates , surfacerate(1:3:end), surfacerate(2:3:end), surfacerate(3:3:end), ...
    sigma*ones(size(scenario.datasets{1}.coordinates))];
save(gpsfname,'gpscoords','-ASCII');


% %% save "true" model
% 
 truemodelfname='test_data_3d/test_true_modelVector.dat';
 truemodel = [patchx,patchy,sliprate(1:length(patchx)),sliprate(length(patchx)+1:2*length(patchx))];
 save(truemodelfname,'truemodel','-ASCII'); 
 
%
Ix0=find(scenario.sources{1}.geom.xc(:,1)==min(scenario.sources{1}.geom.xc(:,1)));

% save coupling ratio
truecouplingfname='test_data_3d/test_true_coupling_ratio.dat';
truecoupling=[scenario.sources{1}.geom.xc(:,1)/1e3,scenario.sources{1}.geom.xc(:,2)/1e3,coupling];
%append duplicate values along x=0
truecoupling=[truecoupling; 0*Ix0, scenario.sources{1}.geom.xc(Ix0,2)/1e3,coupling(Ix0)];
save(truecouplingfname,'truecoupling','-ASCII'); 

% save fault depth
depthfname='test_data_3d/test_fault_depth.dat';
depthout=[scenario.sources{1}.geom.xc(:,1)/1e3,scenario.sources{1}.geom.xc(:,2)/1e3,scenario.sources{1}.geom.xc(:,3).^2/1e3];

depthout=[depthout; 0*Ix0, scenario.sources{1}.geom.xc(Ix0,2)/1e3,scenario.sources{1}.geom.xc(Ix0,3).^2/1e3];

save(depthfname,'depthout','-ASCII'); 


%
figure,clf
scatter3(scenario.sources{1}.geom.xc(:,1),scenario.sources{1}.geom.xc(:,2),1-coupling)
ylim([-5.1e5,5.1e5])
%xlim([0,1e4])
daspect([1 1 1e-5])



% plot stresses
stress_ss_ds = scenario.sources{1}.KK*(100/cosd(10) - sliprate);
stress_ss = stress_ss_ds(1:length(sliprate)/2);
stress_ds = stress_ss_ds(length(sliprate)/2+1:end);

figure,clf
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

outfname = 'test_data_3d/stress_synthetic.dat';
save_jointinv_model_trench0(scenario, -stress_ds, outfname)



