% run a 2D forward model - use the stress kernel (computed as smoothing matrix) to do the stress shadow computation
clear all

% setup paths if needed
%set_jointinv_path()

% write Unicycle patch file
dip = 10;
v_plate = 1;
fault_width = 312*1e3; %meters - along-dip width
npatch = 104;
patch_width = fault_width/npatch;
patchfname='test_data_2d/ramp_2d.seg';
fileID = fopen(patchfname,'w');
fprintf(fileID,'%s\n','#patch file generated automatically - 2D ramp model, constant patch size');
fprintf(fileID,'%s\n',        '# n  Vpl    x1      x2   x3   Length  Width   Strike  Dip  Rake      L0     W0    qL qW');
fprintf(fileID,'%s %.9f %s %.9f %s %.9f %s %.9f %s\n', '  1  ', v_plate, ...
    ' -50e7       0    0    100e7 ', fault_width, ' 0 ', dip, ' 90   100e7  ',patch_width,'  1  1 ');
fclose(fileID);

% synthetic noise level
sigma=0;

% generate GPS coordinates
gpsx=linspace(120e3,300e3,40)';
gpscoords=[gpsx, repmat(0*gpsx,1,8)];

% write polygon file
poly=[0 -1e9; 0 1e9; 1e9 1e9; 1e9 -1e9; 0 -1e9];
polyfname='test_data_2d/polygon_simple_east.dat';
save(polyfname,'poly','-ASCII');

% create the master object, and load matrices
expNumber = 0;
scenario = Jointinv(expNumber);
scenario.run_setup();

% now assign some locked patches and directly compute the slip rates.
% Note, this step is not done inside jointinv.

% assign locked patches based on patch (center) coordinate
%lockstart=0e3;
lockend=150e3;
creepstart=200e3;
lockstart=lockend/3;
patchx = scenario.sources{1}.geom.xc(:,1);
Ilocked =  (patchx > lockstart  &  patchx < lockend) | patchx > creepstart; 

% compute the driving (long-term) slip rate on the fault patches
Vpl = [100; 0; 0]; % include a zero z component for ease of dot products
Vds = -scenario.sources{1}.geom.dv * Vpl ./ cosd(scenario.sources{1}.geom.dip).^2;
%Vss = scenario.sources{1}.geom.sv * Vpl; % strike not used for 2D case
%in: 0 or 1 if on the moving plate
in=scenario.datasets{1}.coordinates(:,1) > 0;

% compute loading stress rate
% make sure to use only the stress kernel part of the smoothing matrix.
KK=scenario.smoothingMatrix(1:length(scenario.modelVector),1:length(scenario.modelVector));
sig0 = KK * Vds;

sliprate=zeros(size(sig0));
% set the deep slip rate and compute the stress caused by it
sliprate(patchx>creepstart)=Vds(patchx>creepstart);
deepcreepstressrate = -KK * sliprate;

% extract the stress kernel and stress rate for only free patches
subK = KK(~Ilocked,~Ilocked); % note, this will be more complicated in 3D
substressrate=sig0(~Ilocked) + deepcreepstressrate(~Ilocked); %note: for 3D case add strike

% use this subset of patches to balance the stress rate
subsliprate=subK\substressrate;

% compute the final slip rate
sliprate(~Ilocked)=subsliprate; %note: for 3D case add strike

% compute the GPS velocities
surfacerate=scenario.designMatrix * (Vds-sliprate) + sigma*randn(size(gpsx));

figure(1)
subplot(2,1,1)
plot(scenario.datasets{1}.coordinates(:,1), surfacerate-Vpl(1)*(1-heaviside(gpscoords(:,1))), '-b.'),hold on

subplot(2,1,2)
plot(scenario.sources{1}.geom.xc(:,1), sliprate,'-k.'), hold on

% save output gps data
%gpsfname='test_data_2d/gps_profile.dat'; % generated once, with one noise realization
gpsfname='test_data_2d/zero_noise_data.dat';

% columns: E N U vE vN vU sigE sigN sigU
% for this case, only columns 1,4,7 are nonzero.
gpscoords=[scenario.datasets{1}.coordinates , surfacerate, repmat(0*surfacerate,1,2), 0*surfacerate + 1, repmat(0*surfacerate,1,2)];
save(gpsfname,'gpscoords','-ASCII');

% save "true" model
truemodelfname='test_data_2d/test_true_modelVector.dat';
truemodel = [scenario.sources{1}.geom.xc(:,1), sliprate];
save(truemodelfname,'truemodel','-ASCII');    

% % save output gps data in another format
% gpsfname=['synthetic_gps_x',num2str(lockstart/1e3) , '.dat'];
% gpscoords=[scenario.datasets{1}.coordinates(:,1)/1e3 , surfacerate-Vpl(1)*(1-heaviside(gpscoords(:,1)))];
% save(gpsfname,'gpscoords','-ASCII');
% 
% truemodelfname=['synthetic_slip_x', num2str(lockstart/1e3), '.dat'];
% truemodel = [scenario.sources{1}.geom.xc(:,1)/1e3, sliprate];
% truemodel(1,1)=0;
% save(truemodelfname,'truemodel','-ASCII');  





