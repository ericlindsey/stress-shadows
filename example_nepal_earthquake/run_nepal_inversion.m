%% set jointinv and unicycle paths automatically
% run this cell from inside the 'example_nepal_earthquake' folder.
here = pwd;
cd ../functions
set_jointinv_path()
cd(here)


%% Run naiive GPS-only inversion

expNumber=401; 
scenario = Jointinv(expNumber);

% read the input file and set up matrices
scenario.run_setup();

% Run inversion with input parameters
tic
scenario.run_inversion()
toc

% plot the results
ax=figure(1); clf
Mw = plot_jointinv_slipmodel(scenario,ax,6e3);
plot_jointinv_dataset(scenario,ax,1e5, 200);
title(['inversion and data fits, magnitude ' num2str(round(Mw,2))])
xlim([-1e5,2e5])
ylim([-2e5,1e5])
daspect([1 1 1])

%% use ABIC to find the optimal smoothing weight

wgts = logspace(-4,-2,5); % logarithmically spaced from 1e-4 to 1e0.
abic = zeros(1,length(wgts));

for i=1:length(wgts)
    scenario.userParams.smoothingWeights{1}{1}=wgts(i);
    % run the inversion
    scenario.run_inversion()
    % get ABIC value
    abic(i) = abic_smoothingonly(scenario);
end

figure(10),clf
semilogx(wgts,abic-min(abic));
xlabel('smoothing weight')
ylabel('Delta ABIC, relative to minimum')

%% Run GPS-only inversion with a penalty on slip in low-resolution areas

expNumber=402; 
scenario = Jointinv(expNumber);

% read the input file and set up matrices
scenario.run_setup();

% Run inversion with input parameters
tic
scenario.run_inversion()
toc

% plot the results
ax=figure(2); clf
Mw = plot_jointinv_slipmodel(scenario,ax,6e3);
plot_jointinv_dataset(scenario,ax,1e5, 200);
title(['inversion and data fits, magnitude ' num2str(round(Mw,2))])
xlim([-1e5,2e5])
ylim([-2e5,1e5])
daspect([1 1 1])


%% Run InSAR only inversion, with smoothing and slip magnitude penalty


expNumber=403; 
scenario = Jointinv(expNumber);

% read the input file and set up matrices
scenario.run_setup();

%
% Run and plot inversion with input parameters
tic
scenario.run_inversion()
toc

%
figure(3); clf
subplot(1,2,1)
scenario.sources{1}.geom.plotPatch(scenario.modelVector(1:scenario.sources{1}.geom.N))
colorbar
title('rake-parallel slip')
view(2)
subplot(1,2,2)
scenario.sources{1}.geom.plotPatch(scenario.modelVector(scenario.sources{1}.geom.N+1:end))
colorbar
title('rake-perpendicular slip')
view(2)

figure(4); clf
ax=subplot(1,2,1)
plot_jointinv_slipmodel(scenario,ax,6e3);
caxis([-1 7])
colormap(bluewhitered)
%plot_jointinv_dataset(scenario,ax,1e5, 200);
title('model')
xlim([-1e5,2e5])
ylim([-2e5,1e5])
daspect([1 1 1])

ax=subplot(1,2,2)
plot_jointinv_slipmodel(scenario,ax,6e3);
caxis([-1 7])
colormap(bluewhitered)
plot_jointinv_dataset(scenario,ax,1e5, 200);
title('model and data fit')
xlim([0,1.5e5])
ylim([-0.9e5,0.4e5])
daspect([1 1 1])


%% Joint inversion

% an exercise for the reader!

% hint: in your exp_ file, you'll need to set, at a minimum:
% params.datasetTypes     = {'Static_GPS_Dataset','Static_LOS_Dataset'}; % cell list of dataset types (matching the corresponding Object filename)
% params.datasetFilenames = {'data/nepal_aria_gps_formatted.dat','data/varres_T048_insar.txt'}; % cell list of data files to read, 1-to-1 correspondence to the above list
      
% next, make sure all the other parameters are still compatible, and then
% try using ABIC to figure out the appropriate smoothing strength

