
% run and plot inversions for Cascadia

expNumber=10; % No stress constraints. With stress smoothing, rake coordinates, plate-rate bound on rake, small rake-perp bound, slip penalty on both. 12 sec
%expNumber=1; % With stress constraints and stress smoothing, rake coordinates, plate-rate bound on rake, small rake-perp bound, slip penalty on both. z sec

% create the master object, and load matrices
scenario = Jointinv(expNumber);
scenario.run_setup(); %high-res fault on desktop: displacement kernels 8sec, stress kernels 205 sec
%
% run the inversion
tic
scenario.run_inversion()
toc

% compute results components
results = calc_coupling_result_components(scenario);
disp(['chi2: ' num2str(results.chi2)])

% make default plots
figoffset=20;
plot_coupling_inversion(scenario,results, figoffset);

% save data
% save_coupling_inversion(scenario,results, 'best');

%%
figure(21)
load coastlines.mat
In=coastlat>39 & coastlat <51 & coastlon > -130 & coastlon < -115 | isnan(coastlat);
[coastXkm,coastYkm]=latlon_to_xy_polyconic(coastlat(In),coastlon(In),scenario.userParams.lat0,scenario.userParams.lon0-360);
plot(coastXkm*1e3,coastYkm*1e3,'b'), hold on


title(['No stress constraints, log_{10}\alpha = ' num2str(log10(scenario.userParams.smoothingWeights{1}{1})) ' , log_{10}\beta = ' num2str(log10(scenario.userParams.smoothingWeights{1}{2}))]);

set(gca, 'fontsize',14)

%%
figure(1)

plot(coastXkm*1e3,coastYkm*1e3,'b'), hold on

%
title(['Coupling - no stress constraints, \alpha = ' num2str(scenario.userParams.smoothingWeights{1}{1}) ' , beta = ' num2str(scenario.userParams.smoothingWeights{1}{2})]);

set(gca, 'fontsize',14)
