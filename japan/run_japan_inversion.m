
% run and plot inversions for Japan

 expNumber=12; % No stress constraints. With stress smoothing, rake coordinates, plate-rate bound on rake, small rake-perp bound, slip penalty on both. 12 sec
%expNumber=1; % With stress constraints and stress smoothing, rake coordinates, plate-rate bound on rake, small rake-perp bound, slip penalty on both. z sec

% create the master object, and load matrices
scenario = Jointinv(expNumber);
scenario.run_setup(); %high-res fault on desktop: displacement kernels 8sec, stress kernels 205 sec
%
% run the inversion
tic
scenario.run_inversion()
toc
%

% scenario.modelVector = 100+0*scenario.modelVector;
% scenario.predVector = scenario.designMatrix * scenario.modelVector;

%
% compute results components
results = calc_coupling_result_components(scenario);
disp(['chi2: ' num2str(results.chi2)])
%
% results.abic = abic_alphabeta(scenario);
% disp(['abic: ' num2str(results.abic)])

figoffset=3*expNumber;
plot_coupling_inversion(scenario,results,figoffset);

% load coastlines.mat
% Incoast=coastlat>32 & coastlat <45 & coastlon > 138 & coastlon < 148 | isnan(coastlat);
% figure(1)
% plot(coastlon(Incoast),coastlat(Incoast))

%%


%save_coupling_inversion(scenario,results);


%% plot laplacian values

figure(40),clf

ipatch = 1201;
%ipatch = 1011;

subplot(1,2,1)

scenario.sources{1}.geom.plotPatch(scenario.smoothingMatrix(ipatch,1:end/2))
scenario.sources{1}.geom.plotPatch()
view(2)
colorbar
%caxis([-3,3])
colormap(gca,bluewhitered)

subplot(1,2,2)

scenario.sources{1}.geom.plotPatch(scenario.smoothingMatrix(scenario.sources{1}.geom.N+ipatch,end/2+1:end))
scenario.sources{1}.geom.plotPatch()
view(2)
colorbar
%caxis([-3,3])
colormap(gca,bluewhitered)

