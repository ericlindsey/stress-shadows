% Example script to run 3D inversion and find optimum value of smoothing
%
% Eric Lindsey, July 2019

%% Create jointinv object

expNumber=102; 
scenario = Jointinv(expNumber);

% read the input file and set up matrices
scenario.run_setup();

% Run and plot inversion with input parameters
tic
scenario.run_inversion()
toc

%
ax=figure(1); clf
Mw=plot_jointinv_slipmodel(scenario,ax,1e5);
plot_jointinv_dataset(scenario,ax,1e5, 200);
title(['inversion and data fits, magnitude ' num2str(round(Mw,2))])
daspect([1 1 1])
view(37.5,37.5)

%% use ABIC to find the optimal smoothing weight

wgts = logspace(-4,1,20); % logarithmically spaced from 1e-4 to 1e0.
abic = zeros(1,length(wgts));

for i=1:length(wgts)
    scenario.userParams.smoothingWeights{1}{1}=wgts(i);
    % run the inversion
    scenario.run_inversion()
    % get ABIC value
    abic(i) = abic_smoothingonly(scenario);
end
bestwgt=wgts(abic==min(abic));

figure(2),clf
semilogx(wgts,abic-min(abic));
xlabel('smoothing weight')
ylabel('Delta ABIC, relative to minimum')
title(['ABIC vs. smoothing weight, minimum at ',num2str(bestwgt)])

% Run and plot inversion with optimized smoothing
scenario.userParams.smoothingWeights{1}{1}=bestwgt;
scenario.run_inversion()

ax=figure(3); clf
Mw=plot_jointinv_slipmodel(scenario,ax,1e5);
plot_jointinv_dataset(scenario,ax,1e5, 200);
title(['inversion and data fits with optimum smoothing, magnitude ' num2str(round(Mw,2))])
daspect([1 1 1])
view(37.5,37.5)