% script to find the optimum value of alpha and assess the range of acceptable models
% using beta
%
% Eric Lindsey, June 2019

% Part 1: optimize alpha

clear all, close all

basedir = '/Users/elindsey/Dropbox/projects/stress_constraint_backslip/Interseismic_StressInv/new_models';

% set parameters for a subduction zone
%model = '2d_models'
%model = 'japan'
model = 'cascadia'
%model = 'sumatra'

cd([basedir '/' model]);

run('../functions/set_jointinv_path.m');


%%
if strcmp(model, '2d_models')
    %guess for optimum alpha, maximum beta
    gA = -3.8; % for laplacian
    %gA = -1.9; % for stressKernel
    %guess for maximum beta:
    gBl = -4;
    gBh = -3.3;
    modeltime_long = .01;
    modeltime_short = .5;
elseif strcmp(model, 'japan')
    %guess for optimum alpha, maximum beta
    gA = -1; % for laplacian
    %gA = -1.9; % for stressKernel
    %guess for maximum beta:
    gBl = -4;
    gBh = -3.3;
    modeltime_long = 60;
    modeltime_short = 25;
elseif strcmp(model,'cascadia')
    %guess for optimum alpha, maximum beta
    gA = -0.15; % for laplacian
    %gA = -1; % for stressKernel
    %guess for maximum beta:
    gBl = -3.4;
    gBh = -2.6;
    modeltime_long = 60;
    modeltime_short = 25;
elseif strcmp(model,'sumatra')
% not implemented yet
     %guess for optimum alpha, maximum beta
     gA = -0.25; % for laplacian
     %gA = -1; % for stressKernel
     %guess for maximum beta:
     gBl = -3.2;
     gBh = -1.7;
     modeltime_long = 34;
     modeltime_short = 12;
end


Alist=[logspace(gA-7,gA-3,5),logspace(gA-2.5,gA-0.9,5),logspace(gA-0.7,gA-0.6,2),logspace(gA-0.5,gA+0.5,22),logspace(gA+0.6,gA+0.7,2),logspace(gA+0.9,gA+2.5,5),logspace(gA+2.8,gA+3,1)];
%Alist = logspace(gA-0.3,gA+0.3,7);
%
Blist = [1e-8,1e-7,1e-6,1e-5,10^-4.5,1e-4,logspace(gBl,gBh,20)];
%Blist = [1e-8,1e-4,logspace(gBl,gBh,5)];

figure(1),clf
plot(Alist, 1:length(Alist),'bx')
set(gca,'xscale','log')

est_time = now + length(Alist)*(modeltime_long/60/60/24);
disp(['Running ' num2str(length(Alist)) ' models. Estimated completion: ' datestr(est_time)]);

expNumber = 0;
%%
for expNumber = 0
    scenario=Jointinv(expNumber); 
    scenario.run_setup();
    
    sigma=.01;
    scenario.dataVector = scenario.dataVector + sigma*randn(size(scenario.dataVector));
%%
    bestmodels = zeros(length(Alist), length(scenario.modelVector));
    abics_alpha = zeros(length(Alist),1);
    chi2s_alpha = zeros(length(Alist),1);

    tic
    for i=1:length(Alist)
        %set alpha, fix beta=0
        scenario.userParams.smoothingWeights{1} = {Alist(i),0};
        %run inversion for this parameter set
        scenario.run_inversion();
        % compute abic and save model vector
        abics_alpha(i) = abic_smoothingonly(scenario);
        chi2s_alpha(i) = scenario.chi2;
        bestmodels(i,:)=scenario.modelVector;
    end
    toc % japan/exp0: ~60 sec/model. ~1hr for full set

    save(['./results/abic_optimize_alpha_exp' num2str(expNumber)], 'abics_alpha','chi2s_alpha','Alist','bestmodels','expNumber','model');

    %
    % find the range of acceptable alphas via ABIC, and find maximum beta allowed by F-test
    delta_abics = abics_alpha - min(abics_alpha);
    I_acceptable = find(delta_abics<10);

    Ibest= find(abics_alpha == min(abics_alpha),1);
    best_chi2 = chi2s_alpha(Ibest);
    bestA = Alist(Ibest);
    best_model_nobeta = bestmodels(Ibest,:)';
    
   %
    [A,B]=meshgrid(Alist(I_acceptable), Blist);
    P = 0*A;
    ftest_min = 0.05;

    est_time = now + length(P(:))*(modeltime_short/60/60/24);
    disp(['Running ' num2str(length(P(:))) ' models. Estimated completion: ' datestr(est_time)]);

    figure(3 + 5*expNumber)
    plot(A,B,'kx')
    xlim([1e-3,1e1])
    set(gca,'xscale','log','yscale','log')

    %%
    bestmodels_beta = zeros(length(Blist),length(I_acceptable),length(scenario.modelVector));

    tic
    for i=1:length(I_acceptable)
        alpha = Alist(I_acceptable(i));
        for j=1:length(Blist)
            scenario.userParams.smoothingWeights{1} = {A(j,i),B(j,i)};
            scenario.run_inversion();
            N = length(scenario.modelVector);
            P(j,i) = f_test_menke(best_chi2,N,scenario.chi2,N);
            bestmodels_beta(j,i,:) = scenario.modelVector;
%             if P(j,i) < ftest_min
%                 % stop running for this value of beta
%                 break
%             end
        end
    end
    toc % japan case / exp 0: 25sec/model. 
    
    save(['./results/ftest_maximize_beta_exp' num2str(expNumber)], 'A','B','P','Blist','bestmodels_beta','expNumber','model');

end

%%

for expNumber = [10,13]
    load(['./results/abic_optimize_alpha_exp' num2str(expNumber)]);
    %
    delta_abics = abics_alpha - min(abics_alpha);
    I_acceptable = find(delta_abics<10);

    %  plot results of alpha search

    styles={'-','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x','-x'};
    colors={'k','k','k','k','k','k','k','k','k','k','b','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k'};
    colors2={'k','k','k','k','k','k','k','k','k','k','c','k','r','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k'};
    figure(25)
    if expNumber ==10
        clf, hold on
    end

    plot(Alist,delta_abics,'-s','color',colors{expNumber+1},'linewidth',1),hold on
    %plot(Alist(I_acceptable),delta_abics(I_acceptable),'-','color','m','linewidth',2,'HandleVisibility','off')

    Ibest= find(abics_alpha == min(abics_alpha),1);
    best_chi2 = chi2s_alpha(Ibest);
    bestA = Alist(Ibest);
    best_model_nobeta = bestmodels(Ibest,:)';
    
    plot(Alist(I_acceptable),delta_abics(I_acceptable),'x','color','r','linewidth',2,'HandleVisibility','off')
    
    miny=0;
    maxy=1000;
    
    minx=min(Alist(I_acceptable));
    maxx=max(Alist(I_acceptable));
    
    patch([minx,minx,maxx,maxx,minx],[miny,maxy,maxy,miny,miny],colors{expNumber+1},'facealpha',0.1,'HandleVisibility','off','linestyle','none')
    
    %plot(log10(alphas),chi2red_alpha,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    xlabel('Smoothing parameter \alpha')
    ylabel('\Delta ABIC')
    set(gca,'fontsize',16)
    ax=gca;
    ax.YRuler.Exponent=0;
    %ylim([35,65])
    %xlim([1e-4,10])
    ylim([miny,maxy])
    set(gca,'xscale','log')
    grid on, box on

    if expNumber == 13
        title(['Model: ' model ])
        legend('No stress constraints','With stress constraints','location','southeast')
    end
end

set(gcf,'position',[800,500,750,400])
saveas2(['suppl_fig_alpha_' model '.pdf'],'-dpdf')


%% plot results of beta search

load(['./results/ftest_maximize_beta_exp' num2str(expNumber)]);

[A,B]=meshgrid(Alist(I_acceptable), Blist);

scenario=Jointinv(expNumber);
scenario.run_setup();

figure(4 + 5*expNumber ),clf


%scatter3(Amat(:),Bmat(:),Pmat(:),10000,Pmat(:),'.'),shading flat, hold on
contourf(A,B,P,100,'linestyle','none'),shading flat, hold on
plot(A(:),B(:),'kx')
colorbar
view(2)
contour(A, B, P,[0.05 0.05],'r','linewidth',2)
colormap(flipud(hot))
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-4,1e-2])
caxis([0 1])
% assess uncertainty

ploti=0;
makePlots=false;
if makePlots
    figure(10),clf
end

lA=log10(A);
lB=log10(B);

for i=1:length(I_acceptable)
    for j=1:length(Blist)
        if P(j,i) > 0.05
            if ploti==0
                max_model = squeeze(bestmodels_beta(j,i,:));
                min_model = max_model;
            end
            ploti = ploti +1;
            min_model=min(squeeze(bestmodels_beta(j,i,:)),min_model);
            max_model=max(squeeze(bestmodels_beta(j,i,:)),max_model);

            if makePlots
                subplot(9,12,ploti)
                scenario.sources{1}.geom.plotPatch(squeeze(bestmodels_beta(j,i,1:end/2))./scenario.sources{1}.Vpl')
                view(2)
                title([num2str(i),num2str(lA(i)) ', ' num2str(lB(i))])
            end
        end
    end
end

max_beta = max(B(P > 0.05));
log_max_beta = log10(max_beta)

modelUncertainty = (max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl;
% note that min/max are switched because the sign of slip is negative. Fix
% it here and use min/max consistently afterward.
max_coupling = -min_model(1:end/2)./scenario.sources{1}.Vpl;
min_coupling = -max_model(1:end/2)./scenario.sources{1}.Vpl;

figure(5 + 5*expNumber),clf
ax = subplot(1,3,3);
scenario.sources{1}.geom.plotPatch(modelUncertainty)
scenario.sources{1}.geom.plotPatch();
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
colorbar
view(2)
title('Difference for P > 5%')
caxis([0 1])
colormap(flipud(hot(10)))
daspect([1 1 1])

ax = subplot(1,3,1);
scenario.sources{1}.geom.plotPatch(max_coupling)
scenario.sources{1}.geom.plotPatch()
colorbar
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Max model')
daspect([1 1 1])

ax = subplot(1,3,2);
scenario.sources{1}.geom.plotPatch(min_coupling)
scenario.sources{1}.geom.plotPatch()
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Min model')
daspect([1 1 1])
colorbar

% save uncertainty and best model

% calc and plot best model
scenario.modelVector = best_model_nobeta;
%%

scenario.predVector = scenario.designMatrix * scenario.modelVector;

%%
results = calc_coupling_result_components(scenario);

figoffset=40;
plot_coupling_inversion(scenario,results, figoffset);

% save best model
save_coupling_inversion(scenario,results,'best_model_nobeta')




lat0=scenario.userParams.lat0;
lon0=scenario.userParams.lon0;

% save stresses
output_filename=['./results/rake_stress_points_' num2str(scenario.expNumber) '_' 'best_model_nobeta' '.gmt'];
save_geom_for_interp_GMT(output_filename,scenario.sources{1}.geom, results.rakeStress, lat0, lon0) 

% save uncertainty
output_filename=['./results/uncertainty_points_' num2str(scenario.expNumber) '.gmt'];
save_geom_for_interp_GMT(output_filename,scenario.sources{1}.geom, modelUncertainty, lat0, lon0) 


% save min model
output_filename=['./results/min_model_points_' num2str(scenario.expNumber) '.gmt'];
save_geom_for_interp_GMT(output_filename,scenario.sources{1}.geom, min_coupling, lat0, lon0) 


% save max model
output_filename=['./results/max_model_points_' num2str(scenario.expNumber) '.gmt'];
save_geom_for_interp_GMT(output_filename,scenario.sources{1}.geom, max_coupling, lat0, lon0) 



%% get full uncertainty from two models

if mod(expNumber,2) == 1
    true_min_coupling = load(['./results/min_model_points_' num2str(scenario.expNumber - 1) '.gmt']);
    full_uncertainty = max_coupling - true_min_coupling(:,3);
    
    figure(5 + 5*expNumber)
    ax=subplot(1,3,3)
    scenario.sources{1}.geom.plotPatch(full_uncertainty)
    scenario.sources{1}.geom.plotPatch();
    ax.Children(1).EdgeColor=[0.8 0.8 0.8];
    colorbar
    view(2)
    title('Max. difference for two models, P > 5%')
    caxis([0 1])
    colormap(flipud(hot(10)))
    daspect([1 1 1])
    
    % save uncertainty
    output_filename=['./results/full_uncertainty_points_' num2str(scenario.expNumber - 1) '_' num2str(scenario.expNumber) '.gmt'];
    save_geom_for_interp_GMT(output_filename,scenario.sources{1}.geom, full_uncertainty, lat0, lon0) 
    
end

