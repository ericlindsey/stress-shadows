%% ABIC determination of best-fitting hyperparameters for Sumatra

% first find the most likely best-fitting pair

scenario = Jointinv(0);
scenario.run_setup();

%[optvals] = find_best_alphabeta_abic(scenario)

%% create a variable-resolution slice along alpha and beta to make the figure

fix_alpha=10^-1;
fix_beta=1e-8;

colors={'r','b'};
styles={'s','x'};

%betas=[1e-8,1e-6,1e-5,10^-4.5,10^-4.25,logspace(-4,-2.5,25),logspace(-2.25,0,10)];
%alphas=[1e-8,1e-5,1e-4,logspace(-3.75,-3.35,3),logspace(-3.25,-0.25,39),logspace(0,1,5)];
betas=[logspace(-4.5,-2.5,15)];%,logspace(-6,-4,5),logspace(-3.9,-1.2,24),logspace(-1,1,10)];
alphas=[logspace(-1,-0.5,1)];

nalpha=length(alphas);
nbeta=length(betas);
figure(5),clf, hold on
plot(log10(alphas),1:length(alphas),'rs')
plot(log10(betas),(1:length(betas)),'bx')
%ylim([-10,20])
grid on
expNumber=0;
%%
for expNumber = 0:0
    
    scenario = Jointinv(expNumber);
    scenario.run_setup();

    abics_alpha=zeros(nalpha,1);
    abics_beta=zeros(nbeta,1);

    %fix the value of beta
    scenario.userParams.smoothingWeights{1}{2} = fix_beta;
    
%     % add data noise
%     sigma=1;
%     scenario.dataVector = scenario.dataVector + sigma*randn(size(scenario.dataVector));

    for ialpha=1:nalpha
        %set alpha
        scenario.userParams.smoothingWeights{1}{1} = alphas(ialpha);
        %run inversion for this parameter set
        scenario.run_inversion();
        % compute abic
        abics_alpha(ialpha) = abic_alphabeta_sum(scenario);

    end
%%
    %fix the value of alpha
    scenario.userParams.smoothingWeights{1}{1} = fix_alpha;

    for ibeta=1:nbeta

        %set beta
        scenario.userParams.smoothingWeights{1}{2} = betas(ibeta);
        %run training inversion for this parameter set
        scenario.run_inversion();
        % compute abic
        abics_beta(ibeta) = abic_alphabeta_sum(scenario);
    end

    %%
    figure(50)
    % make test plots
    subplot(2,1,1)
    plot(log10(alphas),abics_alpha,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    xlabel('log_{10}(alpha)')
    ylabel('ABIC with fixed beta')
    set(gca,'fontsize',14)
%%
    subplot(2,1,2)
    plot(log10(betas),abics_beta,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    xlabel('log_{10}(beta)')
    ylabel('ABIC with fixed alpha')
    set(gca,'fontsize',14)
%%
    % save the results to a file
    alphafname=['./results/test_alpha_fixbeta_exp' num2str(expNumber) '.dat'];
    alphatest=[alphas', abics_alpha];
    save(alphafname,'alphatest','-ASCII')
    %
    % save the results to a file
    betafname=['./results/test_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    betatest=[betas', abics_beta];
    save(betafname,'betatest','-ASCII')
    
end

%% make final plots

figure(10),clf
styles={'-','--'};
colors={'r','b'};

for expNumber = 0:1
    
    % load the results
    alphafname=['./results/test_alpha_fixbeta_exp' num2str(expNumber) '.dat'];
    alphatest = load(alphafname);
    alphas = alphatest(:,1);
    abics_alpha = alphatest(:,2);
    
    % load the results
    betafname=['./results/test_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    betatest = load(betafname);
    betas = betatest(:,1);
    abics_beta = betatest(:,2);
    
    %make plots
    subplot(2,1,1)
    [~,Isort]=sort(alphas);
    plot(alphas(Isort),abics_alpha(Isort),styles{expNumber+1},'color',colors{expNumber+1}),hold on
    %plot(log10(alphas),chi2red_alpha,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    xlabel('Smoothing parameter \alpha')
    ylabel('ABIC on test dataset')
    set(gca,'fontsize',14)
    %ylim([35,65])
    xlim([1e-4,1])
    set(gca,'xscale','log')
    grid on

    subplot(2,1,2)
    
    [~,Isort]=sort(betas);
    plot(betas(Isort),abics_beta(Isort),styles{expNumber+1},'color',colors{expNumber+1}),hold on
    %plot(log10(betas),chi2red_beta,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion', 'Best model','location','NorthWest')
    xlabel('Moment deficit penalty \beta')
    ylabel('ABIC on test dataset')
    set(gca,'fontsize',14)
    %ylim([35,65])
    xlim([9.999e-5,1])
    set(gca,'xscale','log')
    grid on

end

%select best points
ibeta=17;
ialpha=15;

subplot(2,1,1)
plot(alphas(ialpha),abics_alpha(ialpha),'ks','markersize',10)
%legend('Unconstrained inversion','Stress-constrained inversion','Best model','location','southwest')
ylim([11000 13000]) 

subplot(2,1,2)
plot(betas(ibeta),abics_beta(ibeta),'ks','markersize',10)
legend('Unconstrained inversion','Stress-constrained inversion','location','northwest')
ylim([11000 13000])

%%
% save figure
fig=gcf;
set(fig, 'Position',  [1500, 1900, 800, 500])
% doesn't work well with different types of axes
stretch_fig_no_whitespace(fig);
%%
print(gcf,'figures/suppl_fig_s9_cv_cascadia', '-dpdf')


