%% plot misfit vs beta for fixed alpha

clear all

fix_alpha=0.029;
fix_beta=1.5e-3;

colors={'r','b'};
styles={'s','x'};

% try a set of beta values
nbeta=200;
betas=logspace(-5.2,-0.5,nbeta);
alphas=[1e-4,logspace(-3.5,0.5,nbeta-1)];

figure(9),clf, hold on

%%

for expNumber = 0:1
    
    scenario = Jointinv(expNumber);
    scenario.run_setup();

%%
    chi2red_alpha=zeros(nbeta,1);
    chi2red_beta=zeros(nbeta,1);

    %fix the value of beta
    scenario.userParams.smoothingWeights{1}{2} = fix_beta;

    for ialpha=1:nbeta
        %set alpha
        scenario.userParams.smoothingWeights{1}{1} = alphas(ialpha);
        %run inversion for this parameter set
        scenario.run_inversion();
        % compute abic
        chi2red_alpha(ialpha) = chi2red_alpha(ialpha) + abic_alphabeta(scenario);

    end

    %fix the value of alpha
    scenario.userParams.smoothingWeights{1}{1} = fix_alpha;

    for ibeta=1:nbeta

        %set beta
        scenario.userParams.smoothingWeights{1}{2} = betas(ibeta);
        %run training inversion for this parameter set
        scenario.run_inversion();
        % compute abic
        chi2red_beta(ibeta) = chi2red_beta(ibeta) + abic_alphabeta(scenario);
    end

    %make test plots
    subplot(2,1,1)
    plot(log10(alphas),chi2red_alpha,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    xlabel('log_{10}(alpha)')
    ylabel('reduced \chi^2 on test dataset with fixed beta')
    set(gca,'fontsize',14)

    subplot(2,1,2)
    plot(log10(betas),chi2red_beta,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    xlabel('log_{10}(beta)')
    ylabel('reduced \chi^2 on test dataset with fixed alpha')
    set(gca,'fontsize',14)

    % save the results to a file
    alphafname=['test_data_2d/test3_alpha_fixbeta_exp' num2str(expNumber) '.dat'];
    alphatest=[alphas', chi2red_alpha];
    save(alphafname,'alphatest','-ASCII')
    %
    % save the results to a file
    betafname=['test_data_2d/test3_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    betatest=[betas', chi2red_beta];
    save(betafname,'betatest','-ASCII')
    
end

%% make final plots

figure(10),clf
styles={'-','--'};
colors={'r','b'};

for expNumber = 0:1
    
    % load the results
    alphafname=['test_data_2d/test3_alpha_fixbeta_exp' num2str(expNumber) '.dat'];
    alphatest = load(alphafname);
    alphas = alphatest(:,1);
    chi2red_alpha = alphatest(:,2);
    
    % load the results
    betafname=['test_data_2d/test3_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    betatest = load(betafname);
    betas = betatest(:,1);
    chi2red_beta = betatest(:,2);
    
    %make plots
    subplot(3,1,2)
    plot(alphas,chi2red_alpha,styles{expNumber+1},'color',colors{expNumber+1}),hold on
    %plot(log10(alphas),chi2red_alpha,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    xlabel('Smoothing parameter \alpha')
    ylabel('ABIC on test dataset')
    set(gca,'fontsize',14)
    %ylim([35,65])
    xlim([1e-4,1])
    set(gca,'xscale','log')
    grid on

    subplot(3,1,3)
    plot(betas,chi2red_beta,styles{expNumber+1},'color',colors{expNumber+1}),hold on
    %plot(log10(betas),chi2red_beta,'color',colors{expNumber+1},'marker',styles{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion', 'Best model','location','NorthWest')
    xlabel('Moment deficit penalty \beta')
    ylabel('ABIC on test dataset')
    set(gca,'fontsize',14)
    %ylim([35,65])
    xlim([9.999e-6,1e-1])
    set(gca,'xscale','log')
    grid on

end

subplot(3,1,1)

plot(scenario.datasets{1}.coordinates(:,1)/1e3, 1+scenario.dataVector/100,'kx'), hold on

%set(get(get(m,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% synthetic model
ibeta=94;
ialpha=98;

% set alpha
scenario.userParams.smoothingWeights{1}{1} = alphas(ialpha);
%set beta
scenario.userParams.smoothingWeights{1}{2} = betas(ibeta);
%run training inversion for this parameter set
scenario.run_inversion();

plot(scenario.datasets{1}.coordinates(:,1)/1e3, 1+scenario.predVector/100,'-r')


xlim([100,305])
ylim([0,1])
legend('Synthetic data','Best fit model','location','southeast')
%legend('Group 1','Group 2','Group 3','Group 4','Group 5','Model (excl. group 1)','location','southeast')
set(gca,'fontsize',14)
xlabel('Distance from trench (km)')
ylabel('GPS velocity / V_{pl}')

subplot(3,1,2)
plot(alphas(ialpha),chi2red_alpha(ialpha),'ro')
legend('Unconstrained inversion','Stress-constrained inversion','Best model','location','southwest')
    

subplot(3,1,3)
plot(betas(ibeta),chi2red_beta(ibeta),'ro')
legend('Unconstrained inversion','Stress-constrained inversion','Best model','location','northwest')

%%
% save figure
fig=gcf;
set(fig, 'Position',  [1500, 1900, 800, 700])
% doesn't work well with different types of axes
%stretch_fig_no_whitespace(fig);

%print(gcf,'suppl_fig_s1_cv_2d_example', '-dpdf')
saveas2('suppl_fig_s1_cv_2d_example.pdf')


