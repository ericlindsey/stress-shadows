% plot misfit vs beta for variable alpha

colors={'r','b'};

% try a set of beta values
nbeta=51;
betas=[1e-9,1e-7,1e-6,1e-5,2e-5,4e-5,6e-5,8e-5,logspace(-3,-2.5,nbeta-8)];

% nbeta=1;
% betas=1e-4;

no_noise_file = load('test_data_2d/zero_noise_data.dat');
zero_noise_data = no_noise_file(:,4);

figure(8),clf

for expNumber=0:0
    
    scenario = Jointinv(expNumber);
    scenario.run_setup();
    
    % add data noise
    sigma=1;
    scenario.dataVector = zero_noise_data + sigma*randn(size(scenario.dataVector));

    %fix the value of alpha
    fix_alpha=1e-2;
    scenario.userParams.smoothingWeights{1}{1} = fix_alpha;

    abics=zeros(nbeta,1);
    best_alphas=zeros(nbeta,1);
    
    for ibeta=1:nbeta
        %set beta
        scenario.userParams.smoothingWeights{1}{2} = betas(ibeta);
        %k-fold CV to find best alpha
        best_alpha = find_best_alpha_abic(scenario);
        scenario.userParams.smoothingWeights{1}{1} = 10^best_alpha;
        %run inversion for this parameter set
        scenario.run_inversion();
        %save results
        abics(ibeta) = abic_alphabeta(scenario);
        best_alphas(ibeta)=best_alpha;
    end
    %
    % save the results to a file
%     betafname=['test_data_2d/test_abic_beta_exp' num2str(expNumber) '.dat'];
%     betatest=[betas', best_alphas, chi2red];
%     save(betafname,'betatest','-ASCII')
 
    subplot(2,1,1)
    plot(log10(betas),best_alphas,'-x','color',colors{expNumber+1}), hold on
    xlabel('log_{10}(beta)')
    ylabel('ABIC-determined log_{10}(alpha)')
    %legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    line([min(log10(betas)),max(log10(betas))],[mean(best_alphas),mean(best_alphas)],'color',colors{expNumber+1})
    set(gca,'fontsize',14)
    ylim([-2,-1])
    
    subplot(2,1,2)
    plot(log10(betas),abics,'x','color',colors{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    xlabel('log_{10}(beta)')
    ylabel('reduced \chi^2 with ABIC-determined alpha')
    set(gca,'fontsize',14)
    
end


%% plot ABIC vs beta instead, for fixed alpha

% global_best_alpha = -1.85;
fix_alpha=0.029;

colors={'r','b'};

% try a set of beta values
nbeta=200;
betas=[1e-8,1e-5,logspace(-4,-2.5,nbeta-2)];

figure(9),clf

for expNumber=0:1
    
    scenario = Jointinv(expNumber);
    scenario.run_setup();

    %fix the value of alpha
    scenario.userParams.smoothingWeights{1}{1} = fix_alpha;

    chi2red=nan*zeros(nbeta,1);
    best_alphas=zeros(nbeta,1);
    
    for ibeta=1:nbeta
        %set beta
        scenario.userParams.smoothingWeights{1}{2} = betas(ibeta);
        %run inversion for this parameter set
        scenario.run_inversion();
        %save results
        %chi2red(ibeta) = (scenario.predVector - scenario.dataVector)'*scenario.datasets{1}.covarianceMatrix*(scenario.predVector - scenario.dataVector)/length(scenario.dataVector);
        chi2red(ibeta) = abic_alphabeta(scenario);
    end
    %
    % save the results to a file
    betafname=['test_data_2d/test_abic_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    betatest=[betas', chi2red];
    save(betafname,'betatest','-ASCII')

    plot(log10(betas),chi2red,'color',colors{expNumber+1}),hold on
    legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
    xlabel('log_{10}(beta)')
    ylabel('ABIC with fixed alpha')
    set(gca,'fontsize',14)
    
       ylim([50,250])
end


%% make plot

figure(8),clf

for expNumber = 0:1
    % load results for CV-determined alpha
    betafname=['test_data_2d/test_abic_beta_exp' num2str(expNumber) '.dat'];
    betatest=load(betafname);
    betas=betatest(:,1);
    best_alphas=betatest(:,2);
    chi2red=betatest(:,3);
    
    mean(best_alphas)
    subplot(2,1,1)
    plot(betas,10.^best_alphas,'-','color',colors{expNumber+1}), hold on
     
    % load results for fixed alpha
    betafname=['test_data_2d/test_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    betatest=load(betafname);
    betas=betatest(:,1);
    chi2red=betatest(:,2);

    subplot(2,1,2)
    plot(betas,chi2red,'-','color',colors{expNumber+1}),hold on
    

end

% edit figure styles

subplot(2,1,1)
l=line([1e-9,10^-2.5],[0.029,0.029],'color','k','linewidth',1);
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(0.35e-2,0.029,'0.029','fontsize',14)
xlabel('Moment deficit penalty \beta')
ylabel('CV-determined \alpha')
legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
set(gca,'fontsize',14,'xscale','log')
bestalphalist(expNumber+1)=mean(best_alphas);
ylim([0,.06])
grid on

subplot(2,1,2)
legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
xlabel('Moment deficit penalty \beta')
ylabel('red. \chi^2 with \alpha = 0.029')
set(gca,'fontsize',14,'xscale','log')
grid on
ylim([0.9,1.1])

%%
fig=gcf;
% set figure/axes size with no whitespace padding
set(fig, 'Position',  [1500, 1900, 700, 500])
%
stretch_fig_no_whitespace(fig);

print(gcf,'supplement_cv_beta_search', '-dpdf')
