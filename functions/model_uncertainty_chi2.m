% plot chi-squared grid search results for a given model and experiment
clear all, close all
%run('../functions/set_jointinv_path.m');
%%
basedir = '/Users/elindsey/Dropbox/projects/stress_constraint_backslip/Interseismic_StressInv/new_models';
modeldir = 'japan';
expNumber = 0;

cd([basedir '/' modeldir]);

scenario = Jointinv(expNumber);
scenario.run_setup();
scenario.run_inversion();

%%
load(['./results/gridsearch_CHI2s_exp' num2str(expNumber) '.mat']);

ratioCHI2s = ( CHI2s ) ./ min(CHI2s(:));

%%
figure(6),clf
contourf(lA,lB,ratioCHI2s,5000,'LineStyle','none'),shading flat, hold on
%contour(lA,lB,deltaCHI2s,[10 10],'r')
%text(-1.9,-7,'\Delta \chi^2 < 10','color','r','fontsize',14,'rotation',90)
hcb=colorbar;
xlabel('log_{10}(\alpha)')
ylabel('log_{10}(\beta)')
%%
caxis([1 1.11])
%
title(hcb,'\Delta \chi^2 / min \chi^2')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
daspect([ 1 1 1])
view(2)
%ylim([-8 0])

% save figure
fig=gcf;
%set(fig, 'Position',  [1500, 1900, 600, 500])
% doesn't work well with different types of axes
%stretch_fig_no_whitespace(fig,2);


%%
% japan case
%Imodels = find(lA>=-2.3 & lA<= -1.8 & lB>=-4 & lB<=-2);

%cascadia case
Imodels = find(lA>=-1.15 & lA<= -0.85 & lB>=-4 & lB<=-2);

newABICs=0*lA;
%%
for i=1:length(Imodels)
    % set up model with known solution
    scenario.modelVector = bestmodels(Imodels(i),:)';
    scenario.userParams.smoothingWeights{1} = {10^lA(Imodels(i)),10^lB(Imodels(i))};

    %abics_old(i)=abic_alphabeta(scenario);
    newABICs(Imodels(i))=abic_alphabeta_sum(scenario);

end
%%
figure(7),clf
plot(logbetas,abics_old,'b'),hold on
plot(logbetas,abics_sum-min(abics_sum),'--r')
plot(logbetas,ratio_chi2s,'k')
legend('old','sum')




%%

DnewABICs = newABICs - min(newABICs(newABICs>0));
figure(8),clf
%contourf(lA,lB,DnewABICs,100,'LineStyle','none'),shading flat, hold on
pcolor(lA,lB,DnewABICs),shading flat, hold on
plot(lA,lB,'rx')

%contour(lA,lB,deltaCHI2s,[10 10],'r')
%text(-1.9,-7,'\Delta \chi^2 < 10','color','r','fontsize',14,'rotation',90)
hcb=colorbar;
xlabel('log_{10}(\alpha)')
ylabel('log_{10}(\beta)')
%
caxis([0 10])
%
title(hcb,'\Delta ABIC')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
daspect([ 1 1 1])
view(2) 


%% manual test
alpha=10^-2.1;

scenario.modelVector = bestmodels(Imodels(i),:)';
scenario.userParams.smoothingWeights{1} = {alpha,10^-8};
abic=abic_alphabeta_sum(scenario);

figure(9)
plot(log10(alpha),abic,'bx'), hold on
