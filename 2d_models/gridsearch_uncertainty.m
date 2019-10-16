% search over alpha ane beta and plot the ABIC value at each point
clear all, close all

run('../functions/set_jointinv_path.m');

expNumber = 0;

scenario=Jointinv(expNumber);
scenario.run_setup();

%
betas=[1e-8,logspace(-7,-5.5,3),logspace(-4.9,-4,3),logspace(-3.75,-0.2,20),logspace(0.1,1,3)];
alphas=[1e-8,logspace(-7,-4.5,5),logspace(-4.1,-0.5,23),logspace(0.1,1,3)];

nalpha=length(alphas);
nbeta=length(betas);
figure(5),clf, hold on
plot(log10(alphas),1:length(alphas),'rs')
plot(log10(betas),(1:length(betas)),'bx')
grid on

[A,B] = meshgrid(alphas,betas);
[lA,lB] = meshgrid(log10(alphas),log10(betas));
C = (lA +2).^2 + (lB+3).^2;

figure(6),clf
pcolor(lA,lB,rand(size(C))),shading flat, hold on
colorbar

%

bestmodels = zeros(length(A(:)), length(scenario.modelVector));
abics = zeros(length(A(:)),1);
chi2s = abics;
sigma=1;
rng(2)
scenario.dataVector = scenario.dataVector + sigma*randn(size(scenario.dataVector));

tic
for i=1:length(A(:))
    %set alpha,beta
    scenario.userParams.smoothingWeights{1} = {A(i),B(i)};
    %run inversion for this parameter set
    scenario.run_inversion();
    % compute abic and save model vector
    abics(i) = abic_alphabeta(scenario);
    chi2s(i) = scenario.chi2;
    bestmodels(i,:)=scenario.modelVector;
end
toc
ABICs=reshape(abics,nbeta,nalpha);
CHI2s=reshape(chi2s,nbeta,nalpha);
%
save('results/gridsearch_ABICs', 'ABICs','lA','lB','bestmodels');

%% plot the results

deltaAbics = ABICs-min(ABICs(:));

figure(6),clf
contourf(lA,lB,deltaAbics,50,'LineStyle','none'),shading flat, hold on
contour(lA,lB,deltaAbics,[10 10],'r')
plot(lA(:),lB(:),'r.')
text(-1.7,-6,'\Delta ABIC < 10','color','r','fontsize',14,'rotation',90)
hcb=colorbar;
xlabel('log_{10}(\alpha)')
ylabel('log_{10}(\beta)')
%caxis([0 300])
title(hcb,'\Delta ABIC')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
daspect([ 1 1 1])
view(2)
ylim([-8 1])


%%
figure(7),clf 
minabic = min(abics);
imodel=0;
for i=1:nalpha*nbeta
  if abics(i) - minabic < 10
        plot(scenario.sources{1}.geom.xc(:,1),squeeze(bestmodels(i,:))),hold on
        imodel=imodel+1;
  end
end

imodel