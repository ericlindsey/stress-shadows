% search over alpha ane beta and plot the ABIC value at each point
clear all, close all
run('../functions/set_jointinv_path.m');


%%
expNumber = 1;

% for Japan, with laplacian

%betas=[logspace(-8,-5.5,5),logspace(-5,-4.0,9),logspace(-3.9,-3.3,10),logspace(-3.1,-2.3,4),logspace(-2,0,4)];
%alphas=[1e-8,logspace(-7,-4.5,3),logspace(-4.1,-2.,5),logspace(-1.8,-1.3,5),logspace(-1.2,-0.7,9),logspace(-0.6,-0.3,3),logspace(-0.1,1,5)];

% for Cascadia, with laplacian

betas=[logspace(-8,-5.5,5),logspace(-5,-4.0,9),logspace(-3.9,-3.3,10),logspace(-3.2,-2.5,4),logspace(-2.1,0,4)];
alphas=[1e-8,logspace(-7,-4.7,3),logspace(-4.1,-0.9,8),logspace(-0.7,0.5,16),logspace(0.6,1.2,5),logspace(1.4,2,3)];


nalpha=length(alphas);
nbeta=length(betas);

[A,B] = meshgrid(alphas,betas);
[lA,lB] = meshgrid(log10(alphas),log10(betas));

figure(6),clf
plot(lA,lB,'k.'),shading flat, hold on
title(['Nalpha: ', num2str(nalpha),', Nbeta: ', num2str(nbeta), ', Nmodels: ', num2str(nalpha*nbeta)]);

%%

scenario=Jointinv(expNumber);
scenario.run_setup();

bestmodels = zeros(length(A(:)), length(scenario.modelVector));
abics = zeros(length(A(:)),1);
chi2s = zeros(length(A(:)),1);

tic
for i=1:length(A(:))
    %set alpha,beta
    scenario.userParams.smoothingWeights{1} = {A(i),B(i)};
    %run inversion for this parameter set
    scenario.run_inversion();
    % compute abic and save model vector
    abics(i) = abic_alphabeta_sum(scenario);
    chi2s(i) = chi2_jointinv_singledataset(scenario);
    bestmodels(i,:)=scenario.modelVector;
end
toc

%%
ABICs=reshape(abics,nbeta,nalpha);
CHI2s=reshape(chi2s,nbeta,nalpha);

save(['./results/gridsearch_laplacian_exp' num2str(expNumber)], 'ABICs','CHI2s','lA','lB','bestmodels');

%% plot the results

deltaAbics = ABICs-min(ABICs(find(ABICs(:)>0)));

figure(6),clf
contourf(lA,lB,deltaAbics,500,'LineStyle','none'),shading flat, hold on
contour(lA,lB,deltaAbics,[10 10],'r')
text(-1.9,-7,'\Delta ABIC < 10','color','r','fontsize',14,'rotation',90)
hcb=colorbar;
xlabel('log_{10}(\alpha)')
ylabel('log_{10}(\beta)')
caxis([0 20000])
title(hcb,'\Delta ABIC')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
daspect([ 1 1 1])
view(2)
%ylim([-8 0])

%%
% save figure
fig=gcf;
set(fig, 'Position',  [1500, 1900, 600, 500])
% doesn't work well with different types of axes
stretch_fig_no_whitespace(fig,2);


%%
deltaAbics = ABICs-min(ABICs(:));

figure(6),clf
contourf(A,B,log10(deltaAbics),50,'LineStyle','none'),shading flat, hold on
%pcolor(A,B,log10(deltaAbics)),shading flat, hold on
%plot(A(:),B(:),'k.','markersize',1) 


% [newA,newB] = meshgrid(newalphas,newbetas);
% plot(newA(:),newB(:),'k.','markersize',1) 

set(gca,'xscale','log','yscale','log')

contour(A,B,deltaAbics,[10 10],'k','linewidth',2)

%text(10^-0.3,10^-6,'\Delta ABIC < 10','color','b','fontsize',14,'rotation',90)

fig=gcf;
set(fig, 'Position',  [1500, 1900, 600, 500])
%stretch_fig_no_whitespace(fig,20);

xlabel('\alpha')
ylabel('\beta')

logColor = true;
if logColor == true
    c1=1;
    c2=10000;
    % set limits for the caxis
    caxis([log10(c1) log10(c2)]);
    % preallocate Ticks and TickLabels

    Ticks      = [0,1,2,3,4];
    num_of_ticks = length(Ticks);
    TickLabels = zeros(1,num_of_ticks);
    % distribute Ticks and TickLabels
    for n = 1:1:num_of_ticks

        %Ticks(n)      = log10(round(c2)/num_of_ticks*n);
        %TickLabels(n) = round(c2)/num_of_ticks*n;
        TickLabels(n)=10^Ticks(n);
    end
    % set Ticks and TickLabels
    hcb = colorbar('Ticks',Ticks,'TickLabels',TickLabels)
else
    hcb = colorbar
    caxis([0 300])
end


title(hcb,'\Delta ABIC')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
set(gca,'layer','top')

xlim([1e-5,1e1])
%daspect([ 1 1 1])
view(2)
colormap(flipud(hot))

%
level=1;
h_axes = axes('position', hcb.Position, 'ylim', hcb.Limits, 'color', 'none', 'visible','off');
line(h_axes.XLim, level*[1 1], 'color', 'black', 'parent', h_axes,'linewidth',2);


%%
print(gcf,'figures/suppl_fig_abic_japan', '-dpdf')


%% method 1: chi2 ratio test

CHI2_min_acceptable = min(CHI2s(find(deltaAbics<10)))
alpha_min_acceptable = min(A(find(deltaAbics<10)))

ratioCHI2s = (CHI2s - CHI2_min_acceptable)./CHI2_min_acceptable;

ploti=0;
makePlots=false;
if makePlots
    figure(1),clf
end
for i=1:length(lA(:))
    
    if A(i)> alpha_min_acceptable & ratioCHI2s(i) < 0.05
        if ploti==0
            max_model = bestmodels(i,:);
            min_model = max_model;
        end
        ploti = ploti +1;
        min_model=min(bestmodels(i,:),min_model);
        max_model=max(bestmodels(i,:),max_model);
        
        if makePlots
            subplot(9,12,ploti)
            scenario.sources{1}.geom.plotPatch(bestmodels(i,1:end/2)./scenario.sources{1}.Vpl')
            view(2)
            title([num2str(i),num2str(lA(i)) ', ' num2str(lB(i))])
        end
    end
end
%
max_beta = max(B(find(A> alpha_min_acceptable & ratioCHI2s < 0.05)));
log_max_beta = log10(max_beta)
%
figure(2),clf
ax = subplot(1,3,3);
%scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch();
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
colorbar
view(2)
title('Difference for \Delta \chi^2 < 5%')
caxis([0 1])
colormap(flipud(hot(10)))
daspect([1 1 1])

ax = subplot(1,3,1);
scenario.sources{1}.geom.plotPatch(-min_model(1:end/2)./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch()
colorbar
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Max model')
daspect([1 1 1])

ax = subplot(1,3,2);
scenario.sources{1}.geom.plotPatch(-max_model(1:end/2)./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch()
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Min model')
daspect([1 1 1])
colorbar
colorbar


%% Method 2: F test - from Menke p. 112

load('./results/gridsearch_laplacian_exp0.mat');

%%
CHI2_min_acceptable = CHI2s(find(ABICs == min(ABICs(find(lB==-8)))))
%CHI2_min_acceptable = min(CHI2s(find(deltaAbics<10)))
deltaAbics = ABICs-min(ABICs(:));

%CHI2_min_acceptable = min(CHI2s(find(deltaAbics<10)))
alpha_min_acceptable = min(A(find(deltaAbics<1)))



A=10.^lA;
B=10.^lB;

bestA = A(find(ABICs == min(ABICs(find(lB==-8)))))
bestB = B(find(ABICs == min(ABICs(find(lB==-8)))))

n = 2*scenario.sources{1}.geom.N;
EA = CHI2_min_acceptable * n;

P = 0*lA;

for i=1:length(lA(:))
    EB = CHI2s(i)*n;
    Fobs = (EA/n) / (EB/n);
    if( Fobs<1)
        Fobs = 1/Fobs;
    end
    P(i) = 1 - (fcdf(Fobs,n,n)-fcdf(1/Fobs,n,n));
end
figure(80),clf
contourf(A,B,P, 100, 'linestyle','none'),shading flat, hold on
contour(A,B,P,[0.05 0.05],'r','linewidth',2)
contour(A,B,deltaAbics,[10 10],'g','linewidth',2)
plot(bestA,bestB,'rx', 'markersize',10)
set(gca,'xscale','log','yscale','log')
hcb=colorbar;
xlabel('log_{10}(\alpha)')
ylabel('log_{10}(\beta)')
caxis([0 1])
title(hcb,'F-test P value')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
%daspect([ 1 1 1])
view(2)
ylim([1e-4,1e-3])
xlim([1e-2,1e-0])
level=0.05;
h_axes = axes('position', hcb.Position, 'ylim', hcb.Limits, 'color', 'none', 'visible','off');
line(h_axes.XLim, level*[1 1], 'color', 'r', 'parent', h_axes,'linewidth',2);
colormap(flipud(hot))

%%

ploti=0;
makePlots=false;
if makePlots
    figure(1),clf
end
for i=1:length(lA(:))
    
    if abs(lA(i) - log10(bestA))<0.2  & P(i) > 0.05
        if ploti==0
            max_model = bestmodels(i,:);
            min_model = max_model;
        end
        ploti = ploti +1;
        min_model=min(bestmodels(i,:),min_model);
        max_model=max(bestmodels(i,:),max_model);
        
        if makePlots
            subplot(9,12,ploti)
            scenario.sources{1}.geom.plotPatch(bestmodels(i,1:end/2)./scenario.sources{1}.Vpl')
            view(2)
            title([num2str(i),num2str(lA(i)) ', ' num2str(lB(i))])
        end
    end
end
%
max_beta = max(B(find(A> alpha_min_acceptable & P < 0.05)));
log_max_beta = log10(max_beta)
%
figure(2),clf
ax = subplot(1,3,3);
%scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch();
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
colorbar
view(2)
title('Difference for P > 5%')
caxis([0 1])
colormap(flipud(hot(10)))
daspect([1 1 1])

ax = subplot(1,3,1);
scenario.sources{1}.geom.plotPatch(-min_model(1:end/2)./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch()
colorbar
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Max model')
daspect([1 1 1])

ax = subplot(1,3,2);
scenario.sources{1}.geom.plotPatch(-max_model(1:end/2)./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch()
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Min model')
daspect([1 1 1])
colorbar
colorbar






%% method 3: pure ABIC statistics

ploti=0;
makePlots=false;
if makePlots
    figure(1),clf
end
for i=1:length(lA(:))
    
    if deltaAbics(i) < 10
        if ploti==0
            max_model = bestmodels(i,:);
            min_model = max_model;
        end
        ploti = ploti +1;
        min_model=min(bestmodels(i,:),min_model);
        max_model=max(bestmodels(i,:),max_model);
        
        if makePlots
            subplot(9,12,ploti)
            scenario.sources{1}.geom.plotPatch(bestmodels(i,1:end/2)./scenario.sources{1}.Vpl')
            view(2)
            title([num2str(i),num2str(lA(i)) ', ' num2str(lB(i))])
        end
    end
end
%
figure(2),clf
ax = subplot(1,3,3);
%scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch();
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
colorbar
view(2)
title('Difference for \Delta ABIC < 10')
caxis([0 1])
colormap(flipud(hot(10)))
daspect([1 1 1])

ax = subplot(1,3,1);
scenario.sources{1}.geom.plotPatch(-min_model(1:end/2)./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch()
colorbar
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Max model')
daspect([1 1 1])

ax = subplot(1,3,2);
scenario.sources{1}.geom.plotPatch(-max_model(1:end/2)./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch()
ax.Children(1).EdgeColor=[0.8 0.8 0.8];
view(2)
title('Min model')
daspect([1 1 1])
colorbar
