% search over alpha ane beta and plot the ABIC value at each point
clear all, close all

run('../functions/set_jointinv_path.m');

expNumber = 0;

%%

% betas=[1e-8,logspace(-7,-5,4),logspace(-4.5,-3,5),logspace(-2.8,-1.2,10),logspace(-1,0,3)];
% alphas=[1e-8,logspace(-7,-5.5,3),logspace(-4.9,-2.2,12),logspace(-1.95,-1.2,7),logspace(-1.15,-0.85,6),logspace(-0.8,-0.5,4),logspace(-0.3,1,6)];

betas=[1e-8,logspace(-7,-5,4),logspace(-4.6,-3.9,4),logspace(-3.75,-2.5,11),logspace(-2.3,0,8)];
alphas=[1e-8,logspace(-7,-5.5,3),logspace(-4.9,-2.2,7),logspace(-1.95,-1.5,5),logspace(-1.4,-0.65,12),logspace(-0.55,0.4,6),logspace(0.7,1,2)];


nalpha=length(alphas);
nbeta=length(betas);
figure(5),clf, hold on
plot(log10(alphas),1:length(alphas),'rs')
plot(log10(betas),(1:length(betas)),'bx')
grid on

[A,B] = meshgrid(alphas,betas);
[lA,lB] = meshgrid(log10(alphas),log10(betas));

figure(6),clf
pcolor(lA,lB,rand(size(A))),shading flat, hold on
colorbar

%%

scenario=Jointinv(expNumber);
scenario.run_setup();

bestmodels = zeros(length(A(:)), length(scenario.modelVector));
%abics = zeros(length(A(:)),1);
chi2s = zeros(length(A(:)),1);

tic
for i=1:length(A(:))
    %set alpha,beta
    scenario.userParams.smoothingWeights{1} = {A(i),B(i)};
    %run inversion for this parameter set
    scenario.run_inversion();
    % compute abic and save model vector
    %abics(i) = abic_alphabeta(scenario);
    chi2s(i) = chi2_jointinv_singledataset(scenario);
    bestmodels(i,:)=scenario.modelVector;
end
toc
%ABICs=reshape(abics,nbeta,nalpha);
CHI2s=reshape(chi2s,nbeta,nalpha);

save(['results/gridsearch_CHI2s_exp' num2str(expNumber)], 'CHI2s','lA','lB','bestmodels');

%% plot the results

deltaAbics = ABICs-min(ABICs(:));

figure(6),clf
contourf(lA,lB,deltaAbics,5000,'LineStyle','none'),shading flat, hold on
contour(lA,lB,deltaAbics,[10 10],'r')
text(-1.4,-7,'\Delta ABIC < 10','color','r','fontsize',14,'rotation',90)
hcb=colorbar;
xlabel('log_{10}(\alpha)')
ylabel('log_{10}(\beta)')
caxis([0 5000])
title(hcb,'\Delta ABIC')
set(gca,'fontsize',14)
set(gca,'layer','top')
box on
grid on 
%daspect([ 1 1 1])
view(2)
%ylim([-8 0])
%xlim([-4,1])


%%

deltaAbics = ABICs-min(ABICs(:));

figure(6),clf
contourf(A,B,log10(deltaAbics),50,'LineStyle','none'),shading flat, hold on
%pcolor(A,B,log10(deltaAbics)),shading flat, hold on
plot(A(:),B(:),'k.','markersize',1) 

set(gca,'xscale','log','yscale','log')

contour(A,B,deltaAbics,[10 10],'b','linewidth',2)

%text(10^-0.3,10^-6,'\Delta ABIC < 10','color','b','fontsize',14,'rotation',90)

fig=gcf;
set(fig, 'Position',  [1500, 1900, 600, 500])
%stretch_fig_no_whitespace(fig,20);

xlabel('\alpha')
ylabel('\beta')

logColor = true;
if logColor == true
    c1=1;
    c2=20000;
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
    hcb = colorbar('Ticks',Ticks,'TickLabels',TickLabels);

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

%%
level=1;
h_axes = axes('position', hcb.Position, 'ylim', hcb.Limits, 'color', 'none', 'visible','off');
line(h_axes.XLim, level*[1 1], 'color', 'blue', 'parent', h_axes,'linewidth',2);


%% save figure

% save figure
fig=gcf;
set(fig, 'Position',  [1500, 1900, 600, 500])
% doesn't work well with different types of axes
stretch_fig_no_whitespace(fig,3);
%%
print(gcf,'figures/suppl_fig_abic_cascadia', '-dpdf')

%%

ploti=0;
%figure(1),clf
for i=1:length(lA(:))
    
    if deltaAbics(i)<100
        if ploti==0
            max_model = bestmodels(i,:);
            min_model = max_model;
        end
        ploti = ploti +1;
%         subplot(8,12,ploti)
%         scenario.sources{1}.geom.plotPatch(bestmodels(i,1:end/2)./scenario.sources{1}.Vpl')
%         view(2)
%         title([num2str(i),num2str(lA(i)) ', ' num2str(lB(i))])
        min_model=min(bestmodels(i,:),min_model);
        max_model=max(bestmodels(i,:),max_model);
    end
end

figure(2),clf
subplot(1,3,3)
%scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
scenario.sources{1}.geom.plotPatch((max_model(1:end/2) - min_model(1:end/2))./scenario.sources{1}.Vpl')
colorbar
view(2)
title('Difference for \Delta ABIC < 10')
caxis([0 1])
colormap(flipud(hot(10)))
subplot(1,3,1)
scenario.sources{1}.geom.plotPatch(-min_model(1:end/2)./scenario.sources{1}.Vpl')
colorbar
view(2)
title('Max model')
subplot(1,3,2)
scenario.sources{1}.geom.plotPatch(-max_model(1:end/2)./scenario.sources{1}.Vpl')
view(2)
title('Min model')


colorbar

