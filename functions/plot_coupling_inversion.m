function plot_coupling_inversion(scenario, results, figoffset)
% function to plot:
% figure(1): coupling and data fit, along with slip vectors on the non-coupled patches
% figure(2): slip rate components
% figure(3): stress rate components
%
% Eric Lindsey, June 2019

disp(['chi2: ' num2str(results.chi2)])

% plot coupling, slip vectors, and data/model vectors
figure(1 + figoffset),clf

scenario.sources{1}.geom.plotPatch(results.coupling), hold on
view([0,90])
daspect([1 1 1])
colorbar
vecScale=2e3;
plot_coupling_vectors(scenario,vecScale)
scenario.sources{1}.geom.plotSlipVectors(results.Vpl_ss + results.strikeSlip, results.Vpl_ds + results.dipSlip,100,'k')
title('Coupling')
% caxis([-1 1])
%colormap(gca,parula(10))
cm=[255 255 255;        
242.5 244.5 132;
246 214 68.5;
251.5 181 25.002;
250.5 149 6.5009;
243 119 26;
229.5 92.501 47.499;
210.5 70.499 67.501;
187 54.5 84.5;
161 42.501 98.499;
134 33 107;
107 23.5 110]/255;
colormap(gca,cm)


% plot slip components in a second figure
figure(2 + figoffset),clf

subplot(2,2,1)
scenario.sources{1}.geom.plotPatch(-results.strikeSlip)
view([0,90])
daspect([1 1 1])
colorbar
title('Strike Slip')

subplot(2,2,2)
scenario.sources{1}.geom.plotPatch(-results.dipSlip)
view([0,90])
daspect([1 1 1])
colorbar
title('Dip Slip')

subplot(2,2,3)
scenario.sources{1}.geom.plotPatch(results.rakeSlip)
view([0,90])
daspect([1 1 1])
colorbar
title('Rake Slip')

subplot(2,2,4)
scenario.sources{1}.geom.plotPatch(results.rakePerpSlip)
view([0,90])
daspect([1 1 1])
colorbar
title('Rake-perp Slip')


%plot stress components in a third figure
figure(3 + figoffset),clf

subplot(2,2,1)
scenario.sources{1}.geom.plotPatch(results.strikeStress)
view([0,90])
daspect([1 1 1])
colorbar
title('Strike stress')
colormap(gca,bluewhitered)

subplot(2,2,2)
scenario.sources{1}.geom.plotPatch(results.dipStress)
view([0,90])
daspect([1 1 1])
colorbar
title('Dip stress')
colormap(gca,bluewhitered)

subplot(2,2,3)
scenario.sources{1}.geom.plotPatch(results.rakeStress)
view([0,90])
daspect([1 1 1])
colorbar
title('Rake stress')
colormap(gca,bluewhitered)

subplot(2,2,4)
scenario.sources{1}.geom.plotPatch(results.rakePerpStress)
view([0,90])
daspect([1 1 1])
colorbar
title('Rake-perp stress')
colormap(gca,bluewhitered)
end
