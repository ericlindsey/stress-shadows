function plot_jointinv_dataset_3panel_residual(dataset, ax, ptScale, redchi2_i, npts_i,datasetname,fltx,flty)

    ax; hold on
    
    chi2 = ((dataset.predVector - dataset.dataVector)'*inv(dataset.covarianceMatrix)*(dataset.predVector - dataset.dataVector))/length(dataset.dataVector);
    
    rms = (dataset.predVector - dataset.dataVector)'*(dataset.predVector - dataset.dataVector);
    
    subplot(1,3,1)
    % LOS datasets are single-component, so the plotting is much simpler
    scatter(dataset.coordinates(:,1)/1e3,dataset.coordinates(:,2)/1e3,ptScale,dataset.dataVector,'filled')
    cdata=ceil(max(abs(dataset.dataVector)));
    clim([-cdata,cdata])
    colormap(bluewhitered)
    colorbar
    title(datasetname)
    hold on, grid on, box on
    set(gca,'color',[0.9,0.9,0.9])
    plot(fltx,flty,'-k')
    xlabel('Distance East (km)')
    ylabel('Distance North (km)')
    xlim([-100,100])

    subplot(1,3,2)
    scatter(dataset.coordinates(:,1)/1e3,dataset.coordinates(:,2)/1e3,ptScale,dataset.predVector,'filled')
    clim([-cdata,cdata])
    colormap(bluewhitered)
    colorbar
    hold on, grid on, box on
    set(gca,'color',[0.9,0.9,0.9])
    plot(fltx,flty,'-k')
    xlabel('Distance East (km)')
    ylabel('Distance North (km)')
    xlim([-100,100])
    title(['Npts = ',num2str(npts_i)])
    
    subplot(1,3,3)
    scatter(dataset.coordinates(:,1)/1e3,dataset.coordinates(:,2)/1e3,ptScale,dataset.dataVector-dataset.predVector,'filled')
    %cresid=ceil(max(2*abs(dataset.dataVector-dataset.predVector)))/2;
    clim([-cdata,cdata])
    colormap(bluewhitered)
    colorbar
    hold on, grid on, box on
    set(gca,'color',[0.9,0.9,0.9])
    plot(fltx,flty,'-k')
    title(['red. \chi^2 ',num2str(round(redchi2_i,2))])
    xlabel('Distance East (km)')
    ylabel('Distance North (km)')
    xlim([-100,100])

end