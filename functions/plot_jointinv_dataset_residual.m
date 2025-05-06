function plot_jointinv_dataset_residual(scenario, ax, vecScale, ptScale)

    ax, hold on

    for i=1:length(scenario.datasets)
        dataset = scenario.datasets{i};
        if isa(dataset,'Static_GPS_Dataset')
            numComponents = length(scenario.userParams.dataComponents);
            if numComponents == 3
                scatter(dataset.coordinates(:,1),dataset.coordinates(:,2),ptScale,dataset.dataVector(3:3:end)-dataset.predVector(3:3:end),'filled','MarkerEdgeColor',[0 0 0 ])
            end
            scaled_quiver(dataset.coordinates(:,1),dataset.coordinates(:,2),dataset.dataVector(1:numComponents:end)-dataset.predVector(1:numComponents:end), ...
                dataset.dataVector(2:numComponents:end)-dataset.predVector(2:numComponents:end),vecScale,{'k','linewidth',2})
            elseif isa(dataset,'Static_LOS_Dataset')
            % LOS datasets are single-component, so the plotting is much simpler
            scatter(dataset.coordinates(:,1),dataset.coordinates(:,2),ptScale,dataset.dataVector-dataset.predVector,'filled')
        else
            error('Dataset type not recognized for plotting');
        end
    end
    colormap(bluewhitered)
    colorbar
end