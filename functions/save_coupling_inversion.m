function save_coupling_inversion(scenario,results,descriptor)
    % save result structures and some basic outputs to be used for GMT
    % 'descriptor' is a string to be added to the filenames to distinguish
    % multiple runs
    % Eric Lindsey, June 2019

    if ~isdir('./results')
        mkdir('results');
    end
    
    % save result objects? Not used currently, large file
    % save(['results/results_' num2str(scenario.expNumber) '.mat']);
    
    %get lat/lon for projection center from userParams
    lat0=scenario.userParams.lat0;
    lon0=scenario.userParams.lon0;
    
    % save coupling values in two different formats: polygons and points
    output_filename=['./results/coupling_polygons_' num2str(scenario.expNumber) '_' descriptor '.gmt'];
    save_geom_for_GMT(output_filename,scenario.sources{1}.geom, results.coupling, lat0, lon0) 
    output_filename=['./results/coupling_points_' num2str(scenario.expNumber) '_' descriptor '.gmt'];
    save_geom_for_interp_GMT(output_filename,scenario.sources{1}.geom, results.coupling, lat0, lon0) 

    % save predicted GPS vectors
    [gpslat,gpslon] = xy_to_latlon_polyconic(scenario.datasets{1}.coordinates(:,1)/1e3, scenario.datasets{1}.coordinates(:,2)/1e3, lon0, lat0);
    model_filename=['./results/predvectors_' num2str(scenario.expNumber) '_' descriptor '.gmt'];
    model_out=[gpslon gpslat -scenario.predVector(1:2:end) -scenario.predVector(2:2:end) [0 0 0].*gpslon ]; %#ok<NASGU>
    save(model_filename,'model_out','-ASCII')
    
    disp(['Wrote output for experiment ' num2str(scenario.expNumber) '.'])
    
end

