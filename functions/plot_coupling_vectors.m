function plot_coupling_vectors(scenario,vecScale)
    % simple script to plot the data and model vectors on the same scale
    %
    % Eric Lindsey, June 2019
    xgps=scenario.datasets{1}.coordinates(:,1);
    ygps=scenario.datasets{1}.coordinates(:,2);
    
    gps_E=scenario.dataVector(1:2:end);
    gps_N=scenario.dataVector(2:2:end);
    
    scaled_quiver(xgps,ygps,gps_E,gps_N,vecScale, {'color','k'})
    
    if ~isempty(scenario.predVector)
        model_E=scenario.predVector(1:2:end);
        model_N=scenario.predVector(2:2:end);
        scaled_quiver(xgps,ygps,model_E,model_N,vecScale, {'color','r'})
    end
    
end