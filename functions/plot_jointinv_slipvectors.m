function plot_jointinv_slipvectors(scenario,source,ax,quiverScale, vecColor)

    ax;
    hold on

    % get the model vector strike/dip components
    if isfield(scenario.userParams, 'faultOptions')
        if strcmp(scenario.userParams.faultOptions, 'rakeCoordinates')
            m = source.Rmat * source.modelVector;
        elseif strcmp(scenario.userParams.faultOptions, 'rakeFixed')
            m = source.Rmat * source.modelVector;
        end
    else
        if scenario.userParams.slipComponents == 1
            m = [source.modelVector; 0*source.modelVector];
        elseif scenario.userParams.slipComponents == 2
            m = [0*source.modelVector; source.modelVector];
        else
            m = source.modelVector;
        end
    end

    % convert to slip rate magnitude and coupling
    strikeSlip = m(1:end/2);
    dipSlip = m(end/2+1:end);

    source.geom.plotSlipVectors(strikeSlip,dipSlip,quiverScale,vecColor)

    
end