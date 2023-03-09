function Mw=plot_jointinv_slipmodel(scenario,ax,quiverScale)

    ax; hold on

    for i=1:length(scenario.sources)
        source = scenario.sources{i};
        if isa(source, 'Static_Halfspace_Fault_Source')
            if length(scenario.userParams.slipComponents) == 1 || (isfield(scenario.userParams, 'rakeOptions') && strcmp(scenario.userParams.rakeOptions, 'rakeFixed'))
                slipmag = abs(source.modelVector);
            else
                slipmag = sqrt(source.modelVector(1:end/2).^2 +source.modelVector(end/2+1:end).^2);
            end
            source.geom.plotPatch(slipmag); shading flat, hold on
            %source.geom.plotPatch()

            [~,Mw] = get_moment_and_magnitude(source.geom, slipmag);
            disp(['Fault ' num2str(i) ' magnitude: ' num2str(Mw)]);
             
            plot_jointinv_slipvectors(scenario,source,ax,quiverScale/5,'k')

        end
       
    end

    colormap(parula)
    colorbar
end