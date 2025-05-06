function Mw=plot_jointinv_slipmodel(scenario,ax,quiverScale)

    ax; hold on

    for i=1:length(scenario.sources)
        source = scenario.sources{i};
        N=length(source.Vpl);
        if isa(source, 'Static_Halfspace_Fault_Source')
            if length(scenario.userParams.slipComponents) == 1
                slip = source.modelVector;
                source.geom.plotPatch(slip); shading flat, hold on
            elseif  (isfield(scenario.userParams, 'faultOptions') && strcmp(scenario.userParams.faultOptions, 'rakeFixed'))
                slip = source.modelVector(1:N);
                source.geom.plotPatch(slip); shading flat, hold on
            else
                slip1 = source.modelVector(1:N/2);
                slip2 = source.modelVector(N/2+1:N);
            
                subplot(2,1,1)
                source.geom.plotPatch(slip1); shading flat, hold on
                subplot(2,1,2)
                source.geom.plotPatch(slip2); shading flat, hold on

            end
        end
       
    end

    colormap(parula)
    colorbar
end