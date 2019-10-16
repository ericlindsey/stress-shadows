classdef Jointinv_Source < handle
    % SOURCE contains the general properties required by jointinv for any 'source'

    % It is inherited by:

    %   Static_Halfspace_Fault_Source

    
    properties
        
        fileName        % name of the file used to load the source
        modelVector     % M x 1  - model parameters, unraveled however you like

    end
    
    methods
        
        % Any 'Source' class inheriting this class must, at a minimum, implement the
        % following:
        %  Source(source_input, userParams) % class constructor
        %  calc_design_matrix(dataset, userParams)
        %  calc_smoothing_matrix(dataset, userParams)
        %  calc_constraint_matrix(dataset, userParams)
        
    end
    
end

