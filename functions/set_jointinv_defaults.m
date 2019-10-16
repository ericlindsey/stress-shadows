function userParams = set_jointinv_defaults(userParams)
    % set defaults for various Jointinv modules

    % for Static_GPS_Dataset
    if ~isfield(userParams, 'coordType')
        userParams.coordType = 'geographic';
    end
    if ~isfield(userParams, 'dataComponents')
        userParams.dataComponents = [1,2,3];
    end
    if ~isfield(userParams, 'minGPSError')
        userParams.minGPSError = 0.5;
    end

    % for Static_Halfspace_Fault_Source
    if ~isfield(userParams, 'shearModulus')
        userParams.shearModulus = 30e3;
    end
    if ~isfield(userParams, 'poissonsRatio')
        userParams.poissonsRatio = 0.25;
    end
    if ~isfield(userParams, 'slipComponents')
        userParams.slipComponents = [1,2];
    end

    % for inversion
    if ~isfield(userParams, 'smoothingType')
        % all smoothing set to 'None'
        [userParams.smoothingType{1:length(userParams.sourceTypes)}]=deal('None');
    end
    if ~isfield(userParams, 'constraints')
        userParams.constraints={};
    end
    if ~isfield(userParams, 'stressKernelFolder')
        userParams.stressKernelFolder = '.';
    end
    if ~isfield(userParams, 'frictionCoeff')
        userParams.frictionCoeff = 0.6;
    end
    
    
end
