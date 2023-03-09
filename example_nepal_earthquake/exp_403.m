function params = exp_403()

    % This parameter file configures jointinv for:
    % GPS-only model for Nepal

    % note: defaults are set in the script defaults_Jointinv.m

    %%% Dataset information %%%
    params.datasetTypes     = {'Static_LOS_Dataset'}; % cell list of dataset types (matching the corresponding Object filename)
    params.datasetFilenames = {'data/varres_T048_insar.txt'}; % cell list of data files to read, 1-to-1 correspondence to the above list
    % Static_LOS_Dataset parameters
    % note, these are approximate only - the fault geometry was projected using a different system
    params.LOSlat0=28.155;
    params.LOSlon0=84.72;
    params.LOSfileType = 'txt'; % 'dat' or 'txt'
    params.LOScoordType = 'geographic'; % 'cartesian' or 'geographic'. Default: geographic
    params.LOSuseWeights = 0; % whether to use weights column from the InSAR resampling
    
    %%% Fault (source) information %%%
    params.sourceTypes      = {'Static_Halfspace_Fault_Source'}; % cell list of fault/source types (matching the corresponding Object filename)
    % cell list of source files to read, 1-to-1 correspondence to the above list
    params.sourceFilenames  = {'faults/qiu+15_1'};
    % Static_Halfspace_Fault_Source parameters
    % required:
    params.unicyclePath = '/Users/elindsey/Dropbox/code/geodesy/unicycle';
    % greensType can be any valid unicycle greens object. currently 'okada92' or 'nikkhoo15'
    params.greensType = 'nikkhoo15';
    % optional:
    params.faultOptions = 'rakeCoordinates';
    %params.faultOptions = 'rakeFixed';
    params.rakeFile = 'faults/nepal_approx_rake.dat';
    
    
    %%% Jointinv inversion parameters %%%
    params.inversionType = 'lsqlin'; %options so far: lsqnonneg, backslash
    params.lsqlinOptions = optimoptions('lsqlin','MaxIter',1e12,'display','off');
    params.bounds ={{  {'rake', 0, 10}, {'rakePerp', -1, 1} }}; %currently defined as one cell list for dip and strike variables (e.g. {'dip', min, max}), with one pair for each source.
    params.smoothingType = {{ {'laplacian', 0} , {'value', 'both', 0} }}; % need one cell list for each source. Current options: {'laplacian_1d', 0}, {'laplacian', 0}, {'stressKernel', 0}, {'sum', 'component', target_sliprate},  {'value', 'component', target_sliprate}
    params.nNearest = 4; % for Laplacian smoothing - how many nearest neighbors to include. Standard is 4, but any number can be used - distance weighting is applied
    params.smoothingWeights = {{ 1e-3 , 5e-3}};
    params.stressKernelFolder = 'faults';

    params.verbose = 0;
end
