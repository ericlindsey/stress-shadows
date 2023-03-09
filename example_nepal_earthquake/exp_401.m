function params = exp_401()

    % This parameter file configures jointinv for:
    % GPS-only model for Nepal

    %%% Dataset information %%%
    params.datasetTypes     = {'Static_GPS_Dataset'}; % cell list of dataset types (matching the corresponding Object filename)
    params.datasetFilenames = {'data/nepal_aria_gps_formatted.dat'}; % cell list of data files to read, 1-to-1 correspondence to the above list
    % Static_GPS_Dataset parameters
    % required:
    params.lon0 = 0; % these set the polyconic projection origin for the cartesian system used in modeling
    params.lat0 = 0;
    params.fileType = 'dat'; % 'dat', 'vel', or 'txt'
    % optional:
    params.coordType = 'cartesian'; % 'cartesian' or 'geographic'. Default: geographic
    params.dataComponents = [1,2,3]; %here we keep only [East] comoponents of the data. Options: single component, or list from [1,2,3]. Default: all 3

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
    params.smoothingType = {{ {'laplacian', 0} }}; % , {'value', 'both', 0} }}; % need one cell list for each source. Current options: {'laplacian_1d', 0}, {'laplacian', 0}, {'stressKernel', 0}, {'sum', 'component', target_sliprate},  {'value', 'component', target_sliprate}
    params.nNearest = 4; % for Laplacian smoothing - how many nearest neighbors to include. Standard is 4, but any number can be used - distance weighting is applied
    params.smoothingWeights = {{ 1e-2 }};
    params.stressKernelFolder = 'faults';

    params.verbose = 0;
    
end
