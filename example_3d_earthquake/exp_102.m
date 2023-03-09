function params = exp_102()

    % This parameter file configures jointinv for:
    % 3D slip model

    %%% Dataset information %%%
    params.datasetTypes     = {'Static_GPS_Dataset'}; % cell list of dataset types (matching the corresponding Object filename)
    params.datasetFilenames = {'./data/synthetic_gps.dat'}; % cell list of data files to read, 1-to-1 correspondence to the above list
    % Static_GPS_Dataset parameters
    % required:
    params.lon0 = 0; % these set the polyconic projection origin for the cartesian system used in modeling
    params.lat0 = 0;
    params.fileType = 'dat'; % 'dat', 'vel', or 'txt'
    % optional:
    params.coordType = 'cartesian'; % 'cartesian' or 'geographic'. Default: geographic
    params.dataComponents = [1,2,3]; %Options: single component, or list from [1,2,3]. Default: all 3

    %%% Fault (source) information %%%
    params.sourceTypes      = {'Static_Halfspace_Fault_Source'}; % cell list of fault/source types (matching the corresponding Object filename)
    % cell list of source files to read, 1-to-1 correspondence to the above list
    params.sourceFilenames  = {'./faults/ramp_3d.seg'};
    % Static_Halfspace_Fault_Source parameters
    % required:
    params.unicyclePath = '/Users/elindsey/Dropbox/code/geodesy/unicycle';
    % greensType can be any valid unicycle greens object. currently 'okada92' or 'nikkhoo15'
    params.greensType = 'okada92';
    % optional:
    %params.shearModulus = 30e3; % All stress units in Unicycle are MPa. Default: 30e3
    %params.poissonsRatio = 0.25; % Default: 0.25
    params.slipComponents = 2; % list of components to keep - either 1 (strike), 2 (dip), or both [1,2]. Default: both

    %%% Jointinv inversion parameters %%%
    params.inversionType = 'lsqlin'; %options so far: lsqnonneg, backslash
    params.lsqlinOptions = optimoptions('lsqlin','MaxIter',1e12,'display','off');
    params.bounds = {{ {'dip', 0, 10} }}; %currently defined as one cell list for dip and strike variables (e.g. {'dip', min, max}), with one pair for each source.
    params.smoothingType = {{ {'laplacian', 0} }}; % need one cell list for each source. Current options: {'laplacian_1d', 0}, {'laplacian', 0}, {'stressKernel', 0}, {'sum', 'component', target_sliprate}
    params.nNearest = 4; % for Laplacian smoothing - how many nearest neighbors to include. Standard is 4, but any number can be used - distance weighting is applied
    params.smoothingWeights = {{ 1e-2 }};
    params.stressKernelFolder = 'faults';

    params.verbose = 0;
end
