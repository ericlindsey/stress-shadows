function userParams = exp_40()

% This script configures jointinv2 for:
% - 3D test with backslip
% note: defaults are set in the script set_jointinv_defaults.m


%%% Dataset information %%%
% cell list of dataset types (matching the corresponding Object filename)
userParams.datasetTypes     = {'Static_GPS_Dataset'};
% cell list of data files to read, 1-to-1 correspondence to the above list
userParams.datasetFilenames = {'test_data_3d/gps_synthetic.dat'};
%
% Static_GPS_Dataset parameters
% required:
userParams.lon0 = 0; % these set the polyconic projection origin for the cartesian system used in modeling
userParams.lat0 = 0;
userParams.fileType = 'dat'; % 'dat', 'vel', or 'txt'
% optional:
userParams.coordType = 'cartesian'; % 'cartesian' or 'geographic'. Default: geographic
userParams.dataComponents = [1,2,3]; %here we keep only [East,North] comoponents of the data. Options: single component, or list from [1,2,3]. Default: all 3


%%% Fault (source) information %%%
% cell list of fault/source types (matching the corresponding Object filename)
userParams.sourceTypes      = {'Static_Halfspace_Fault_Source'};
% cell list of source files to read, 1-to-1 correspondence to the above list
userParams.sourceFilenames  = {'test_data_3d/ramp_3d.seg'};
%
% Static_Halfspace_Fault_Source parameters
% required:
userParams.unicyclePath = '/Users/elindsey/Dropbox/code/geodesy/unicycle';
% greensType can be any valid unicycle greens object. currently 'okada92' or 'nikkhoo15'
userParams.greensType = 'okada92';
% optional:
%userParams.shearModulus = 30e3; % All stress units in Unicycle are MPa. Default: 30e3
%userParams.poissonsRatio = 0.25; % Default: 0.25
%userParams.slipComponents = [2]; % list of components to keep - either 1 (strike), 2 (dip), or both [1,2]. Default: both


% Inversion parameters
userParams.inversionType = 'lsqlin'; %options so far: lsqnonneg, backslash
userParams.lsqlinOptions = optimoptions('lsqlin','MaxIter',1e5);

userParams.bounds = {{ {'dip', 0, 100/cosd(10)}, {'strike', -35, 35} }}; %currently defined as one cell list for dip and strike variables (e.g. {'dip', min, max}), with one pair for each source.

userParams.smoothingType = {{ {'stressKernel', 0}, {'sum','dip', 0/cosd(10)} }}; % need one cell list for each source. Current options: {'stressKernel',0}, {'sum','component',target_sliprate}
userParams.smoothingWeights = {{ 1e-3, 1e-3 }};
userParams.stressKernelFolder = 'test_data_3d';
%userParams.constraintType = {{ {'positiveStress', 200e3, 'test_data_3d/rake_angle.dat' } }}; % Current options: {'positiveStress',depth, Vds or file containing [Vss, Vds] for each patch}
%userParams.constraintType = {{ {'positiveStress', 2000e3, 0} }}; % Current options: {'positiveStress',depth, Vds or file containing [Vss, Vds] for each patch}

userParams.verbose = 0;

end