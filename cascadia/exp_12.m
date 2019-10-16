function userParams = exp_12()

% This script configures jointinv2 for:
% - Cascadia inversion with backslip
% - With stress constraints
% - Rake coordinates with bounds on rake direction, slip penalty on rake
% - McCaffrey data with 'WWW' rake model and block corrections

% note: defaults are set in the script set_jointinv_defaults.m
%%% Dataset information %%%
% cell list of dataset types (matching the corresponding Object filename)
userParams.datasetTypes     = {'Static_GPS_Dataset'};
% cell list of data files to read, 1-to-1 correspondence to the above list
userParams.datasetFilenames = {'gps/mccaffrey_gps_www_corrected.dat'};
%userParams.datasetFilenames = {'gps/michel_gps.dat'};
%
% Static_GPS_Dataset parameters
% required:
userParams.lon0 = 237; % these set the polyconic projection origin for the cartesian system used in modeling
userParams.lat0 = 45;
userParams.fileType = 'dat'; % 'dat', 'vel', or 'txt'
% optional:
userParams.coordType = 'cartesian'; % 'cartesian' or 'geographic'. Default: geographic
userParams.dataComponents = [1,2]; %e.g. [1,2] will keep only [East,North] comoponents of the data. Options: single component, or list from [1,2,3]. Default: all 3
userParams.minGPSError = 0.5; %default: 0.5 (mm/yr)
%%% Fault (source) information %%%
% cell list of fault/source types (matching the corresponding Object filename)
userParams.sourceTypes      = {'Static_Halfspace_Fault_Source'};
% cell list of source files to read, 1-to-1 correspondence to the above list
userParams.sourceFilenames  = {'faults/cascadia_100km_xyz'};
%
% Static_Halfspace_Fault_Source parameters
% required:
userParams.unicyclePath = '/Users/elindsey/Dropbox/code/geodesy/unicycle';
% greensType can be any valid unicycle greens object. currently 'okada92' or 'nikkhoo15'
userParams.greensType = 'nikkhoo15';
% optional:
%userParams.shearModulus = 30e3; % All stress units in Unicycle are MPa. Default: 30e3
%userParams.poissonsRatio = 0.25; % Default: 0.25
%userParams.slipComponents = [2]; % list of components to keep - either 1 (strike), 2 (dip), or both [1,2]. Default: both
userParams.faultOptions = 'rakeCoordinates'; %options: 'rakeFixed', 'rakeCoordinates'

% Inversion parameters
userParams.inversionType = 'lsqlin'; %options: backslash, lsqnonneg, lsqlin
userParams.lsqlinOptions = optimoptions('lsqlin','MaxIter',1e5); %create a matlab-defined optimoptions() object matching your inversionType here, if desired

userParams.bounds = {{ {'rakeSlipRate'}, {'rakePerp', -5, 5} }}; %currently defined as one cell list for dip and strike variables (e.g. {'dip', min, max}), with one pair for each source.

userParams.smoothingType = {{ {'laplacian', 0}, {'sum', 'rake', 0 } }}; %, {'value', 'rakePerp', 0} }}; % {'value', 'rakePar', 0}, {'value', 'rakePerp', 0} }}; % , {'sum','strike', 0}  }}; % need one cell list for each source. Current options: {'stressKernel',0}, {'sum','component',target_sliprate}
userParams.smoothingWeights = {{ 0.6702, 10^-8.05 }};
userParams.stressKernelFolder = 'kernels';
userParams.nNearest = 4; % for Laplacian smoothing
userParams.constraintType = {{ {'positiveStress', 30e3, 'rake'} }}; % Current options: {'positiveStress',depth, Vds}
%userParams.constraintType = {{ {'positiveStress', 200e3, 'rake'}, {'rakeSlipRate'},{'rakePerpendicularSlipRate',1}, {'negativeRakePerpendicularSlipRate',-1} }}; % Current options: {'positiveStress',depth, Vds}
%userParams.constraintType = {{ {'rakeSlipRate'} }}; % Current options: {'positiveStress',depth, Vds}

userParams.rakeFile = 'faults/cascadia_100km_rake_kreemer_li_corrected.dat';
%userParams.rakeFile = 'faults/cascadia_100km_rake_schmalzle.dat';

userParams.verbose = 0;

end
