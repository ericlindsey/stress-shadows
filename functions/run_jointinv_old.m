%
% This is a sample script to help us plan the structure of the
% code for the jointinv version update.
%
% The structure is based on the existing RUN_JOINT_INV.m script and the
% functions it calls, but an attempt has been made to unify the function
% nomenclature and simplify the function calls. 
%
% Philosophy: keep jointinv as general as possible.
%
% Goal: find the least-squares solution to an equation G*m = d, subject to:
% regularization: L*m ~= k, and
% inequality constraints: A*m >= b.
% The matrices G, m, d, L, k, A, and b should have as little assumed
% structure as possible - this is all to be determined at run-time.
% 
% The program is implemented in an object-oriented manner to simplify but
% also constrain the execution.
%
% First version: Eric Lindsey, Nov. 2017
%

% notes: use something besides 'project' for the object name
% use lower camel case for variables
% use lower underscore case for functions
% use upper underscore case for classes

addpath('functions')

%% General parameters

% create the master object
project = Jointinv_Project(234);

%%
% set userParams - these are global parameters
project.userParams.lat0 = -1.0386;
project.userParams.lon0 = 98.9809;
project.userParams.greens_type = 'okada92';
project.userParams.shearModulus = 30e3;
project.userParams.poissonsRatio = 0.25;

% other ideas:
% project.userParams.eqStartTime = '2016-04-27 04:47:15.5';
% project.userParams.velocityModel = 'prem';
% etc.

%% Data lists

% data types specifies the types of data for each input file - ('Static_GPS', 'Teleseismic', etc.)
% each option must correspond to an object inheriting the 'Dataset' type
datasetTypes={'Static_GPS_Dataset'};
% lists the data files to be read
datasetInputs={'test_data/SuGAr_Mentawai.vel'};

project.add_datasets(datasetTypes, datasetInputs);


%% Model lists

% list the geometry or model types ('Planar_Fault', 'Planar_Fault_Layered', 'Mogi', etc.)
source_types={'Static_Halfspace_Fault_Source'};
% list the input files or give fault parameters
source_inputs={'test_data/Ramp2d.seg'};
% for a planar fault, this could be a simple list of [xc,yc,zc,str,dip,L,W,nL,nW]
%source_input={[0,-100e3,10,0,10,200,100,10,10]};


project.add_sources(source_types, source_inputs);


% smoothing types ('None', 'Laplacian', 'Stress', etc.) - can be whatever is implemented by your source type.
smoothing_types={'Laplacian'};
% smoothing parameters, specific to the smoothing type
smoothing_inputs={200};


%% Inversion: to be implemented...


%project.run_inversion();

%make_plots(project);

