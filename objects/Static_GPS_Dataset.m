classdef Static_GPS_Dataset < Jointinv_Dataset
    % Static_GPS_Dataset defines the basic general-purpose GPS offsets.
    % Note it can be equally-well used for GPS velocities, or anything else
    % with 3 components of observed data at a single point.
    
    properties
        
% Inherited from Jointinv_Dataset:
%         
%         coordinates      % n x [2, 3, or 4] - independent variables (locations, times) for each data point
%         
%         dataVector       % N x 1 - column vector containing the dependent variables (data)
%                          % the data should be unraveled by point, not by component.
%                          % eg. (pt1_E, pt1_N, pt1_U, pt2_E, pt2_N, pt2_U, ...)
%                          % note that N = n*k
%         
%         covarianceMatrix % N x N - covariance matrix containing the data uncertainties/weights
%         
%         numComponents    % dimension k of the data (note that data is
%                          % implicitly unraveled; this is needed to determine
%                          % how the coordinates match to the data)
        
        % Properties specific to Static_GPS_Dataset
        lat0
        lon0
        name
        lat
        lon
    end
    
    methods
        
        %
        % constructor - this is the only method for this class
        %
        function obj = Static_GPS_Dataset(dataFile, userParams)
            %
            % userParams values required: lon0, lat0, fileType, coordType
            %
            % fileType may be 'vel' or 'dat'
            %
            % coordType may be 'cartesian' or 'geographic' and controls
            % whether the conversion using latlon_to_xy_polyconic is used.
            %
            assert(nargin==2, 'Error: Static_GPS_Dataset did not receive any arguments');
            assert(logical(exist(dataFile,'file')), ['Error: Static_GPS_Dataset cannot find file: ', dataFile]);
            
            obj.numComponents=length(userParams.dataComponents);  % number of dependent data components (GPS data is assumed to contain 3D displacements)
            obj.name={};
                        
            if strcmp(userParams.fileType,'dat')
                % data is a matlab-readable .dat file of size N x 6
                % containing simple locations and displacements:
                % [lon,lat,hgt,uE,uN,uZ,sigE,sigN,sigZ,[wgt]]
                
                data=load(dataFile);
                
                %each point has 3 components
                numData=3*size(data,1);
                
                % load coordinates
                obj.lon = data(:,1);
                obj.lat = data(:,2);
                height = data(:,3);
                
                % reshape data and uncertainties
                % note, when reshaping, the transpose (') is needed to make sure
                % the data column keeps E,N,U components for each point together.
                obj.dataVector=reshape(data(:,4:6)',numData,1);
                
                % weights are defined by sigma^2/weight. These become the
                % diagonal elements of the covariance matrix
                if size(data,2) == 10
                    weights=reshape(repmat(data(:,10),3)',numData,1);
                else
                    weights=ones(numData,1);
                end
                variances=max(userParams.minGPSError,reshape(data(:,7:9)',numData,1)).^2./weights;
                
                % generate default names
                for i=1:length(height)
                    obj.name{i}=sprintf('%04d',i);
                end
                
            elseif strcmp(userParams.fileType,'txt')
                % data is a matlab-readable .txt file of size N x 8 with
                % named columns (in any order):
                % required: NAME Long. Lat. V_E V_N  Sig_E Sig_N  rho
                % optional: [V_Z] [Sig_Z] [Elev] [weight]

                data=readtable(dataFile);
               
                % load coordinates
                obj.name = data.NAME;
                obj.lon = data.Long_;
                obj.lat = data.Lat_;
                if ismember('Elev', data.Properties.VariableNames)
                    height = data.Elev;
                else
                    height = 0*data.Long_;
                end
                
                %if V_Z does not exist, use zeros
                if ~ismember('V_Z', data.Properties.VariableNames)
                    data.V_Z = 0*data.V_E;
                    data.Sig_Z = 0*data.Sig_E;
                end
                
                %each point has 3 components
                numData=3*size(data.V_E,1);
                
                % reshape data and uncertainties
                % note, when reshaping, the transpose (') is needed to make sure
                % the data column keeps E,N,U components for each point together.
                obj.dataVector=reshape([data.V_E data.V_N data.V_Z]',numData,1);
                
                % weights are defined by sigma^2/weight. These become the
                % diagonal elements of the covariance matrix
                if ismember('weight', data.Properties.VariableNames)
                    weights=reshape(repmat(data.weight,3)',numData,1);
                else
                    weights=ones(numData,1);
                end
                variances=max(userParams.minGPSError,reshape([data.Sig_E data.Sig_N data.Sig_Z]',numData,1)).^2./weights;
                
            elseif strcmp(userParams.fileType,'vel')
                % second case: data is a file with .vel format:
                % [Sta Lon Lat Height OE ON OU ErE ErN ErU Weight T0 T0D T1 T1D Cne]
                % This code is modified from GTDef
                %
                % The full 'vel' format is given in the following sample header:
                %
                % # EQ coseismic deformation extracted from <synthetic_model>
                % # 1    2   3   4    5   6    7    8    9          10    11  12
                % # Year Mon Day Hour Min Sec  Lon  Lat  Depth[km]  Mag   ID  Decyr
                % # 2016 04 27 04 47 15.5000  98.9809   -1.0386  30.28  8.90    201604270000     2016.320128
                % # 1   2   3   4      5  6  7  8   9   10  11      12 13  14 15    16
                % # Sta Lon Lat Height OE ON OU ErE ErN ErU Weight  T0 T0D T1 T1D   Cne
                % # Height [m] Disp [mm] error [mm]
                % ABGS 99.387520914 0.220824642 236.2533 -521.5 -1106.8 -195.8 1.0 1.0 3.0 1.0 20160427 2016.32030 20160427 2016.32030 0.0
                % BABI 96.67120337 2.11679941 -21.3099 -7.0 -11.0 -48.0 1.0 1.0 3.0 1.0 20160427 2016.32030 20160427 2016.32030 0.0
                %
                % note that only the first 11 columns are read in this code.
                
                obj.lat=[];
                obj.lon=[];
                height=[];
                obj.dataVector=[];
                variances=[];
                
                ln=0;
                data_fid = fopen(dataFile);
                while(1)
                   
                    % read in one line
                    tline = fgetl(data_fid);
                    ln=ln+1;

                    % test if it is the end of file; exit if yes
                    if ischar(tline)~=1, break; end

                    % read the 1st term: either a '#' or the station name; could be some white-spaces before the 1st term
                    [name,remain] = strtok(tline);		% the default delimiter is white-space

                    % omit '#,%, and blank' comment lines
                    if (strncmpi(name,'#',1)||strncmpi(name,'%',1)||isempty(name)), continue; end

                    % append name
                    obj.name{ln} = name;

                    % read rest of the numeric parameters
                    linedata=str2num(remain); %#ok<ST2NM>
                    obj.lon = [obj.lon; linedata(1)];
                    obj.lat = [obj.lat; linedata(2)];
                    height = [height; linedata(3)];

                    obj.dataVector=[obj.dataVector; linedata(4:6)'];

                    % weights are defined by sigma^2/weight. These become the
                    % diagonal elements of the covariance matrix
                    if length(linedata) >= 10
                        weights=repmat(linedata(10),3,1);
                    else
                        weights=ones(3,1);
                    end
                    variances = [variances; max(userParams.minGPSError,linedata(7:9)').^2./weights];

                end
                fclose(data_fid);
          
            else
                error('Error: Static_GPS_Dataset did not recognize fileType') 
            end
            
            %
            % format and project the coordinates as needed
            %
            
            if isfield(userParams, 'dataComponents')
               %keep only some of the data components (E,N,U) and weights
               obj.dataVector = keep_matrix_rows(obj.dataVector,userParams.dataComponents);
               variances = keep_matrix_rows(variances,userParams.dataComponents);
            end
            
            % TODO: include non-diagonal correlations
            obj.covarianceMatrix = diag(variances);
            
            %set (lat0,lon0) reference point: assumes these are present in userParams.
            obj.lat0 = userParams.lat0;
            obj.lon0 = userParams.lon0;
            
            % convert coordinates to cartesian
            if strcmp(userParams.coordType,'geographic')
                %convert x,y coordinates from lat, lon to meters relative to the reference point
                [x,y] = latlon_to_xy_polyconic(obj.lat, obj.lon, obj.lat0, obj.lon0);
                obj.coordinates=[x*1e3,y*1e3,height];
            elseif strcmp(userParams.coordType,'cartesian') 
                obj.coordinates=[obj.lon,obj.lat,height];  
            else
                error('Error: Static_GPS_Dataset did not recognize coordType') 
            end
            
        end % end constructor
        
    end
    
end

