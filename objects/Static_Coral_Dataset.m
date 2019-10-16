classdef Static_Coral_Dataset < Jointinv_Dataset
    % Static_Coral_Dataset defines the basic Coral uplift-rate dataset.
    % Note it can be equally-well used for anything else
    % with 1 component of observed data at a single point.
    
    properties
        
% Inherited from Jointinv_Dataset:
%         
%         coordinates      % n x [2, 3, or 4] - independent variables (locations, times) for each data point
%         
%         dataVector       % n x 1 - column vector containing the dependent variables (data)
%                          % By definition, there is only one component.
%                          % eg. (pt1_z, pt2_z, pt3_z ...)
%         
%         covarianceMatrix % n x n - covariance matrix containing the data uncertainties/weights
%         
%         numComponents    % dimension k of the data (note that data is
%                          % implicitly unraveled; this is needed to determine
%                          % how the coordinates match to the data)
        
        % Properties specific to Static_Coral_Dataset
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
        function obj = Static_Coral_Dataset(dataFile, userParams)
            %
            % userParams values required: lon0, lat0, fileType, coordType
            %
            % fileType may be 'vel' or 'dat'
            %
            % coordType may be 'cartesian' or 'geographic' and controls
            % whether the conversion using latlon_to_xy_polyconic is used.
            %
            assert(nargin==2, 'Error: Static_Coral_Dataset did not receive any arguments');
            assert(logical(exist(dataFile,'file')), ['Error: Static_Coral_Dataset cannot find file: ', dataFile]);
            
            obj.numComponents=1;  % number of dependent data components (Coral data is assumed to contain 1D displacements)
            obj.name={};
                        
            if strcmp(userParams.fileType,'dat')
                % data is a matlab-readable .dat file of size N x 5
                % containing simple locations and displacements:
                % [lon,lat,hgt,uZ,sigZ,[wgts]]
                data=load(dataFile);
                
                %each point has 1 components
                numData=size(data,1);
                
                % load coordinates
                obj.lon = data(:,1);
                obj.lat = data(:,2);
                height = data(:,3);
                
                % read data
                obj.dataVector=data(:,4);
                
                % weights are defined by sqrt(wgt)/sigma. These become the
                % diagonal elements of the covariance matrix
                if size(data,2) == 6
                    weights=data(:,6);
                else
                    weights=ones(numData,1);
                end
                weightVector=sqrt(weights)./data(:,5);
                
                % generate default names
                for i=1:length(height)
                    obj.name{i}=sprintf('%04d',i);
                end
                
            end
            
            % format and project the coordinates as needed
            
            % TODO: we do not load true correlations (yet)
            obj.covarianceMatrix = diag(weightVector);
            
            %set (lat0,lon0) reference point: assumes these are present in userParams.
            obj.lat0 = userParams.lat0;
            obj.lon0 = userParams.lon0;
            
            % convert coordinates to cartesian
            if strcmp(userParams.coordType,'geographic')
                %convert x,y coordinates from lat, lon to meters relative to the reference point
                [x,y] = latlon_to_xy_polyconic(obj.lat, obj.lon, obj.lat0, obj.lon0);
                obj.coordinates=[x,y,height];
            elseif strcmp(userParams.coordType,'cartesian') 
                % do not convert anything, just use the coordinates that
                % were read in ("lon" = X, "lat" = Y)
                obj.coordinates=[obj.lon,obj.lat,height];  
            else
                error('Error: Static_Coral_Dataset did not recognize coordType') 
            end
            
        end % end constructor
        
    end
    
end

