classdef Jointinv < handle
    
    % this class contains the generic methods and properties common to all
    % geophysical inversion types. The methods correctly assemble the 
    % matrices given a multitude of different data types
    %
    % Latest version: Eric Lindsey, June 2019
    
    properties
        
        expNumber        % unique identifier specifying the filename
                         %  containing the experiment-specific parameters
        datasets         % cell list of dataset objects
        sources          % cell list of source (model) objects
        userParams       % struct containing all other shared information
                         %  that may vary with each dataset and model type.
                         %  The only required elements are 'datasetType',
                         %  'datasetFilename', 'sourceType' and
                         %  'sourceFilename'. May also contain zero times /
                         %  locations, material properties, dataset- and
                         %  model-specific visualization options, inversion
                         %  options, etc.
        dataVector       % combined data vector from all datasets (currently, this duplicates memory...)
        modelVector      % combined model vector from all fault sources (currently, this duplicates memory...)
        predVector       % designMatrix * modelVector
        designMatrix     % Green's function matrix
        dataCovarianceMatrix % Data covariance matrix, block diagonal for independent datasets
        smoothingMatrix  % for regularization, left side
        smoothingVector  % for regularization, right side
        smoothingLengths % keep track of the number of smoothing parameters
        smoothingWidths  % keep track of the number of smoothing parameters
        constraintMatrix % for inequality constraints, left side
        constraintVector % for inequality constraints, right side
        lowerBounds      % for bounds
        upperBounds      % for bounds
        chi2             % misfit statistic, not normalized by N
        
    end % properties
    
    methods
        
        % object constructor.
        function obj = Jointinv(expNumber)
            disp(['Creating jointinv experiment: ' num2str(expNumber) '.']) 
            obj.expNumber = expNumber;
            obj.datasets={};
            obj.sources={};
            obj.dataVector=[];
            obj.modelVector=[];
            obj.predVector=[];
            obj.designMatrix=[];
            obj.dataCovarianceMatrix=[];
            obj.smoothingMatrix=[];
            obj.smoothingLengths={};
            obj.smoothingWidths={};
            obj.smoothingVector=[];
            obj.constraintMatrix=[];
            obj.constraintVector=[];
            obj.lowerBounds=[];
            obj.upperBounds=[];
            obj.chi2=[];
        end
        
        
        % convenience function to run all setup functions
        function run_setup(obj, varargin)
            disp(['Setting up jointinv experiment: ' num2str(obj.expNumber) '.']) 
            %read user parameters and dataset lists from 'exp_#.m'
            % experiment file may optionally be specified
            if (nargin == 2)
                obj.read_user_params(varargin{1});
            else
                obj.read_user_params();
            end
            % create the datasets
            disp('Loading datasets.')
            obj.add_datasets();
            obj.calc_data_covariance_matrix(); %assemble covariance matrices as block diagonal
            disp('')
            % create the sources
            disp('Loading deformation sources.')
            obj.add_sources();
            disp('')
            % compute the G matrix
            disp('Calculating design matrix.')
            obj.calc_design_matrix();
            disp('')
            % compute the Laplacian or stress kernel, depending on options
            if isfield(obj.userParams, 'smoothingType')
                disp('Calculating smoothing matrix.')
                obj.calc_smoothing_matrix();
            end
            % setup any specified model constraints, depending on options
            if isfield(obj.userParams, 'constraintType')
                disp('Creating constraints.')
                obj.calc_constraint_matrix();
                disp('')
            end
            % setup any specified model bounds, depending on options
            if isfield(obj.userParams, 'bounds')
                disp('Creating bounds.')
                obj.calc_bounds();
                disp('')
            end
            disp(['Done setting up jointinv experiment ' num2str(obj.expNumber) '!'])
            disp(obj)
        end % function run_setup

        % read information from the experiment file
        function read_user_params(obj, varargin)
            % experiment file may optionally be specified
            if (nargin == 2)
                expfile = varargin{1};
                disp(['Reading user-specified parameter file: ' expfile '.'])   
            else
                expfile = ['exp_' num2str(obj.expNumber)];
                disp(['Reading default parameter file: ' expfile '.'])    
            end
            if exist(expfile,'file')
                % The experiment file is a function. This function should return
                % one object, a struct containing the user-specified params
                expfunction = str2func(expfile);
                obj.userParams = expfunction();
            else
                error(['Error: parameter file ' expfile ' not found.'])
            end
            % set defaults
            obj.userParams = set_jointinv_defaults(obj.userParams);
        end % function read_user_params
        
        
        % load a list of datasets
        function add_datasets(obj)
            % require that input list lengths are the same
            assert(length(obj.userParams.datasetTypes) == length(obj.userParams.datasetFilenames), 'Error: length of dataset type and input cells must be the same in experiment file.');
            for i=1:length(obj.userParams.datasetTypes)
                % test whether this file was loaded already
                load_this=1;
                for j=1:length(obj.datasets)
                    if strcmp(obj.userParams.datasetFilenames{i}, obj.datasets{j}.fileName)
                        disp(['File ' obj.userParams.datasetFilenames{i} ' was loaded already, skipping.'])
                        load_this=0;
                        break
                    end
                end
                if load_this==1
                    disp(['loading dataset ' obj.userParams.datasetFilenames{i}])
                    % convert the string datasetType into a callable function
                    dataset_handle = str2func(obj.userParams.datasetTypes{i});
                    % call this function, passing the associated data file and any other userParams that were set
                    obj.datasets{end+1} = dataset_handle(obj.userParams.datasetFilenames{i}, obj.userParams);
                    obj.datasets{end}.fileName = obj.userParams.datasetFilenames{i};
                    % copy the data into the jointinv scenario dataVector
                    obj.dataVector = [obj.dataVector; obj.datasets{end}.dataVector];
                    % print a summary of what was loaded
                    disp(obj.datasets{end})
                end
            end    
        end % function add_datasets
        
        
        % load a list of sources
        function add_sources(obj)
            % require that input list lengths are the same
            assert(length(obj.userParams.sourceTypes) == length(obj.userParams.sourceFilenames), 'Error: length of source types and input cells must be the same in experiment file.');
            for i=1:length(obj.userParams.sourceFilenames)
                % test whether this file was loaded already
                load_this=1;
                for j=1:length(obj.sources)
                    if strcmp(obj.userParams.sourceFilenames{i}, obj.sources{j}.fileName)
                        disp(['File ' obj.userParams.sourceFilenames{i} ' was loaded already, skipping.'])
                        load_this=0;
                        break
                    end
                end
                if load_this==1
                    % convert the string sourceType into a callable function
                    source_handle = str2func(obj.userParams.sourceTypes{i});
                    % call this function, passing the associated source file and any other userParams that were set
                    obj.sources{end+1} = source_handle(obj.userParams.sourceFilenames{i}, obj.userParams);
                    % save the filename that this object was created from
                    obj.sources{end}.fileName = obj.userParams.sourceFilenames{i};
                    % copy the model vector into the jointinv scenario modelVector
                    obj.modelVector = [obj.modelVector; obj.sources{end}.modelVector];
                    % print a summary of what was loaded
                    disp(obj.sources{end})
                end
            end
        end % function add_sources
        
        
        function calc_data_covariance_matrix(obj)
            % setup data covariance matrix
            obj.dataCovarianceMatrix=[];
            for i=1:length(obj.datasets)
                obj.dataCovarianceMatrix=blkdiag(obj.dataCovarianceMatrix,obj.datasets{i}.covarianceMatrix);
            end
        end % function calc_data_covariance_matrix
        
        
        function calc_design_matrix(obj)
            %loop over each source and dataset
            thisCol=1;
            for i=1:length(obj.sources)
                thisRow=1;
                for j=1:length(obj.datasets)
                    % call the greens function computation
                    [thisG,thisModelVector] = obj.sources{i}.calc_design_matrix(obj.datasets{j}, obj.userParams, j);
                    nextRow = thisRow + size(thisG,1);
                    nextCol = thisCol + size(thisG,2);
                    %assign the correct block of the design matrix.
                    %unfortunately this re-dimensions the matrix each time.
                    obj.designMatrix(thisRow:nextRow-1, thisCol:nextCol-1) = thisG;
                    obj.modelVector(thisCol:nextCol-1,1) = thisModelVector;
                    %track the location to place the next block
                    thisRow = nextRow;
                end
                thisCol = nextCol;
            end            
        end % function calc_design_matrix

        
        function calc_smoothing_matrix(obj)
            %loop over each source, and over list of smoothing kernels for that source
            thisCol=1;
            for i=1:length(obj.sources)
                thisRow=1;
                for j=1:length(obj.userParams.smoothingType{i})
                    % call the actual smoothing matrix computation
                    [thisL,thisLvector] = obj.sources{i}.calc_smoothing_matrix(obj.userParams.smoothingType{i}{j}, obj.userParams);
                    nextRow = thisRow + length(thisLvector);
                    nextCol = thisCol + size(thisL,2);
                    obj.smoothingLengths{i}{j}=length(thisLvector);
                    obj.smoothingWidths{i}{j}=size(thisL,2);
                    %assign the correct block of the smoothing matrix and vector
                    %unfortunately this re-dimensions the matrix each time.
                    obj.smoothingMatrix(thisRow:nextRow-1, thisCol:nextCol-1) = thisL;
                    obj.smoothingVector(thisRow:nextRow-1,1) = thisLvector;
                    %track the location to place the next block
                    thisRow = nextRow;
                end
                thisCol = nextCol;
            end            
        end % function calc_smoothing_matrix


        function calc_constraint_matrix(obj)
            %loop over each source and list of constraints
            thisCol=1;
            for i=1:length(obj.sources)
                thisRow=1;
                for j=1:length(obj.userParams.constraintType{i})
                    % call the actual smoothing matrix computation
                    [thisK,thisKvector] = obj.sources{i}.calc_constraint_matrix(obj.userParams.constraintType{i}{j}, obj.userParams);
                    nextRow = thisRow + length(thisKvector);
                    nextCol = thisCol + size(thisK,2);
                    %assign the correct block of the constraint matrix and vector
                    %unfortunately this re-dimensions the matrix each time.
                    obj.constraintMatrix(thisRow:nextRow-1, thisCol:nextCol-1) = thisK;
                    obj.constraintVector(thisRow:nextRow-1,1) = thisKvector;
                    %track the location to place the next block
                    thisRow = nextRow;
                end
                thisCol = nextCol;
            end
        end % function calc_constraint_matrix
        
        function calc_bounds(obj)
            %loop over each source
            thisRow=1;
            for i=1:length(obj.sources)
                [l,u] = obj.sources{i}.calc_bounds(obj.userParams.bounds{i}, obj.userParams);
                nextRow = thisRow + length(l);
                obj.lowerBounds(thisRow:nextRow-1) = l;
                obj.upperBounds(thisRow:nextRow-1) = u;
                thisRow = nextRow;
            end
        end % function calc_bounds
        
        
        function run_inversion(obj)
            
            % Apply weights: multiply G and d by W (W = cov^-0.5)
            augmentedDesignMatrix = obj.dataCovarianceMatrix^(-0.5) * obj.designMatrix;
            augmentedDataVector = obj.dataCovarianceMatrix^(-0.5) * obj.dataVector;
            
            % use in case weight matrix breaks (e.g. zero errors?):
            %augmentedDesignMatrix = obj.designMatrix;
            %augmentedDataVector = obj.dataVector;
            
            % setup smoothing matrix and add to end of design matrix and
            % data vector
            if ~isempty(obj.smoothingMatrix)
                %loop over each source, and over list of smoothing kernels for that source
                addRow=size(obj.designMatrix,1);
                thisCol=1;
                for i=1:length(obj.sources)
                    thisRow=1;
                    for j=1:length(obj.userParams.smoothingType{i})
                        nextRow = thisRow +  obj.smoothingLengths{i}{j};
                        nextCol = thisCol + obj.smoothingWidths{i}{j};
                        %assign the correct block of the smoothing matrix and vector to the augmented matrices.
                        augmentedDesignMatrix(addRow+thisRow:addRow+nextRow-1, thisCol:nextCol-1) = obj.userParams.smoothingWeights{i}{j} * obj.smoothingMatrix(thisRow:nextRow-1, thisCol:nextCol-1);
                        augmentedDataVector(addRow+thisRow:addRow+nextRow-1,1) = obj.userParams.smoothingWeights{i}{j} * obj.smoothingVector(thisRow:nextRow-1,1);
                        %track the current location
                        thisRow = nextRow;
                    end
                    % note, column skip here is not the same as smoothingWidths, if some parameters
                    % like trend do not have smoothing applied.
                    thisCol = thisCol + length(obj.sources{i}.modelVector); 
                end            
            end
            
            if(obj.userParams.verbose), disp(['Running inversion: using method ' obj.userParams.inversionType]), end
            
            if strcmp(obj.userParams.inversionType, 'backslash')
                assert(isempty(obj.constraintMatrix), 'Error: cannot use constraints with this inversion type.');
                obj.modelVector = augmentedDesignMatrix \ augmentedDataVector;

            elseif strcmp(obj.userParams.inversionType, 'lsqnonneg')
                assert(isempty(obj.constraintMatrix), 'Error: cannot use constraints with this inversion type.');
                obj.modelVector = lsqnonneg(augmentedDesignMatrix, augmentedDataVector);
                
            elseif strcmp(obj.userParams.inversionType, 'lsqlin')
                % lsqlin arguments: G, d, ineqC, ineqCd, [eqC], [eqCd], lowerB, upperB, [X0], options
                obj.modelVector = lsqlin(augmentedDesignMatrix, augmentedDataVector, obj.constraintMatrix, obj.constraintVector, [], [], obj.lowerBounds, obj.upperBounds, [], obj.userParams.lsqlinOptions);
            end
            
            % compute predicted data
            obj.predVector = obj.designMatrix * obj.modelVector; 
            %distribute the predicted vector back to the datasets
            thisRow=1;
            for i=1:length(obj.datasets)
                nextRow = thisRow + length(obj.datasets{i}.dataVector);
                obj.datasets{i}.predVector = obj.predVector(thisRow:nextRow-1);
                thisRow = nextRow;
            end
            % compute the chi-squared misfit statistic
            obj.chi2 = ((obj.predVector - obj.dataVector)'*inv(obj.dataCovarianceMatrix)*(obj.predVector - obj.dataVector)); 
    
            %distribute the model vector back to the sources
            thisRow=1;
            for i=1:length(obj.sources)
                nextRow = thisRow + length(obj.sources{i}.modelVector);
                obj.sources{i}.modelVector = obj.modelVector(thisRow:nextRow-1);
                thisRow = nextRow;
            end
            if(obj.userParams.verbose), disp('Inversion finished. Model results are in scenario.modelVector; predicted data in scenario.predVector'), end
        end  % function run_inversion
        
    end % methods
    
end % classdef Jointinv
