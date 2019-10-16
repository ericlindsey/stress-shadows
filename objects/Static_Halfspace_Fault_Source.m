classdef Static_Halfspace_Fault_Source < Jointinv_Source
    %Unicycle_Fault_Source provides an interface to greens function and
    % smoothing/stress matrix computations from Unicycle.
    % Valid only for static displacements. Has an option to project onto a
    % Radar LOS or return only a single component
    % 
    % by Eric Lindsey
    % Latest version: June 2019
    
    properties

        % Inherited from Jointinv_Source:
        %         
        % modelVector     % M x 1  - model parameters, unraveled however you like
        
        % Specific to Unicycle_Fault_Source
        %
        earthModel    % Unicycle greens.okada92 object
        geom          % Unicycle geometry.source object
        KK            % stress kernel matrix, [Kss, Kds; Ksd, Kdd]
        Ksn           % stress kernel matrix, strike/normal component
        Kdn           % stress kernel matrix, dip/normal component
        rakeFile      % file with columns: [rake(deg), Vss, Vds, Ve, Vn]. Only first 3 columns are used.
        Rmat          % standard rotation matrix, when using Rake-oriented coordinates
        Vpl           % plate motion rate on each fault patch, computed from rakeFile
    end
    
    methods
        
        % % % % % % % % % % % % % % % % % % %
        %                                   %
        %      c o n s t r u c t o r        %
        %                                   %
        % % % % % % % % % % % % % % % % % % %
        function obj = Static_Halfspace_Fault_Source(patchFile, userParams)
            
            %setup Unicycle path, only if it wasn't already done
            set_unicycle_path(userParams.unicyclePath)
            
            % greens_type can be any valid unicycle greens object.
            % currently this includes 'okada92' or 'nikkhoo15'
            greensHandle = str2func(['unicycle.greens.' userParams.greensType]);
            obj.earthModel = greensHandle(userParams.shearModulus, userParams.poissonsRatio);

            % read the patch file and save as unicycle object
            if strcmp(userParams.greensType,'nikkhoo15')
                %note, in this case we cannot use the 'source' class
                %because that code is old and assumes rectangular patches
                obj.geom = unicycle.geometry.triangleReceiver(patchFile, obj.earthModel);
                % check signs on the cross products
                check_triangleReceiver_signs(obj.geom);
            else
                obj.geom = unicycle.geometry.source(patchFile, obj.earthModel);
            end
            
            % load the rake file, if needed
            if isfield(userParams, 'rakeFile')
                obj.rakeFile = load(userParams.rakeFile);
                [obj.Rmat, obj.Vpl] = get_rake_rotation_matrix(obj.rakeFile);
            end
            
        end
                
        % Design matrix computation - this function name is required
        function [G,modelVector] = calc_design_matrix(obj, dataset, userParams, Idataset)
            
            % keep different components of the matrix, depending on the dataset type
            if isa(dataset,'Static_GPS_Dataset')
                % simplest case - 3 component displacement vector
                coords = [dataset.coordinates(:,1),dataset.coordinates(:,2)];
                G = unicycle_displacement_kernel(obj.geom,coords,userParams.slipComponents,userParams.stressKernelFolder);
                
                if isfield(userParams, 'dataComponents')
                   %keep only some of the data rows (E,N,U) in the matrix
                   G = keep_matrix_rows(G,userParams.dataComponents);
                end
                
            elseif isa(dataset,'Static_Coral_Dataset')
                % next simplest case - 1 component displacement vector
                coords = [dataset.coordinates(:,1),dataset.coordinates(:,2)];
                G = unicycle_displacement_kernel(obj.geom,coords,userParams.slipComponents,userParams.stressKernelFolder);
                
                % keep only the uplift component
                G = keep_matrix_rows(G,3);
                
            elseif isa(dataset,'Static_LOS_Dataset')
                coords = [dataset.coordinates(:,1),dataset.coordinates(:,2)];
                G_3d = unicycle_displacement_kernel(obj.geom,coords,userParams.slipComponents,userParams.stressKernelFolder);
                %project onto radar LOS
                for i = 1:length(dataset.dataVector)
                    G(i,:) = dataset.losVector(i,:) * G_3d(1+3*(i-1):3*i,:);
                end
                
                 % extra options
                if isfield(userParams, 'trendParameters')
                    assert(length(userParams.trendParameters) >= Idataset,'Error: trendParameters list is too short');
                    numtrend=sum([userParams.trendParameters{:}]);
                    Ioffset=sum([userParams.trendParameters{1:Idataset-1}]);
                    M=size(G,2);
                    G=[G,zeros(size(G,1),numtrend)]; %add columns for all offset parameters
                    foundtrend=0;
                    if userParams.trendParameters{Idataset} >= 1
                        disp('Adding one parameter for offset')
                        G(:,M+Ioffset+1)=1; %add a constant offset
                        foundtrend=1;
                    end
                    if userParams.trendParameters{Idataset} >= 3
                        disp('Adding two parameters for x and y trend')
                        G(:,M+Ioffset+2)=dataset.coordinates(:,1); %add x trend
                        G(:,M+Ioffset+3)=dataset.coordinates(:,2); %add y trend
                    end
                    if userParams.trendParameters{Idataset} >= 6
                        disp('Adding three parameters for quadratic trend')
                        G(:,M+Ioffset+4)=dataset.coordinates(:,1).*dataset.coordinates(:,1); %add x.^2 trend
                        G(:,M+Ioffset+5)=dataset.coordinates(:,2).*dataset.coordinates(:,2); %add y.^2 trend
                        G(:,M+Ioffset+6)=dataset.coordinates(:,2).*dataset.coordinates(:,1); %add x.*y trend
                    end
                    
                    if foundtrend==0
                        disp('Did not recognize trendParameters value, doing nothing.')
                    end
                end
                
            else
                disp(['Unknown dataset type ' class(dataset) ' for ' class(obj) ' method get_design_matrix. Returning zeros.'])
                G = zeros(length(dataset.dataVector), 1);
            
            end
            
           
            
            if isfield(userParams, 'faultOptions')
                if strcmp(userParams.faultOptions, 'rakeFixed')
                    % rotate the strike/dip coordinates into the plate-rake
                    % direction and keep only the first half
                    
                    % rotate G matrix into rake direction, keeping only the
                    % rake component
                    G = G * obj.Rmat(:,1:end/2); % is this correct?
                    %G = - G * [diag(cosd(rakeangle)); diag(sind(rakeangle))]; % is this correct?
                elseif strcmp(userParams.faultOptions, 'rakeCoordinates')
                    % rotate the strike/dip coordinates into the plate-rake
                    % direction and keep both coordinates
                    
                    % rotate G matrix
                    G = G * obj.Rmat; % is this correct?
                    %G = G * [diag(sind(rakeangle)), diag(-cosd(rakeangle)); diag(cosd(rakeangle)), diag(sind(rakeangle))]; % is this correct?
                end
            end
            
            % allocate modelVector.
            obj.modelVector = zeros( size(G,2), 1);
            modelVector=obj.modelVector;
            
        end
        
        % Smoothing matrix computation - this function name is required
        function [L,lVector] = calc_smoothing_matrix(obj, smoothingType, userParams)
            disp(['Smoothing type is: ' smoothingType{1} '.'])
            if strcmp(smoothingType{1},'laplacian')
                N = userParams.nNearest;
                disp(['Using nearest ' num2str(N) ' patches to generate finite difference Laplacian for each slip component.'])
                x=obj.geom.xc(:,1);
                y=obj.geom.xc(:,2);
                z=obj.geom.xc(:,3);
                L=compute_laplacian(x,y,z,N);
                if length(userParams.slipComponents) == 2 && ~(isfield(userParams, 'faultOptions') && strcmp(userParams.faultOptions, 'rakeFixed'))
                    % include separate laplacians for each slip component
                    L= [L, 0*L; 
                        0*L, L];
                end
                % adjust units
                L = 1e8 .* L;
                % lVector is the length of L
                lVector=smoothingType{2}*ones(size(L,1),1);
                
            elseif strcmp(smoothingType{1},'laplacian_1d')
                % special, simplified case for 1D faults. The function used for 3D faults doesn't work for this case
                L=compute_laplacian_1d(obj.geom.N);
                if length(userParams.slipComponents) == 2 && ~(isfield(userParams, 'rakeOptions') && strcmp(userParams.rakeOptions, 'rakeFixed'))
                    % include separate laplacians for each slip component
                    L= [L, 0*L; 
                        0*L, L];
                end
                % lVector is the length of L
                lVector=smoothingType{2}*ones(size(L,1),1); 
                
            elseif strcmp(smoothingType{1},'stressKernel')
                % stress kernel matrices are stored as part of this object
                % for re-use in the constraint matrix
                [obj.KK, obj.Ksn, obj.Kdn]=unicycle_stress_kernel(obj.geom, userParams.slipComponents,userParams.stressKernelFolder);
                L=obj.KK;
                
                % extra options
                if isfield(userParams, 'faultOptions')
                    if strcmp(userParams.faultOptions, 'rakeFixed')
                        % rotate L matrix into rake direction, keeping only the first half
                        L = L * obj.Rmat(:,1:end/2);
                    elseif strcmp(userParams.faultOptions, 'rakeCoordinates')
                        % rotate L matrix into rake direction, keeping both
                        L = obj.Rmat' * L * obj.Rmat;
                    end
                end
                % lVector is the length of L
                lVector=smoothingType{2}*ones(size(L,1),1);
                
            elseif strcmp(smoothingType{1},'sum')
                disp(['Creating one smoothing row for specified slip component (', smoothingType{2}, ')  for slip rate sum.'])
                L=zeros(1,length(obj.modelVector));
                if strcmp(smoothingType{2},'rake')
                    % multiply scaling factor by sum of plate rate on all patches
                    %lVector = smoothingType{3} * sum(obj.Vpl)/sqrt(obj.geom.N);
                    lVector = smoothingType{3} * sum(obj.Vpl);
                else
                    % fall back to unit weighting
                    lVector=smoothingType{3}*obj.geom.N;
                end
                if length(userParams.slipComponents) == 1
                    L(:)=1; % in this case, the modelVector contains only one component, so all entries are 1
                else
                    if strcmp(smoothingType{2},'strike')
                        L(1:obj.geom.N)=1;
                    elseif strcmp(smoothingType{2},'dip')
                        L(obj.geom.N+1:2*obj.geom.N)=1;
                    elseif strcmp(smoothingType{2},'both')
                        L(1:end)=1;
                    elseif strcmp(smoothingType{2},'rake')
                        if isfield(userParams, 'faultOptions') && ( strcmp(userParams.faultOptions, 'rakeCoordinates') || strcmp(userParams.faultOptions, 'rakeFixed') )
                            % treat this like the strike case - rake direction is the first coordinate
                            L(1:obj.geom.N)=1;
                        else
                            % not rake coordinates, so convert strike/dip
                            rakeangle=obj.rakeFile(:,1);
                            % rotate slip into rake direction
                            L=[cosd(rakeangle)', sind(rakeangle)'];
                        end
                    elseif strcmp(smoothingType{2}, 'rakeFixed') 
                        L(1:obj.geom.N)=1;
                    else
                        disp('Slip component not recognized. Returning zero row.')
                    end
                end
                
            elseif strcmp(smoothingType{1},'value')
                if ~isnumeric(smoothingType{3}) && strcmp(smoothingType{3}, 'Vpl')
                    lVector = - obj.Vpl; %sign convention for backslip
                elseif isnumeric(smoothingType{3})
                    lVector = smoothingType{3} .* ones(obj.geom.N,1);
                else
                    error('Error: smoothingType{3} must be numeric or "Vpl"')
                end
                
                disp(['Creating penalties for specified slip component (', smoothingType{2}, ') relative to value (', num2str(smoothingType{3}), ').'])
                if length(userParams.slipComponents) == 1 || (isfield(userParams, 'faultOptions') && strcmp(userParams.faultOptions, 'rakeFixed'))
                    % in this case, all entries are 1
                    L=diag(ones(1,obj.geom.N));
                else
                    if strcmp(smoothingType{2},'strike') || ( strcmp(smoothingType{2},'rake') && isfield(userParams, 'faultOptions') && strcmp(userParams.faultOptions, 'rakeCoordinates') )
                        Lpart=diag(ones(1,obj.geom.N));
                        L=[Lpart,zeros(obj.geom.N)];
                    elseif strcmp(smoothingType{2},'dip') || ( strcmp(smoothingType{2},'rakePerp') && isfield(userParams, 'faultOptions') && strcmp(userParams.faultOptions, 'rakeCoordinates') )
                        Lpart=diag(ones(1,obj.geom.N));
                        L=[zeros(obj.geom.N),Lpart];
                    elseif strcmp(smoothingType{2},'both')
                        L=diag(ones(1,2*obj.geom.N));
                        % in this case, we need a vector twice as long
                        lVector=[lVector; lVector];
                    elseif strcmp(smoothingType{2},'rake')
                        assert( ~(isfield(userParams, 'faultOptions') && strcmp(userParams.faultOptions, 'rakeCoordinates')), 'Error: unknown value smoothing arguments');
                        % not rake coordinates, so convert strike/dip
                        L=obj.Rmat(:,1:end/2)';
                    elseif strcmp(smoothingType{2},'rakePerp')
                        assert( ~(isfield(userParams, 'faultOptions') && strcmp(userParams.faultOptions, 'rakeCoordinates')), 'Error: unknown value smoothing arguments');
                        % not rake coordinates, so convert strike/dip
                        L=obj.Rmat(:,end/2+1:end)';
                    else
                        disp('smoothingType value: Slip component not recognized. Returning empty matrix.')
                        L=[];
                        lVector=[];
                    end
                end
                
            else
                disp('Smoothing type not recognized. Returning empty matrix.')
                L=[];
                lVector=[];
            end
            
        end %end function calc_smoothing_matrix
        
        % Constraint matrix computation - this function name is required
        function [K,kVector] = calc_constraint_matrix(obj, constraintType, userParams)
            disp(['Constraint type is: ' constraintType{1} '.'])
            if strcmp(constraintType{1},'positiveStress')
                % other entries in the constraint list are: depth above
                % which to apply the constraint, and the reference slip
                % rate (currently not used - for backslip it is zero)
                % load stress kernels if it was not already done
                if isempty(obj.KK)
                    [obj.KK, obj.Ksn, obj.Kdn] = unicycle_stress_kernel(obj.geom, userParams.slipComponents,userParams.stressKernelFolder);
                end
                
                % third parameter is either a dip-slip value (assumes 2D
                % dip-only constraint), or it is a file listing the
                % strike,dip values for every patch.
                
                %syntheticslip = zeros(2*obj.geom.N,1);
                
                if isnumeric(constraintType{3})
                    disp(['Warning, cannot specify dip-slip rate for stress constraint: ignoring entry (' num2str(constraintType{3}) ')'])
                    disp('Using full K matrix as stress constraint.')
                    %Vds=constraintType{3};
                    %syntheticslip = [zeros(obj.geom.N,1); ones(obj.geom.N,1).*Vds];

                    if (length(userParams.slipComponents) == 2)
%                         Kdd=obj.KK(obj.geom.N+1:end,obj.geom.N+1:end); % get the dip-dip stress component
%                         Ksd=obj.KK(obj.geom.N+1:end,1:obj.geom.N); % get the strike-dip stress component
%                         K=[0*Ksd, 0*Ksd; Ksd, Kdd]; % zero out everything except the dip component
                        kVector = zeros(length(userParams.slipComponents)*obj.geom.N,1);
                        K=obj.KK;

                    else
                        kVector = zeros(length(userParams.slipComponents)*obj.geom.N,1);
                        K=obj.KK;
                        
                    end
                    
                elseif strcmp(constraintType{3},'full')
                    disp('Using full K matrix as stress constraint.')
                    kVector = zeros(length(userParams.slipComponents)*obj.geom.N,1);
                    K=obj.KK;
                
                elseif strcmp(constraintType{3},'rake')
                    disp('Imposing stress constraint only in rake direction.')
                    K = - obj.Rmat'*obj.KK;
                    K(end/2+1:end,:)=[]; %keep only the first half, rake-parallel constraints
                   
                    kVector = zeros(obj.geom.N,1);

                else
                    error('Error: Constraint parameter #3 must be a single number (V_ds) or a file with V_ss and V_ds as separate columns for every patch.')
                end
                
                % extra options
                if isfield(userParams, 'faultOptions')
                    if strcmp(userParams.faultOptions, 'rakeFixed')
                        % rotate K matrix into rake direction
                        K = K * obj.Rmat(:,1:end/2);
                        % this option has fewer constraints:
                        kVector = zeros(obj.geom.N,1);
                    elseif strcmp(userParams.faultOptions, 'rakeCoordinates')
                        % rotate K matrix into rake direction, keeping both
                        K = K * obj.Rmat;
                    end
                end
                
                %find un-constrained patches to ignore
                depth=constraintType{2};
                Idepth=obj.geom.xc(:,3) < -depth; %note, input depth should be positive, opposite to Unicycle's convention
                
                if(strcmp(constraintType{3},'rake') || length(userParams.slipComponents) == 1)
                    Iunconstrained=Idepth;
                else
                    Iunconstrained=[Idepth; Idepth];
                end
                
                K(Iunconstrained,:)=[];
                %K(:,Iunconstrained)=0; %not needed?
                kVector(Iunconstrained)=[];
                
                % TODO: add normal stress component
                %sig0=Kdd*(-Vpl*ones(length(rcv.Vpl),1))-userParams.frictionCoeff*Kdn*(-Vpl*ones(length(rcv.Vpl),1));
                
                %TODO: fix the Unconstrained option
                
                %stressing rate constraint is based on the assumed plate slip rate. 
                %stressrate = -K * (syntheticslip); % + userParams.frictionCoeff*obj.Kdn*ones(length(obj.modelVector),1)*Vds; %TODO: use correct (variable) dip to get Vds! Also note tricky negative sign
                

                %disp(minmax(K-obj.KK))
                
            elseif strcmp(constraintType{1},'rakeSlipRate')
                % rakeSlipRate option requires that the slip rate in the
                % rake direction is not greater than the input block motion
                
                % K matrix is the projection of the strike,dip components
                % onto the rake direction
                rakeangle=obj.rakeFile(:,1);
                K = - [diag(cosd(rakeangle)), diag(sind(rakeangle))]; 
                
                % K vector is the predicted plate slip rate for each patch 
                % (stored as strike, dip components in columns 2, 3)
                kVector = obj.Vpl;

            elseif strcmp(constraintType{1},'rakePerpendicularSlipRate')
                % rakePerpendicularSlipRate option requires that the slip rate in the
                % rake-perpendicular direction is LESS than a certain value
                
                % K matrix is the projection of the strike,dip components
                % onto the rake direction
                rakeangle=obj.rakeFile(:,1);
                K = [diag(-sind(rakeangle)), diag(cosd(rakeangle))];
                
                % K vector is the predicted plate slip rate for each patch 
                % (stored as strike, dip components in columns 2, 3)
                kVector = constraintType{2} * ones(obj.geom.N,1);

            elseif strcmp(constraintType{1},'negativeRakePerpendicularSlipRate')
                % negativeRakePerpendicularSlipRate option requires that the slip rate in the
                % rake-perpendicular direction is GREATER than a certain value
                
                % K matrix is the projection of the strike,dip components
                % onto the rake direction
                rakeangle=obj.rakeFile(:,1);
                 K = - [diag(-sind(rakeangle)), diag(cosd(rakeangle))];
                
                % K vector is the predicted plate slip rate for each patch 
                % (stored as strike, dip components in columns 2, 3)
                kVector= - constraintType{2} * ones(obj.geom.N,1);
                
            else
                disp('Constraint type not recognized. Returning empty matrix.')
                K=[];
                kVector=[];
            end
        end %end function calc_smoothing_matrix
        
        function [lowerBounds,upperBounds] = calc_bounds(obj,bounds,userParams)
            if userParams.slipComponents == 1
                strikerange=1:obj.geom.N;
                diprange=[];
                
            elseif userParams.slipComponents == 2
                strikerange=[];
                diprange=1:obj.geom.N;
                
            elseif length(userParams.slipComponents) == 2
                strikerange=1:obj.geom.N;
                diprange=obj.geom.N+1:2*obj.geom.N;
                
            end
            
            lowerBounds=zeros(max(length(strikerange),length(diprange)),1);
            upperBounds=lowerBounds;
            
            for i=1:length(bounds)
                if (strcmp(bounds{i}{1},'strike'))
                    lowerBounds(strikerange)=bounds{i}{2};
                    upperBounds(strikerange)=bounds{i}{3};
                    
                elseif (strcmp(bounds{i}{1},'dip') || strcmp(bounds{i}{1},'rakePerp') )
                    lowerBounds(diprange)=bounds{i}{2};
                    upperBounds(diprange)=bounds{i}{3};
                    
                elseif strcmp(bounds{i}{1},'rakeSlipRate')
                    % Place bounds on the fault slip, based on slip rates specified in the rake file.
                    % Requires the use of 'rakeCoordinates' - if not using rake coordinates, 
                    % use the 'rakeSlipRate' constraint instead of bounds.
                    lowerBounds(strikerange)= - obj.Vpl; %note, negative slip convention
                    upperBounds(strikerange)= zeros(obj.geom.N,1); 
                    
                else
                    error('Error: Bounds type not recognized. Returning empty bounds.')
                end
            end
            
        end %end function calc_bounds
        
    end %end methods
    
end %end classdef

