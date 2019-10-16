classdef Backslip_Translation_Source < Jointinv_Source
    %Backslip_Translation_Source provides an interface to a simple
    %translation source for points inside a specified region. Will add 2
    %elements to the model vector, and 2 columns to the design matrix.
    % Eric Lindsey, April 2018
    
    properties

        % Inherited from Jointinv_Source:
        %         
        % modelVector     % M x 1  - model parameters, unraveled however you like
        %
        % Specific to this class:
        poly_region       % coordinates of the polygon
        
    end
    
    methods
        
        % % % % % % % % % % % % % % % % % % %
        %                                   %
        %      c o n s t r u c t o r        %
        %                                   %
        % % % % % % % % % % % % % % % % % % %
        function obj = Backslip_Translation_Source(polygonFile, ~)
            
            % set number of parameters and preallocate model vector
            numParams=2;
            obj.modelVector = zeros( numParams, 1);
            
            %load polygon file specifying the region of the block
            obj.poly_region=load(polygonFile);
             
        end
        
        % Design matrix computation - this function name is required
        function G = calc_design_matrix(obj, dataset, userParams)
            
            % get points inside the polygon
            in_region=inpolygon(dataset.coordinates(:,1),dataset.coordinates(:,2),obj.poly_region(:,1),obj.poly_region(:,2));
            
            % add translation elements to design matrix
            % added components are the rest of the design matrix multiplied by uniform (negative) slip, then by 1/0 depending on whether the site is in the plate.
            % add a bunch of zeroes for the other data components.
            x_col=reshape([in_region,0*in_region,0*in_region]',[3*length(in_region),1]);
            y_col=reshape([0*in_region,in_region,0*in_region]',[3*length(in_region),1]);
            G = [x_col y_col]; 
            
            %keep only some of the data rows (E,N,U) in the matrix
            if isfield(userParams, 'dataComponents')
                G = keep_matrix_rows(G, userParams.dataComponents);
            end
            
        end
        
        % Smoothing matrix computation - this function name is required
        function [L,lVector] = calc_smoothing_matrix(~, ~, ~)
            
            disp('No smoothing for backslip source.')
            L=[];
            lVector=[];
            
        end %end function calc_smoothing_matrix
        
        % Constraint matrix computation - this function name is required
        function [K,kVector] = calc_constraint_matrix(~, ~, ~)
            
            %TODO
            disp('Constraints not yet implemented for backslip source.')
            K=[];
            kVector=[];
            
        end %end function calc_smoothing_matrix
        
    end
    
end
