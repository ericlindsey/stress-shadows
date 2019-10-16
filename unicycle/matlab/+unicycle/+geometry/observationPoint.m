classdef observationPoint < handle
    properties
        % index
        index;
        % position
        xc;
        % solution time
        t;
        % solution vector time series
        y;
        % number of time steps
        N;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = observationPoint(index)
            % OBSERVATIONPOINT is a class representing a patch and its
            % evolution.
            %
            %   opt = geometry.observationPoint(index)
            %
            % where index is the patch coordinate on a receiver fault.
            %
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle.
            
            obj.index=index;
            obj.N=0;
            
        end % constructor
        
    end % methods
    
end
