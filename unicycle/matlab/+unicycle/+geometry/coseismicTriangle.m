classdef coseismicTriangle < unicycle.geometry.triangleSource
    properties
        % time to event
        t0;
        % name of event
        name;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = coseismicTriangle(varargin)
            % PATCH is a class representing the geometry and physical
            % properties of fault patches, including position, orientation,
            % dimension and friction properties.
            %
            %   src = geometry.coseismicTriangle('path/to/source/model',t0)
            %
            % creates a instance of fault patches with a slip distribution
            % defined in 'path/to/source/model.{ned,tri}' to occur at time t0.
            %
            % SEE ALSO: unicycle
            
            obj@unicycle.geometry.triangleSource(varargin{1},varargin{3});
            
            if isempty(varargin)
                return
            end
            
            obj.name=varargin{1};
            obj.t0=varargin{2};

        end % constructor
        
    end % methods
   
end