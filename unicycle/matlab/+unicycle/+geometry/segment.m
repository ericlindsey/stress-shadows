classdef segment < unicycle.geometry.patch
    properties
        % name
        name;
        % index of first patches
        starti;
        % number of sub-patches
        nPatch;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = segment(x,y,z,L,W,strike,dip,rake,starti,nPatch)
            % PATCH is a class representing the geometry and physical
            % properties of fault patches, including position, orientation,
            % dimension and friction properties.
            %
            %   src = geometry.patch('path/to/source/model.flt')
            %
            % creates a instance of fault patches with a slip distribution.
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            % segment properties
            obj.starti=starti;
            obj.nPatch=nPatch;
            
            % patch properties
            obj.x=[x,y,z];
            obj.L=L;
            obj.W=W;
            obj.strike=strike;
            obj.dip=dip;
            obj.rake=rake;
            obj.N=length(x);
            
            % unit vectors in the strike direction
            obj.sv=[sind(obj.strike),cosd(obj.strike),zeros(obj.N,1)];
            
            % unit vectors in the dip direction
            obj.dv=[...
                -cosd(obj.strike).*cosd(obj.dip), ...
                +sind(obj.strike).*cosd(obj.dip), ...
                +sind(obj.dip)];
            
            % unit vectors in the normal direction
            obj.nv=[...
                +cosd(obj.strike).*sind(obj.dip), ...
                -sind(obj.strike).*sind(obj.dip), ...
                +cosd(obj.dip)];
            
            % center of fault patch
            obj.xc=[...
                obj.x(:,1)+obj.L/2.*obj.sv(:,1)-obj.W/2.*obj.dv(:,1),...
                obj.x(:,2)+obj.L/2.*obj.sv(:,2)-obj.W/2.*obj.dv(:,2),...
                obj.x(:,3)+obj.L/2.*obj.sv(:,3)-obj.W/2.*obj.dv(:,3)];
        end % constructor
    end % end methods
end