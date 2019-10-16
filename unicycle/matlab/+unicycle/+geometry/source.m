classdef source < unicycle.geometry.patch
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = source(varargin)
            % SOURCE is a class representing the geometry and kinematics of
            % source patches.
            %
            %   src = geometry.source('path/to/source/model.flt')
            %
            % creates a instance of fault patches with a velocity 
            % distribution. 
            %
            % OUTPUT:
            %
            %   src.slip   : slip distribution (column vector)
            %   src.x      : patch upper left corner position
            %   src.L      : patch length
            %   src.W      : patch width
            %   src.strike : patch strike
            %   src.dip    : patch dip
            %   src.rake   : patch rake
            %   src.sv     : patch strike unit vector (:,3)
            %   src.dv     : patch dip unit vector (:,3)
            %   src.nv     : patch normal unit vector (:,3)
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            if nargin<1 || isempty(varargin)
                obj.N=0;
                return
            end

            filename=varargin{1};
            
            if ~iscell(filename)
                filename={filename};
            end
            obj.earthModel=varargin{2};
           
             values={};
            
             fm=[];
             seg={};
             for k=1:length(filename)
                 assert(2==exist(filename{k},'file'),['unicycle:geometry:patch:file not found ' filename{k}]);
                 [~,~,ext]=fileparts(filename{k});
                 switch ext
                     case '.seg'
                         [lfm,lseg]=obj.seg2flt(filename{k});
                         if 0~=numel(seg)
                             for i=1:length(lseg)
                                 lseg{i}.starti=lseg{i}.starti+seg{end}.starti+seg{end}.nPatch;
                             end
                         end
                         values=[values,length(seg)+(1:length(lseg))];
                         fm=[fm;lfm];
                         seg=[seg;lseg];
                         
                         switch size(fm,2)
                             case 9
                                 obj.N=size(fm,1);
                                 obj.slip=fm(:,1);
                                 obj.x=[fm(:,[3,2]),-fm(:,4)];
                                 obj.L=fm(:,5);
                                 obj.W=fm(:,6);
                                 obj.strike=fm(:,7);
                                 obj.dip=fm(:,8);
                                 obj.rake=fm(:,9);
                             otherwise
                                 error('unicycle:geometry:patch:source:incorrect file format');
                         end
                             case '.flt'
                     
                         [~,obj.slip,ys,xs,zs,obj.L,obj.W,obj.strike,obj.dip,obj.rake]=...
                             textread(filename{k},'%u %f %f %f %f %f %f %f %f %f','commentstyle','shell');
                         obj.N=length(obj.slip);
                         obj.x=[xs,ys,-zs];
                         
                     otherwise
                         error(['unicycle:geometry:patch:wrong file type ' filename{k}]);
                 end
             end
            
            
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
        
    end % methods
    
end