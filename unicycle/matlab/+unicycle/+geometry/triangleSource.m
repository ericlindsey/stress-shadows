classdef triangleSource < unicycle.geometry.triangle
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = triangleSource(varargin)
            % TRIANGLESOURCE is a class representing the geometry and kinematics of
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
            %   src.x      : coordinates of vertices
            %   src.xc     : center coordinates
            %   src.rake   : patch rake
            %   src.sv     : patch strike unit vector (:,3)
            %   src.dv     : patch dip unit vector (:,3)
            %   src.nv     : patch normal unit vector (:,3)
            %
            % SEE ALSO: unicycle
              
            if nargin<1 || isempty(varargin)
                obj.N=0;
                return
            end
            
            basename=varargin{1};
            obj.earthModel=varargin{2};
            
            if ~iscell(basename)
                basename={basename};
            end
            
            obj.x=[];
            obj.vertices=[];
            obj.rake=[];
            obj.slip=[];
            for k=1:length(basename)
                fname=[basename{k} '.ned'];
                assert(2==exist(fname,'file'));
                [~,x1,x2,x3]=...
                    textread(fname,'%u %f %f %f','commentstyle','shell');
                obj.x=[obj.x;[x2,x1,-x3]];
                
                fname=[basename{k} '.tri'];
                assert(2==exist(fname,'file'));
                [~,slip,A,B,C,rake]=...
                    textread(fname,'%u %f %u %u %u %f','commentstyle','shell');
                obj.vertices=[obj.vertices;[A(:),B(:),C(:)]];
                obj.rake=[obj.rake;rake(:)];
                obj.slip=[obj.slip;slip(:)];
            end
            
            % triangle properties
            obj.N=size(obj.vertices,1);
            obj.id=1:obj.N;
            
            % center of fault patch
            obj.xc=[(obj.x(obj.vertices(:,1),1)+ ...
                     obj.x(obj.vertices(:,2),1)+ ...
                     obj.x(obj.vertices(:,3),1))/3, ...
                    (obj.x(obj.vertices(:,1),2)+ ...
                     obj.x(obj.vertices(:,2),2)+ ...
                     obj.x(obj.vertices(:,3),2))/3, ...
                    (obj.x(obj.vertices(:,1),3)+ ...
                     obj.x(obj.vertices(:,2),3)+ ...
                     obj.x(obj.vertices(:,3),3))/3];
            
            % vertices
            A=[obj.x(obj.vertices(:,1),1),obj.x(obj.vertices(:,1),2),obj.x(obj.vertices(:,1),3)];
            B=[obj.x(obj.vertices(:,2),1),obj.x(obj.vertices(:,2),2),obj.x(obj.vertices(:,2),3)];
            C=[obj.x(obj.vertices(:,3),1),obj.x(obj.vertices(:,3),2),obj.x(obj.vertices(:,3),3)];
            
            % normal vector (B-A) x (C-A)
            obj.nv=[(B(:,2)-A(:,2)).*(C(:,3)-A(:,3))-(B(:,3)-A(:,3)).*(C(:,2)-A(:,2)), ...
                    (B(:,3)-A(:,3)).*(C(:,1)-A(:,1))-(B(:,1)-A(:,1)).*(C(:,3)-A(:,3)), ...
                    (B(:,1)-A(:,1)).*(C(:,2)-A(:,2))-(B(:,2)-A(:,2)).*(C(:,1)-A(:,1))];
            norm=sqrt(obj.nv(:,1).^2+obj.nv(:,2).^2+obj.nv(:,3).^2);
            obj.nv=obj.nv./repmat(norm,1,3);
            
            % choose upward-pointing normal vectors
            obj.nv(obj.nv(:,3)<0,:)=-obj.nv(obj.nv(:,3)<0,:);
            
            % strike-direction vector
            obj.sv=[-sin(atan2(obj.nv(:,2),obj.nv(:,1))), ...
                     cos(atan2(obj.nv(:,2),obj.nv(:,1))), ...
                     zeros(obj.N,1)];
                 
            % dip-direction vector
            obj.dv=[obj.nv(:,2).*obj.sv(:,3)-obj.nv(:,3).*obj.sv(:,2), ...
                    obj.nv(:,3).*obj.sv(:,1)-obj.nv(:,1).*obj.sv(:,3), ...
                    obj.nv(:,1).*obj.sv(:,2)-obj.nv(:,2).*obj.sv(:,1)];

        end % constructor
        
    end % methods
    
end
