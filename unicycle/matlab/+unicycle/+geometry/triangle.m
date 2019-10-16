classdef triangle < handle
    properties
        % number of triangles
        N;
        % identification number
        id;
        % triangular vertices
        x;
        % Delaunay mesh information (vertices of triangles)
        vertices;
        % triangular center position
        xc;
        % slip (amplitude, rake)
        slip;
        rake;
        % unit vectors in strike, dip and normal directions
        sv;
        dv;
        nv;
        % area of triangles
        area;
        % structure of segments
        segments;
        % hash table directing segment name to segment
        segment;
        % earth model
        earthModel;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = triangle()
            % TRIANGLE is a meta class representing the geometry of fault 
            % triangular patches, defined in terms of the three vertices,
            % as illustrated below:
            %
            %                 N (x1)
            %                /
            %               /
            %   x1,x2,x3 ->A ------------- C
            %              |\             /
            %              : \           /
            %              |  \         /
            %              :   \       /
            %              |    \     /
            %              :     \   /
            %              |       B
            %              :
            %              Z (x3)
            %
            % SEE ALSO: unicycle, unicycle.geometry.triangleReceiver
            
            if (0==nargin)
                return
            end
            
        end % constructor
        
        
        function [varargout]=tractionKernels(obj,rcv)
            % TRACTIONKERNELS computes the traction on receiver faults due
            % to motion of dislocations. 
            %
            % rcv - receiver fault
            
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.tractionKernels(obj,rcv);
        end
                
        function [varargout]=stressKernels(obj,rcv)
            % STRESSKERNELS computes the stress on receiver faults due to 
            % motion of dislocations. 
            %
            % rcv - receiver fault
            
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.stressKernels(obj,rcv);
        end
        
        function [varargout]=displacementKernels(obj,x,vecsize)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of dislocations.
            %
            % x       - observations points coordinates
            % vecsize - length of displacement vector
       
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.displacementKernels(obj,x,vecsize);
        end
        
        function plotPatch(obj,varargin)
            % PLOTPATCH plot contour of slip patches in 3d
            %
            % SEE ALSO: unicycle
            if 1>obj.N
                fprintf('unicycle.geometry.patch: nothing to plot\n');
                return
            end
            if nargin > 1
                trisurf(obj.vertices,obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                    varargin{1}), shading flat;
            else
                hold on
                trisurf(obj.vertices,obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                    0.*obj.x(:,3),'FaceColor','None','LineWidth',1);
            end
        end
        
        function [xp,yp,zp,dim]=computeVertexPosition(obj)
            % COMPUTEVERTEXPOSITION computes the position of vertices
            %
            % SEE ALSO: unicycle
            
            xp=[obj.x(obj.vertices(:,1),1),obj.x(obj.vertices(:,2),1),obj.x(obj.vertices(:,3),1)]';
            yp=[obj.x(obj.vertices(:,1),2),obj.x(obj.vertices(:,2),2),obj.x(obj.vertices(:,3),2)]';
            zp=[obj.x(obj.vertices(:,1),3),obj.x(obj.vertices(:,2),3),obj.x(obj.vertices(:,3),3)]';
            dim=3;
        end
        
        function plotUnitVectors(obj,sc)
            % PLOTUNITVECTORS plot normal, strike-slip and dip-slip unit
            % vectors.
            %
            % SEE ALSO: unicycle
            
            hold on
            trisurf(obj.vertices,obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                0.*obj.x(:,3),'FaceColor','None','LineWidth',1);
            
            plot3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),'+');
            
            plot3(obj.x(obj.vertices(:,1),1),obj.x(obj.vertices(:,1),2),obj.x(obj.vertices(:,1),3),'+');
            plot3(obj.x(obj.vertices(:,2),1),obj.x(obj.vertices(:,2),2),obj.x(obj.vertices(:,2),3),'+');
            plot3(obj.x(obj.vertices(:,3),1),obj.x(obj.vertices(:,3),2),obj.x(obj.vertices(:,3),3),'+');
            
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),sc*obj.nv(:,1),sc*obj.nv(:,2),sc*obj.nv(:,3),0,'r');
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),sc*obj.sv(:,1),sc*obj.sv(:,2),sc*obj.sv(:,3),0,'g');
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),sc*obj.dv(:,1),sc*obj.dv(:,2),sc*obj.dv(:,3),0,'b');
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %              plot patch index in 3d              %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function plotPatchIndex(obj,varargin)
            if 1>obj.N
                fprintf('unicycle.geometry.triangle: nothing to plot\n');
                return
            end
            trisurf(obj.vertices,obj.x(:,1),obj.x(:,2),obj.x(:,3), ...
                    0.*obj.x(:,3),'EdgeAlpha',0.2,'FaceColor','None','LineWidth',1);
            for k=1:length(obj.vertices)
                text(obj.xc(k,1),obj.xc(k,2),obj.xc(k,3),num2str(k));
            end
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %            plot slip vectors in 3d               %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function plotSlipVectors(obj,ss,ds,scale,varargin)
            if 1>obj.N
                return
            end
            if nargin > 4
                s=repmat(ss,[1,3]).*obj.sv+repmat(ds,[1,3]).*obj.dv;
                quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),...
                    scale*s(:,1),scale*s(:,2),scale*s(:,3),0,varargin{1});
            else
                s=repmat(ss,[1,3]).*obj.sv+repmat(ds,[1,3]).*obj.dv;
                quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),...
                    scale*s(:,1),scale*s(:,2),scale*s(:,3),0,'k');
            end
        end
        
    end % methods
    
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                    %
        %  compute normal vectors from position of vertices  %
        %                                                    %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [nv,sv,dv,area] = computeUnitVectors(x,vertices)
            % COMPUTEUNITVECTORS compute the triangle unit normal, strike
            % and dip direction vectors.
            %
            %   [nv,sv,dv]=triangle.computeUnitVectors(x,vertices)
            %
            % INPUT:
            %   x        - position of mesh points
            %   vertices - list of vertices forming triangles
            %
            % SEE ALSO: unicycle, unicycle.geometry.triangle
            
            % number of points
            nVertices=size(vertices,1);
            
            % vertices
            A=[x(vertices(:,1),1),x(vertices(:,1),2),x(vertices(:,1),3)];
            B=[x(vertices(:,2),1),x(vertices(:,2),2),x(vertices(:,2),3)];
            C=[x(vertices(:,3),1),x(vertices(:,3),2),x(vertices(:,3),3)];
            
            % normal vector (B-A) x (C-A)
            nv=[(B(:,2)-A(:,2)).*(C(:,3)-A(:,3))-(B(:,3)-A(:,3)).*(C(:,2)-A(:,2)), ...
                (B(:,3)-A(:,3)).*(C(:,1)-A(:,1))-(B(:,1)-A(:,1)).*(C(:,3)-A(:,3)), ...
                (B(:,1)-A(:,1)).*(C(:,2)-A(:,2))-(B(:,2)-A(:,2)).*(C(:,1)-A(:,1))];
            area=sqrt(nv(:,1).^2+nv(:,2).^2+nv(:,3).^2)/2;
            nv=nv./repmat(2*area,1,3);
            
            % choose upward-pointing normal vectors
            nv(nv(:,3)<0,:)=-nv(nv(:,3)<0,:);
            
            % strike-direction vector
            sv=[-sin(atan2(nv(:,2),nv(:,1))), ...
                 cos(atan2(nv(:,2),nv(:,1))), ...
                 zeros(nVertices,1)];
                 
            % dip-direction vector
            dv=[nv(:,2).*sv(:,3)-nv(:,3).*sv(:,2), ...
                nv(:,3).*sv(:,1)-nv(:,1).*sv(:,3), ...
                nv(:,1).*sv(:,2)-nv(:,2).*sv(:,1)];
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %        load the triangle mesh from file          %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [Vpl,tri,rak] = loadTri(filename)
            % LOADTRI loads the .tri mesh file and detect the availability
            % of Vpl.
            %
            %   [Vpl,tri,rake]=loadTri(filename)
            %
            % where tri is a list Vpl, mesh vertices and rake parameters.
            
            % open and count the number of columns
            fid=fopen(filename);
            line=strtrim(fgetl(fid));             
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % Open the file, get the first line
            end
            fclose(fid);
            nColumn=numel(strsplit(line));
            
            switch nColumn
                case 6 % if Vpl is present
                    [~,Vpl,i1,i2,i3,rak]=...
                        textread(filename,'%u %f %d %d %d %f',...
                        'commentstyle','shell');
                    tri = [i1(:),i2(:),i3(:)];
                case 5 % no Vpl
                    [~,i1,i2,i3,rak]=...
                        textread(filename,'%u %d %d %d %f',...
                        'commentstyle','shell');
                    tri = [i1(:),i2(:),i3(:)];
                    Vpl = zeros(size(i1));
                otherwise
                    error('unicycle:geometry:triangle:invalid file format');
            end
        end
    end % methods (Static)

end
