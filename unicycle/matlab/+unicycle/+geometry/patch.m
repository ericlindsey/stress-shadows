classdef patch < handle
    properties
        N;
        % identification number
        id;
        % patch position (upper corner and centre)
        x;xc;
        % patch dimension (length and width)
        L;W;
        % patch orientation (degrees)
        strike;dip;
        % patch slip (amplitude, rake)
        slip;rake;
        % unit vectors in strike, dip, normal and rake directions
        sv;dv;nv;
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
        function obj = patch()
            % PATCH is a meta class representing the geometry of fault 
            % patches, defined in terms of position, orientation,
            % and slip, as illustrated below:
            %
            %                 N (x1)
            %                /
            %               /| Strike
            %   x1,x2,x3 ->@------------------------      (x2)
            %              |\        p .            \ W
            %              :-\      i .              \ i
            %              |  \    l .                \ d
            %              :90 \  S .                  \ t
            %              |-Dip\  .                    \ h
            %              :     \. | Rake               \
            %              |      -------------------------
            %              :             L e n g t h
            %              Z (x3)
            %
            % SEE ALSO: unicycle, unicycle.geometry.triangle
            
            if (0==nargin)
                return
            end
            
        end % constructor
        
        function [varargout]=tractionKernels(obj,rcv)
            % STRESSKERNELS computes the traction on receiver faults due to
            % motion of dislocations.
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
        
        function [tri]=toTriangle(obj,varargin)
            % TOTRIANGLE provides a triangle class representation of the
            % fault plane with two triangles per rectangle patch.
            %
            %   [triangle] = obj.toTriangle()
            %
            % OUTPUT:
            %   triangle - an object of class triangle.
            
            import unicycle.geometry.triangle;
            
            if isempty(varargin)
                tri=triangle();
            else
                tri=varargin{1};
            end
            
            % number of triangles
            tri.N=2*obj.N;
            
            % identification number
            tri.id=1:tri.N;
            
            % triangular vertices
            tri.x=[ ...
                obj.x; ...
                obj.x+repmat(obj.L,1,3).*obj.sv; ...
                obj.x+repmat(obj.L,1,3).*obj.sv-repmat(obj.W,1,3).*obj.dv; ...
                obj.x-repmat(obj.W,1,3).*obj.dv];
            
            % mesh information (vertices of triangles)
            position=(1:obj.N)';
            tri.vertices=[ ...
                [position,obj.N+position,2*obj.N+position];
                [position,2*obj.N+position,3*obj.N+position]];
            
            % triangular center position
            tri.xc=[(tri.x(tri.vertices(:,1),1)+ ...
                     tri.x(tri.vertices(:,2),1)+ ...
                     tri.x(tri.vertices(:,3),1))/3, ...
                    (tri.x(tri.vertices(:,1),2)+ ...
                     tri.x(tri.vertices(:,2),2)+ ...
                     tri.x(tri.vertices(:,3),2))/3, ...
                    (tri.x(tri.vertices(:,1),3)+ ...
                     tri.x(tri.vertices(:,2),3)+ ...
                     tri.x(tri.vertices(:,3),3))/3];
            
            % unit vectors
            [tri.nv,tri.sv,tri.dv]=triangle.computeUnitVectors(tri.x,tri.vertices);
            
            % slip (amplitude, rake)
            tri.slip=[obj.slip;obj.slip];
            tri.rake=[obj.rake;obj.rake];
        end
        
        function [xp,yp,zp,dim]=computeVertexPosition(obj)
            % COMPUTEVERTEXPOSITION computes the position of vertices
            %
            % SEE ALSO: unicycle
            
            [xp,yp,zp,~]=unicycle.geometry.transform4patch_general(...
                    obj.x(:,1),obj.x(:,2),-obj.x(:,3),0*obj.x(:,3),...
                    obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
            dim=4;
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %       plot contour of slip patches in 3d         %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function plotPatch(obj,varargin)
            if 1>obj.N
                fprintf('unicycle.geometry.patch: nothing to plot\n');
                return
            end
            if nargin > 1
                [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                    obj.x(:,1),obj.x(:,2),-obj.x(:,3),varargin{1},...
                    obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
                patch(xp,yp,-zp,up);
            else
                [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                    obj.x(:,1),obj.x(:,2),-obj.x(:,3),obj.L(:)*0,...
                    obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
                patch(xp,yp,-zp,0*up,'FaceColor','None','LineWidth',1);
                
            end
        end
        
        function plotUnitVectors(obj,sc)
            % PLOTUNITVECTORS plot normal, strike-slip and dip-slip unit
            % vectors.
            %
            % SEE ALSO: unicycle
            
            [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                obj.x(:,1),obj.x(:,2),-obj.x(:,3),obj.L(:)*0,...
                obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
            patch(xp,yp,-zp,0*up,'FaceColor','None','LineWidth',1);
            
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
                fprintf('unicycle.geometry.patch: nothing to plot\n');
                return
            end
            [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                obj.x(:,1),obj.x(:,2),-obj.x(:,3),obj.L(:)*0,...
                obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
            patch(xp,yp,-zp,0*up,'FaceColor','None','LineWidth',1);
            for k=1:length(obj.L)
                text(obj.x(k,1),obj.x(k,2),obj.x(k,3),num2str(k));
            end
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                          %
        %        export geometry to Paraview VTP file              %
        %                                                          %
        % INPUT:                                                   %
        %                                                          % 
        % scale        - scaling factor of geometrical features    %
        % fname        - output file name                          %
        % varargin     - pairs of variable name (string) and       %
        %                attributes.                               %
        %                                                          %
        % EXAMPLE:                                                 %
        %                                                          %
        %   flt.exportVTP(1e-3,'output/flt.vtp','slip',flt.slip)   %
        %                                                          %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportVTP(obj,scale,fname,varargin)
            if 1>obj.N
                return
            end
            
            [xp,yp,zp,dim]=obj.rcv.computeVertexPosition();
            
            if nargin > 1
                unicycle.export.exportvtk_rfaults( ...
                    scale*xp, ...
                    scale*yp, ...
                    scale*zp, ...
                    dim, ...
                    fname,...
                    varargin);
            else
                unicycle.export.exportvtk_rfaults( ...
                    scale*xp, ...
                    scale*yp, ...
                    scale*zp, ...
                    dim, ...
                    fname);
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
            nargin
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
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [varargout] = seg2flt(filename)
            % SEG2FLT converts a list of segments from a file to a list
            % of patches
            %
            %   flt=seg2flt(filename)
            %
            % or
            %
            %   [flt,seg]=seg2flt(filename)
            %
            % where flt is a list of [N,nColumn] patch geometry parameters, and
            % seg is a cell array of class 'geometry.segment'.
            
            % The length of nColumn will tell how many arguments the file
            % has nColumn=14 implies there is Vpl, else nColumn = 13.
            fid=fopen(filename);  % open and count the number of columns
            line=strtrim(fgetl(fid));             
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % line ignores all shell comments and check for the first line that contains data
            end
            fclose(fid);
            nColumn=numel(strsplit(line)); 
            
            switch nColumn
                case 14 % if there is Vpl
                    [~,Vpl,x1,x2,x3,len,width,str,d,rak,lo,wo,alphal,alphaw]=...
                        textread(filename,'%u %f %f %f %f %f %f %f %f %f %f %f %f %f',...
                        'commentstyle','shell');
                case 13 % no Vpl
                    [~,x1,x2,x3,len,width,str,d,rak,lo,wo,alphal,alphaw]=...
                        textread(filename,'%u %f %f %f %f %f %f %f %f %f %f %f %f',...
                        'commentstyle','shell');
                otherwise
                    error('unicycle:geometry:patch:invalid file format');
            end
                      
            fm=[];
            seg=cell(length(x1),1);
            for k=1:length(x1)
                % list of patches for current segment
                flt=unicycle.geometry.flt2flt([x1(k);x2(k);x3(k)],len(k),width(k),...
                    str(k),d(k),rak(k),...
                    lo(k),wo(k),alphal(k),alphaw(k));
                                   
                switch nColumn
                    case 13
                    case 14
                        flt=[Vpl(k)*ones(size(flt,1),1),flt];
                    otherwise
                        error('unicycle:geometry:patch:coding error');
                end
                % segment
                seg{k}=unicycle.geometry.segment(x2(k),x1(k),-x3(k),len(k),width(k),str(k),d(k),rak(k),size(fm,1),size(flt,1));
                % list of patches for all segments
                fm=[fm;flt];
            end
            
            if 1==nargout
                varargout{1}=fm;
            elseif 2==nargout
                varargout{1}=fm;
                varargout{2}=seg;
            end
        end
    end % methods (Static)
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [varargout] = loadFlt(filename)
            % LOADFLT loads the file and creates flt based on  
            % availability of Vpl
            %
            %   flt=loadFlt(filename)
            %
            %
            % where flt is a list of [N,n] patch geometry parameters
            
            fid=fopen(filename);  % open and count the number of columns
            line=strtrim(fgetl(fid));             
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % Open the file, get the first line
            end
            fclose(fid);
            nColumn=numel(strsplit(line));
            
            switch nColumn
                case 10 % if there is Vpl
                    [~,Vpl,x1,x2,x3,len,width,str,d,rak]=...
                        textread(filename,'%u %f %f %f %f %f %f %f %f %f',...
                        'commentstyle','shell');
                    fm = [Vpl,x1,x2,x3,len,width,str,d,rak];
                case 9 % no Vpl
                    [~,x1,x2,x3,len,width,str,d,rak]=...
                        textread(filename,'%u %f %f %f %f %f %f %f %f',...
                        'commentstyle','shell');
                    fm = [x1,x2,x3,len,width,str,d,rak];
                otherwise
                    error('unicycle:geometry:patch:invalid file format');
            end
            
            varargout{1}=fm;
           
            
        end
    end % methods (Static)
end