classdef shearZone < handle
    properties
        N;
        % identification number
        id;
        % shear zone position (upper corner at mid thickness and centre)
        x;xc;
        % shear zone dimension (length and width)
        L;T;W;
        % shear zone orientation (degrees)
        strike;dip;
        % shear zone strain
        eps;
        % unit vectors in strike, dip, normal and rake directions
        sv;dv;nv;
        % structure of levels
        levels;
        % hash table directing level name to level index
        levelHashTable;
        % earth model
        earthModel;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = shearZone()
            % SHEARZONE is a meta class representing the geometry of shear 
            % zones, defined in terms of position, orientation,
            % and strain, as illustrated below:
            %
            %                 N (x1)
            %                /
            %               /| Strike
            %   x1,x2,x3 ->@------------------------      (x2)
            %              |\                       \ W
            %              :-\                       \ i       +
            %              |  \                       \ d     / s
            %              :90 \                       \ t   / s
            %              |-Dip\                       \ h / e
            %              :     \                       \ / n
            %              |      ------------------------+ k   (T)
            %              :             L e n g t h     / c
            %              Z (x3)            (L)        / i
            %                                          / h
            %                                         / t
            %                                        +
            %
            % SEE ALSO: unicycle, unicycle.geometry.triangle
            
            if 0==nargin
                obj.N=0;
                return
            end
            
        end % constructor
        
        function [varargout]=tractionKernels(obj,rcv)
            % TRACTIONKERNELS computes the traction on receiver faults due
            % to strain on shear zones.
            %
            % rcv - receiver fault
            
            varargout = cell(1,nargout);
            [varargout{:}]=obj.earthModel.tractionKernels(obj,rcv);
        end
        
        function [varargout]=stressKernels(obj,rcv)
            % STRESSKERNELS computes the stress on receiver shear zones due
            % to motion of dislocations.
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
        
        function [xp,yp,zp,dim]=computeVertexPosition(o)
            % COMPUTEVERTEXPOSITION computes the position of vertices
            %
            % SEE ALSO: unicycle
            
            xyz=cell(3,1);
            for k=1:3
                xyz{k}=[ ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    ]';
            end
            [xp,yp,zp]=deal(xyz{:});
            dim=24;
            
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %       plot contour of slip patches in 3d         %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function plotShearZoneEdges(obj,varargin)
            % PLOTSHEARZONESEDGES plot the 12 edges of shear zones as
            % follows:
            %
            %          +-------+
            %         /       /|
            %        /       / |
            %       +-------+  +
            %       |       | /
            %       |       |/
            %       +-------+
            %
            % SEE ALSO: unicycle
            
            if 1>obj.N
                fprintf('unicycle.geometry.shearzone: nothing to plot\n');
                return
            end
            
            if 1>=nargin
                
                % 12 segments
                plot3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),'k+');
                
                pm=[1,-1];
                
                % 4 along-strike segments
                for k=1:2
                    for j=1:2
                        x1=[obj.xc(:,1)-obj.sv(:,1).*obj.L/2+pm(j)*obj.dv(:,1).*obj.W/2+pm(k)*obj.nv(:,1).*obj.T/2, ...
                            obj.xc(:,1)+obj.sv(:,1).*obj.L/2+pm(j)*obj.dv(:,1).*obj.W/2+pm(k)*obj.nv(:,1).*obj.T/2, ...
                            NaN(size(obj.sv(:,1)))]';
                        x2=[obj.xc(:,2)-obj.sv(:,2).*obj.L/2+pm(j)*obj.dv(:,2).*obj.W/2+pm(k)*obj.nv(:,2).*obj.T/2, ...
                            obj.xc(:,2)+obj.sv(:,2).*obj.L/2+pm(j)*obj.dv(:,2).*obj.W/2+pm(k)*obj.nv(:,2).*obj.T/2, ...
                            NaN(size(obj.sv(:,2)))]';
                        x3=[obj.xc(:,3)-obj.sv(:,3).*obj.L/2+pm(j)*obj.dv(:,3).*obj.W/2+pm(k)*obj.nv(:,3).*obj.T/2, ...
                            obj.xc(:,3)+obj.sv(:,3).*obj.L/2+pm(j)*obj.dv(:,3).*obj.W/2+pm(k)*obj.nv(:,3).*obj.T/2, ...
                            NaN(size(obj.sv(:,3)))]';
                        
                        plot3(x1(:),x2(:),x3(:),'k-');
                    end
                end
                
                % 4 down-dip segments
                for k=1:2
                    for j=1:2
                        x1=[obj.xc(:,1)+pm(j)*obj.sv(:,1).*obj.L/2+obj.dv(:,1).*obj.W/2+pm(k)*obj.nv(:,1).*obj.T/2, ...
                            obj.xc(:,1)+pm(j)*obj.sv(:,1).*obj.L/2-obj.dv(:,1).*obj.W/2+pm(k)*obj.nv(:,1).*obj.T/2, ...
                            NaN(size(obj.sv(:,1)))]';
                        x2=[obj.xc(:,2)+pm(j)*obj.sv(:,2).*obj.L/2+obj.dv(:,2).*obj.W/2+pm(k)*obj.nv(:,2).*obj.T/2, ...
                            obj.xc(:,2)+pm(j)*obj.sv(:,2).*obj.L/2-obj.dv(:,2).*obj.W/2+pm(k)*obj.nv(:,2).*obj.T/2, ...
                            NaN(size(obj.sv(:,2)))]';
                        x3=[obj.xc(:,3)+pm(j)*obj.sv(:,3).*obj.L/2+obj.dv(:,3).*obj.W/2+pm(k)*obj.nv(:,3).*obj.T/2, ...
                            obj.xc(:,3)+pm(j)*obj.sv(:,3).*obj.L/2-obj.dv(:,3).*obj.W/2+pm(k)*obj.nv(:,3).*obj.T/2, ...
                            NaN(size(obj.sv(:,3)))]';
                        
                        plot3(x1(:),x2(:),x3(:),'k-');
                    end
                end
                
                % 4 normal-direction segments
                for k=1:2
                    for j=1:2
                        x1=[obj.xc(:,1)+pm(j)*obj.sv(:,1).*obj.L/2+pm(k)*obj.dv(:,1).*obj.W/2+obj.nv(:,1).*obj.T/2, ...
                            obj.xc(:,1)+pm(j)*obj.sv(:,1).*obj.L/2+pm(k)*obj.dv(:,1).*obj.W/2-obj.nv(:,1).*obj.T/2, ...
                            NaN(size(obj.sv(:,1)))]';
                        x2=[obj.xc(:,2)+pm(j)*obj.sv(:,2).*obj.L/2+pm(k)*obj.dv(:,2).*obj.W/2+obj.nv(:,2).*obj.T/2, ...
                            obj.xc(:,2)+pm(j)*obj.sv(:,2).*obj.L/2+pm(k)*obj.dv(:,2).*obj.W/2-obj.nv(:,2).*obj.T/2, ...
                            NaN(size(obj.sv(:,2)))]';
                        x3=[obj.xc(:,3)+pm(j)*obj.sv(:,3).*obj.L/2+pm(k)*obj.dv(:,3).*obj.W/2+obj.nv(:,3).*obj.T/2, ...
                            obj.xc(:,3)+pm(j)*obj.sv(:,3).*obj.L/2+pm(k)*obj.dv(:,3).*obj.W/2-obj.nv(:,3).*obj.T/2, ...
                            NaN(size(obj.sv(:,3)))]';
                        
                        plot3(x1(:),x2(:),x3(:),'k-');
                    end
                end
                
            else
                % vertical slice in strike direction
                [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                    obj.x(:,1),obj.x(:,2),-obj.x(:,3),varargin{1},...
                    obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
                patch(xp,yp,-zp,up,'EdgeColor','None');
                
                % vertical slice perpendicular to strike direction
                [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                    obj.x(:,1)+obj.sv(:,1).*obj.L/2-obj.nv(:,1).*obj.T/2, ...
                    obj.x(:,2)+obj.sv(:,2).*obj.L/2-obj.nv(:,2).*obj.T/2, ...
                  -(obj.x(:,3)+obj.sv(:,3).*obj.L/2-obj.nv(:,3).*obj.T/2), ...
                    varargin{1},...
                    obj.L(:),obj.W(:),obj.dip(:),obj.strike(:)+90);
                patch(xp,yp,-zp,up,'EdgeColor','None');
                
                % horizontal slice
                [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                    obj.x(:,1)-obj.dv(:,1).*obj.W/2-obj.nv(:,1).*obj.T/2, ...
                    obj.x(:,2)-obj.dv(:,2).*obj.W/2-obj.nv(:,2).*obj.T/2, ...
                  -(obj.x(:,3)-obj.dv(:,3).*obj.W/2-obj.nv(:,3).*obj.T/2), ...
                    varargin{1},...
                    obj.L(:),obj.T(:),obj.dip(:)-90,obj.strike(:));
                patch(xp,yp,-zp,up,'EdgeColor','None');
            end
            
        end
        
        function plotShearZoneById(obj,shznum)
            % PLOTSHEARZONES plot the 12 edges of shear zones based on the
            % its number:
            % 
            %
            %          +-------+
            %         /       /|
            %        /       / |
            %       +-------+  +
            %       |       | /
            %       |       |/
            %       +-------+
            %
            % EXAMPLE:                                                   
            %   shz.plotShearZoneById(2) or shz.plotShearZoneById([1:3])
            
            if 1>obj.N
                fprintf('unicycle.geometry.shearzone: nothing to plot\n');
                return
            end
            
            % 12 segments
            plot3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),'k+');

            pm=[1,-1];

            % 4 along-strike segments
            for k=1:2
                for j=1:2
                    x1=[obj.xc(shznum,1)-obj.sv(shznum,1).*obj.L(shznum)/2+pm(j)*obj.dv(shznum,1).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,1).*obj.T(shznum)/2, ...
                        obj.xc(shznum,1)+obj.sv(shznum,1).*obj.L(shznum)/2+pm(j)*obj.dv(shznum,1).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,1).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,1)))]';
                    x2=[obj.xc(shznum,2)-obj.sv(shznum,2).*obj.L(shznum)/2+pm(j)*obj.dv(shznum,2).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,2).*obj.T(shznum)/2, ...
                        obj.xc(shznum,2)+obj.sv(shznum,2).*obj.L(shznum)/2+pm(j)*obj.dv(shznum,2).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,2).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,2)))]';
                    x3=[obj.xc(shznum,3)-obj.sv(shznum,3).*obj.L(shznum)/2+pm(j)*obj.dv(shznum,3).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,3).*obj.T(shznum)/2, ...
                        obj.xc(shznum,3)+obj.sv(shznum,3).*obj.L(shznum)/2+pm(j)*obj.dv(shznum,3).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,3).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,3)))]';

                    plot3(x1(:),x2(:),x3(:),'k-');
                end
            end

            % 4 down-dip segments
            for k=1:2
                for j=1:2
                    x1=[obj.xc(shznum,1)+pm(j)*obj.sv(shznum,1).*obj.L(shznum)/2+obj.dv(shznum,1).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,1).*obj.T(shznum)/2, ...
                        obj.xc(shznum,1)+pm(j)*obj.sv(shznum,1).*obj.L(shznum)/2-obj.dv(shznum,1).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,1).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,1)))]';
                    x2=[obj.xc(shznum,2)+pm(j)*obj.sv(shznum,2).*obj.L(shznum)/2+obj.dv(shznum,2).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,2).*obj.T(shznum)/2, ...
                        obj.xc(shznum,2)+pm(j)*obj.sv(shznum,2).*obj.L(shznum)/2-obj.dv(shznum,2).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,2).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,2)))]';
                    x3=[obj.xc(shznum,3)+pm(j)*obj.sv(shznum,3).*obj.L(shznum)/2+obj.dv(shznum,3).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,3).*obj.T(shznum)/2, ...
                        obj.xc(shznum,3)+pm(j)*obj.sv(shznum,3).*obj.L(shznum)/2-obj.dv(shznum,3).*obj.W(shznum)/2+pm(k)*obj.nv(shznum,3).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,3)))]';

                    plot3(x1(:),x2(:),x3(:),'k-');
                end
            end

            % 4 normal-direction segments
            for k=1:2
                for j=1:2
                    x1=[obj.xc(shznum,1)+pm(j)*obj.sv(shznum,1).*obj.L(shznum)/2+pm(k)*obj.dv(shznum,1).*obj.W(shznum)/2+obj.nv(shznum,1).*obj.T(shznum)/2, ...
                        obj.xc(shznum,1)+pm(j)*obj.sv(shznum,1).*obj.L(shznum)/2+pm(k)*obj.dv(shznum,1).*obj.W(shznum)/2-obj.nv(shznum,1).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,1)))]';
                    x2=[obj.xc(shznum,2)+pm(j)*obj.sv(shznum,2).*obj.L(shznum)/2+pm(k)*obj.dv(shznum,2).*obj.W(shznum)/2+obj.nv(shznum,2).*obj.T(shznum)/2, ...
                        obj.xc(shznum,2)+pm(j)*obj.sv(shznum,2).*obj.L(shznum)/2+pm(k)*obj.dv(shznum,2).*obj.W(shznum)/2-obj.nv(shznum,2).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,2)))]';
                    x3=[obj.xc(shznum,3)+pm(j)*obj.sv(shznum,3).*obj.L(shznum)/2+pm(k)*obj.dv(shznum,3).*obj.W(shznum)/2+obj.nv(shznum,3).*obj.T(shznum)/2, ...
                        obj.xc(shznum,3)+pm(j)*obj.sv(shznum,3).*obj.L(shznum)/2+pm(k)*obj.dv(shznum,3).*obj.W(shznum)/2-obj.nv(shznum,3).*obj.T(shznum)/2, ...
                        NaN(size(obj.sv(shznum,3)))]';

                    plot3(x1(:),x2(:),x3(:),'k-');
                end
            end
%             for k=1:length(shznum)
%                 text(obj.xc(shznum(k),1),obj.xc(shznum(k),2),obj.xc(shznum(k),3),num2str(shznum(k)));
%             end
        end
        
        function plotUnitVectors(obj,sc)
            % PLOTUNITVECTORS plot normal, strike-slip and dip-slip unit
            % vectors.
            %
            % SEE ALSO: unicycle
            
            [xp,yp,zp,up]=unicycle.geometry.transform4patch_general(...
                obj.x(:,1),obj.x(:,2),-obj.x(:,3),obj.L(:)*0,...
                obj.L(:),obj.W(:),obj.dip(:),obj.strike(:));
            %patch(xp,yp,-zp,0*up,'FaceColor','None','LineWidth',1);
            
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),sc*obj.nv(:,1),sc*obj.nv(:,2),sc*obj.nv(:,3),0,'r');
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),sc*obj.sv(:,1),sc*obj.sv(:,2),sc*obj.sv(:,3),0,'g');
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),sc*obj.dv(:,1),sc*obj.dv(:,2),sc*obj.dv(:,3),0,'b');
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %            plot shear zone index in 3d           %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function plotShearZoneIndex(obj,varargin)
            if 1>obj.N
                fprintf('unicycle.geometry.shearZone: nothing to plot\n');
                return
            end
            
            for k=1:length(obj.L)
                text(obj.xc(k,1),obj.xc(k,2),obj.xc(k,3),num2str(k));
            end
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %            plot slip vectors in 3d               %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function plotSlipVectors(obj,ss,ds,scale)
            if 1>obj.N
                return
            end
            
            s=repmat(ss,[1,3]).*obj.sv+repmat(ds,[1,3]).*obj.dv;
            quiver3(obj.xc(:,1),obj.xc(:,2),obj.xc(:,3),...
                scale*s(:,1),scale*s(:,2),scale*s(:,3),0);    
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                            %
        %        export geometry to Relax .shz format                %
        %                                                            %
        % INPUT:                                                     %
        %                                                            % 
        % scale        - scaling factor of geometrical features      %
        % fname        - output file name                            %
        % varargin     - pairs of variable name (string) and         %
        %                attributes.                                 %
        %                                                            %
        % EXAMPLE:                                                   %
        %                                                            %
        %   shz.exportVTP(1e-3,'output/shearzone.shz','strain',flt.slip)   %
        %                                                            %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportSHZ(obj,fname,varargin)
            if 1>obj.N
                return
            end
            
            fid=fopen(fname,'wt');
            fprintf(fid,'# export from unicycle\n');
            fprintf(fid,'# n e11 e12 e13 e22 e23 e33 x1 x2 x3 length width thickness strike dip\n');
            if nargin > 5
                fprintf(fid,'%05.5d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', ...
                    [cumsum(ones(size(obj.xc,1),1)) varargin{1} varargin{2} varargin{3} varargin{4} varargin{5} varargin{6} obj.x(:,2) obj.x(:,1) -obj.x(:,3) obj.L obj.W obj.T obj.strike obj.dip]');
            else
                fprintf(fid,'%05.5d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', ...
                    [cumsum(ones(size(obj.xc,1),1)) zeros(size(obj.xc,1),6) obj.x(:,2) obj.x(:,1) -obj.x(:,3) obj.L obj.W obj.T obj.strike obj.dip]');
            end
            fclose(fid);
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                            %
        %        export geometry to Paraview VTP file                %
        %                                                            %
        % INPUT:                                                     %
        %                                                            % 
        % scale        - scaling factor of geometrical features      %
        % fname        - output file name                            %
        % varargin     - pairs of variable name (string) and         %
        %                attributes.                                 %
        %                                                            %
        % EXAMPLE:                                                   %
        %                                                            %
        %   shz.exportVTP(1e-3,'output/shz.vtp','strain',flt.slip)   %
        %                                                            %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportVTP(obj,scale,fname,varargin)
            if 1>obj.N
                return
            end
            
            [xp,yp,zp,dim]=obj.computeVertexPosition();
            
            if nargin > 1
                unicycle.export.exportVTKshearZone( ...
                    scale*xp, ...
                    scale*yp, ...
                    scale*zp, ...
                    fname,...
                    varargin{:});
            else
                unicycle.export.exportVTKshearZone( ...
                    scale*xp, ...
                    scale*yp, ...
                    scale*zp, ...
                    fname);
            end
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                            %
        %        export geometry to GMT .xyz format                  %
        %                                                            %
        % INPUT:                                                     %
        %                                                            % 
        % sideIndex    - side index                                  %
        %                  1:top, 2:bottom,                          % 
        %                  3:left, 4:right,                          % 
        %                  5:front, 6:back                           %
        % scale        - scaling factor of geometrical features      %
        % fname        - output file name                            %
        % value        - attribute to plot                           % 
        %                                                            %
        % EXAMPLE:                                                   %
        %                                                            %
        %   shz.exportSideXYZ(1,1e-3,'output/shz.xyz',value)         %
        %                                                            %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportSideXYZ(o,sideIndex,scale,fname,field)
            if 1>o.N
                return
            end
            
            xyz=cell(3,1);
            for k=1:3
                switch sideIndex
                    case 1
                        xyz{k}=[ ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            ]';
                    case 2
                        xyz{k}=[ ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            ]';
                    case 3
                        xyz{k}=[ ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            ]';
                    case 4
                        xyz{k}=[ ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            ]';
                    case 5
                        xyz{k}=[ ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            ]';
                    case 6
                        xyz{k}=[ ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                            ]';
                    otherwise
                        error('unicycle.geometry.shearZone.exportTopSideXYZ: invalid side index.');
                end
            end
            [xp,yp,zp]=deal(xyz{:});
            
            fid=fopen(fname,'wt');
            fprintf(fid,'# export from unicycle.geometry.shearZone.exportSideXYZ\n');
            fprintf(fid,'# side index %d\n',sideIndex);
            
            for k=1:length(field)
                fprintf(fid,'> -Z%f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n', ...
                    field(k),[scale*xp(1:4,k),scale*yp(1:4,k),scale*zp(1:4,k)]');
            end
            
            fclose(fid);
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                            %
        %        export geometry to GMT .xyz format                  %
        %                                                            %
        % INPUT:                                                     %
        %                                                            % 
        % scale        - scaling factor of geometrical features      %
        % fname        - output file name                            %
        % value        - attribute to plot                           % 
        %                                                            %
        % EXAMPLE:                                                   %
        %                                                            %
        %   shz.exportXYZ(1e-3,'output/shz.xyz',value)   %
        %                                                            %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function exportXYZ(o,scale,fname,field)
            if 1>o.N
                return
            end
            
            xyz=cell(3,1);
            for k=1:3
                xyz{k}=[ ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2-o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2-o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)+o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    o.xc(:,k)-o.sv(:,k).*o.L/2+o.dv(:,k).*o.W/2+o.nv(:,k).*o.T/2, ...
                    ]';
            end
            [xp,yp,zp]=deal(xyz{:});
            
            unicycle.export.exportXYZshearZone( ...
                field,scale*xp,scale*yp,scale*zp,fname);
        end
        
    end % methods
    
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [varargout] = lvl2shz(filename)
            % LVL2FLT converts a list of levels from a file to a list
            % of shear zones
            %
            %   shz=lvl2shz(filename)
            %
            % or
            %
            %   [shz,lvl]=lvl2shz(filename)
            %
            % where shz is a list of [N,9] shear zone geometry parameters,
            % and lvl is a cell array of class 'geometry.segment'.
            [~,x1,x2,x3,len,width,thick,str,d,lo,wo,to,alphal,alphaw,alphat]=...
                textread(filename,'%u %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
                'commentstyle','shell');
            
            fm=[];
            lvl=cell(length(x1),1);
            for k=1:length(x1)
                % list of patches for current segment
                shz=unicycle.geometry.shz2shz([x1(k);x2(k);x3(k)], ...
                    len(k),width(k),thick(k),...
                    str(k),d(k),...
                    lo(k),wo(k),to(k), ...
                    alphal(k),alphaw(k),alphat(k));
                
                % level
                lvl{k}=unicycle.geometry.level(x2(k),x1(k),-x3(k),len(k),width(k),thick(k),str(k),d(k),size(fm,1),size(shz,1));
                
                % list of shear zones for all levels
                fm=[fm;shz];
            end
            
            if 1==nargout
                varargout{1}=fm;
            elseif 2==nargout
                varargout{1}=fm;
                varargout{2}=lvl;
            end
        end
    end % methods (Static)
end
