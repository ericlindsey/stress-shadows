classdef receiver < unicycle.geometry.patch
    properties
        % friction properties
        Vo;a;b;l;
        % shear wave speed
        Vs;
        % initial confining pressure
        sigma;
        % initial shear stress
        tau;
        % friction coefficient at reference velocity
        mu0;
        % plate velocity
        Vpl;
        % rake of plate rate
        Vrake;
        % is rake constrained?
        isRakeConstraint;
        % fixed rake positions
        FixedRakePosition;
        % pinned patch positions (index)
        pinnedPosition;
        % observation points
        observationPoints;
        % event catalogue
        eventCatalogue;
        % dynamic variables
        tMax;
        vMax;
        % degrees of freedom (number of parameters solved in numerical integration)
        dgf;
        %pinned patches modified by Eric - July 25, 2017
        isfree;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = receiver(filename,earthModel)
            % RECEIVER is a class representing the geometry and physical
            % properties of receiver fault patches, including position, 
            % orientation, dimension and friction properties.
            %
            %   src = geometry.receiver('filename')
            %
            % or
            %
            %   src = geometry.receiver({'filename1','filename2'})
            %
            % where 'filename' is short for 'filename.seg' or
            % 'filename.flt' creates a instance of fault patches for 
            % a receiver fault.
            %
            % default friction properties are velocity strengthening:
            % 
            %   V0    = 10 m/yr
            %   a     = 1e-2
            %   b     = 6e-3
            %   l     = 1e-3 m
            %   sigma = 100 MPa.
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            if ~iscell(filename)
                filename={filename};
            end
            
            obj.earthModel=earthModel;
            
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
                        values={values,length(seg)+(1:length(lseg))};
                        fm=[fm;lfm];
                        seg=[seg;lseg];
                        
                        % map segment names to segment indices
                        obj.segment=containers.Map(filename{k},values);
                    case '.flt'
                        [lfm]=obj.loadFlt(filename{k});
                        
                        fm=[fm;lfm];
                        
                    otherwise
                        error(['unicycle:geometry:patch:wrong file type ' filename{k}]);
                end
            end
            
            
            % patch properties     
           
            switch size(fm,2)
                
                case 9 %if there is Vpl in receiver
                    obj.N=size(fm,1);
                    obj.slip=zeros(obj.N,1);
                    obj.x=[fm(:,[3,2]),-fm(:,4)];
                    obj.L=fm(:,5);
                    obj.W=fm(:,6);
                    obj.strike=fm(:,7);
                    obj.dip=fm(:,8);
                    obj.Vrake=fm(:,9);
                    obj.Vpl=fm(:,1);
                    
                case 8
                    obj.N=size(fm,1);
                    obj.slip=zeros(obj.N,1);
                    obj.x=[fm(:,[2,1]),-fm(:,3)];
                    obj.L=fm(:,4);
                    obj.W=fm(:,5);
                    obj.strike=fm(:,6);
                    obj.dip=fm(:,7);
                    obj.Vrake=fm(:,8);
                    obj.Vpl=obj.L*0;
                    
                otherwise
                    error('unicycle:geometry:patch:receiver:coding error');
            end
            
            
            % no +-90º rake constraint by default
            obj.isRakeConstraint=false;
            % no patches with fixed rake by default
            obj.FixedRakePosition=[];
            
            % default friction properties (velocity strengthening)
            obj.a=obj.L*0+1e-2;
            obj.b=obj.a-4e-3;
            obj.Vo=obj.L*0+1e-1;
            obj.l=obj.L*0+1e-3;
            obj.sigma=obj.L*0+1e2;
            obj.mu0=obj.L*0+0.6;
            obj.tau=obj.mu0.*obj.sigma;
            obj.Vs=obj.L*0+3e3*3.1536e7; % (m/yr)
            
            
            obj.Vrake=obj.L*0;

            % segments
            obj.segments=seg;
            
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
        
        function [tri]=toTriangleReceiver(obj,varargin)
            % TOTRIANGLERECEIVER provides a triangle class representation
            % of the fault plane with two triangles per rectangle patch.
            %
            %   [triangleReceiver] = obj.toTriangleReceiver()
            %
            % OUTPUT:
            %   triangleReceiver - an object of class triangleReceiver.
            
            import unicycle.geometry.triangle;
            
            if isempty(varargin)
                tri=unicycle.geometry.triangleReceiver();
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
            tri.Vrake=[obj.Vrake;obj.Vrake];
            
            % friction properties
            tri.Vo=[obj.Vo;obj.Vo];
            tri.a=[obj.a;obj.a];
            tri.b=[obj.b;obj.b];
            tri.l=[obj.l;obj.l];
            tri.mu0=[obj.mu0;obj.mu0];
            tri.Vs=[obj.Vs;obj.Vs];
            tri.Vpl=[obj.Vpl;obj.Vpl];
            
            % confining pressure
            tri.sigma=[obj.sigma;obj.sigma];
            % is rake constrained?
            tri.isRakeConstraint=[obj.isRakeConstraint;obj.isRakeConstraint];
            % fixed rake positions
            tri.FixedRakePosition=[obj.FixedRakePosition;obj.FixedRakePosition];
            
        end
        
        function set.observationPoints(obj,points)
            % function SET.OBSERVATIONPOINTS adds information to the
            % observation points based on receiver elements.
            %
            % SEE ALSO: unicycle.
            
            if ~iscell(points)
                points={points};
            end
            
            for k=1:length(points)
                points{k}.xc=obj.xc(points{k}.index);
            end
            
            obj.observationPoints=points;
        end

        %% % % % % % % % % % % % % % % % % % % % %
        %                                        %
        %  c o n s t r a i n t s   m a t r i x   %
        %                                        %
        %                f o r                   %
        %                                        %
        %   s t a t i c   i n v e r s i o n s    %
        %                                        %
        % % % % % % % % % % % % % % % % % % % % %%
        function [A,B] = rakeConstraintMatrix(obj)
            % function RAKECONSTRAINTMATRIX produces regularization
            % and constraint matrices for non-negative static inversion
            % of the kind 
            %
            %   mt = lsqlin(H,h,A,a) solving the least-squares problem
            %
            %      min  (NORM(H*m-h)).^2   subject to  A*x <= a
            %       m
            %
            % where H=[G;B] and h=[d;b].
            %
            
            % matrix for non-negative constraints
            A=[];
            
            % matrix for penalization of rake-perpendicular slip
            B=zeros(length(obj.FixedRakePosition),2*obj.N);
            for k=1:length(obj.FixedRakePosition)
                B(k,      obj.FixedRakePosition(k))=cosd(obj.rake(obj.FixedRakePosition(k))+90);
                B(k,obj.N+obj.FixedRakePosition(k))=sind(obj.rake(obj.FixedRakePosition(k))+90);
            end
        end
        
    end % methods
    
end
