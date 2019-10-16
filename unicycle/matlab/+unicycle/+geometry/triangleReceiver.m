classdef triangleReceiver < unicycle.geometry.triangle
    properties
        % friction properties
        Vo;a;b;l;mu0;
        % shear wave speed
        Vs;
        % confining pressure
        sigma;
        % plate velocity
        Vpl;
        % rake of plate rate
        Vrake;
        % is rake constrained?
        isRakeConstraint;
        % fixed rake positions (index)
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
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = triangleReceiver(varargin)
            % TRIANGLERECEIVER is a class representing the geometry and 
            % physicalproperties of receiver fault patches, including 
            % position, orientation, dimension and friction properties.
            %
            %   src = geometry.triangleReceiver('basename')
            %
            % or
            %
            %   src = geometry.triangleReceiver({'basename1','basename2'})
            %
            % where 'basename' is short for 'basename_seg.flt' or
            % 'basename_patch.flt' creates a instance of fault patches for 
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
            % SEE ALSO: unicycle, unicycle.geometry.receiver
            
            import unicycle.geometry.triangle;
            
            if isempty(varargin)
                return
            else
                basename=varargin{1};
                obj.earthModel=varargin{2};
            end
            
            if ~iscell(basename)
                basename={basename};
            end
            
            obj.x=[];
            obj.vertices=[];
            obj.rake=[];
            for k=1:length(basename)
                fname=[basename{k} '.ned'];
                assert(2==exist(fname,'file'),sprintf('error: can''t find %s',fname));
                [~,x1,x2,x3]=...
                    textread(fname,'%u %f %f %f','commentstyle','shell');
                obj.x=[obj.x;[x2,x1,-x3]];
                
                assert(0>=max(obj.x(:,3)),'error: all vertices should have positive depth.');
                
                fname=[basename{k} '.tri'];
                if 2==exist(fname,'file')
                    [Vpl,tri,rake]=obj.loadTri(fname);
                    obj.vertices=[obj.vertices;tri];
                    obj.rake=[obj.rake;rake(:)];
                    obj.Vpl=[obj.Vpl;Vpl(:)];
                else
                    fprintf('missing %s file: building Delaunay mesh.\n',fname)
                    tri=delaunayn([x2,x1,-x3],{'QJ'});
                    obj.vertices=[obj.vertices;tri];
                    obj.rake=[obj.rake;zeros(size(tri,1),1)];
                end
            end
            
            % patch properties
            obj.N=size(obj.vertices,1);
            obj.id=1:obj.N;
            obj.slip=zeros(obj.N,1);
            
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
            
            % no +-90º rake constraint by default
            obj.isRakeConstraint=false;
            % no patches with fixed rake by default
            obj.FixedRakePosition=[];
            
            % default friction properties (velocity strengthening)
            obj.a=zeros(obj.N,1)+1e-2;
            obj.b=obj.a-4e-3;
            obj.Vo=zeros(obj.N,1)+1e-1;
            obj.l=zeros(obj.N,1)+1e-3;
            obj.sigma=zeros(obj.N,1)+1e2;
            obj.mu0=zeros(obj.N,1)+0.6;
            obj.Vs=obj.l*0+3e3*3.1536e7; % (m/yr)

            % plate rate
            obj.Vpl=obj.l*0;
            obj.Vrake=obj.l*0;

            % unit vectors
            [obj.nv,obj.sv,obj.dv,obj.area]=triangle.computeUnitVectors(obj.x,obj.vertices);

        end % constructor
        
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
        
        function set.observationPoints(obj,points)
            % function SET.OBSERVATIONPOINTS adds information to the
            % observation points based on receiver elements.
            %
            % SEE ALSO: unicycle.
            
            if ~iscell(points)
                points={points};
            end
            
            for k=1:length(points)
                points{k}.xc=obj.xc(points{k}.index,:);
            end
            
            obj.observationPoints=points;
        end
        
        function set.eventCatalogue(obj,eCat)
            % function SET.OBSERVATIONPOINTS adds information to the
            % observation points based on receiver elements.
            %
            % SEE ALSO: unicycle.
            
            if ~iscell(eCat)
                eCat={eCat};
            end
            
            obj.eventCatalogue=eCat;
        end
        
    end % methods
    
end
