classdef shearZoneReceiver < unicycle.geometry.shearZone
    properties
        % rheological properties
        A;Ak;n;Coh;m;etaK;Gk;etaM;Q;Vstar;P;r;d;p;
        % background confining pressure
        sigma0;
        % background temperature
        T0;
        % strain rate velocity
        epsilonPlate;
        % observation points
        observationPoints;
        % dynamic variables
        tMax;
        eMax;
        % degrees of freedom (number of parameters solved in numerical integration)
        dgf;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = shearZoneReceiver(basename,earthModel)
            % SHEARZONERECEIVER is a class representing the geometry and physical
            % properties of receiver shear zones, including position, 
            % orientation, dimension and rheological properties.
            %
            %   src = geometry.shearZoneReceiver('basename')
            %
            % or
            %
            %   src = geometry.shearZoneReceiver({'basename1','basename2'})
            %
            % where 'basename' is short for 'basename_level.shz' or
            % 'basename.shz' creates a instance of shear zones for 
            % a receiver shear zone.
            %
            % default rheological properties are:
            % 
            %   n     = 1
            %   etaM  = 1e18 Pa /s
            %   etaK  = 1e17 Pa /s
            %   Gk    = 30e3 MPa
            %   Vstar = 11e3 J/mol
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            obj=obj@unicycle.geometry.shearZone();
            
            if 0==nargin
                return
            end
            
            if ~iscell(basename)
                basename={basename};
            end
            
            obj.earthModel=earthModel;
            
            values={};
            
            fm=[];
            levels={};
            for k=1:length(basename)
                
                [wdir,base,~]=fileparts(basename{k});
                basename{k}=strcat(wdir,'/',base);
                
                fname=[basename{k} '_lvl.shz'];
                if 2==exist(fname,'file')
                    [lfm,lLevel]=obj.lvl2shz(fname);
                    if 0~=numel(lLevel)
                        for i=1:length(levels)
                            lLevel{i}.starti=lLevel{i}.starti+levels{end}.starti+levels{end}.nShearZone;
                        end
                    end
                    values=[values,length(levels)+(1:length(lLevel))];
                    fm=[fm;lfm];
                    levels=[levels;lLevel];
                else
                    fname=[basename{k} '.shz'];
                    assert(2==exist(fname,'file'),sprintf('file %s does not exist.',fname));
                    [id,x1,x2,x3,len,width,thick,str,d]=...
                        textread(fname,'%u %f %f %f %f %f %f %f %f',...
                        'commentstyle','shell');
                    obj.id=id;
                    lfm=[x1,x2,x3,len,width,thick,str,d];
                    lLevel={unicycle.geometry.level(mean(x2),mean(x1),mean(-x3),...
                        mean(len),mean(width),mean(thick),mean(str),mean(d),...
                        0,length(x1))};
                    if 0~=numel(levels)
                        lLevel{1}.starti=lLevel{1}.starti+levels{end}.starti+levels{end}.nShearZone;
                    end
                    values=[values,length(levels)+(1:length(lLevel))];
                    fm=[fm;lfm];
                    levels=[levels;lLevel];
                end
            end
            
            % map segment names to segment indices
            obj.levelHashTable=containers.Map(basename,values);
            
            % patch properties
            obj.N=size(fm,1);
            obj.x=[fm(:,[2,1]),-fm(:,3)];
            obj.L=fm(:,4);
            obj.W=fm(:,5);
            obj.T=fm(:,6);
            obj.strike=fm(:,7);
            obj.dip=fm(:,8);
            
            % default viscoelastic properties
            obj.n=obj.L*0+1;
            obj.m=obj.L*0+1;
            obj.etaM=obj.L*0+1e18/(365*24*60*60)/1e6;
            obj.etaK=obj.L*0+1e17/(365*24*60*60)/1e6;
            obj.Gk=obj.L*0+30e3;
            
            % plate rate
            obj.epsilonPlate=obj.L*0+1e-15;

            % segments
            obj.levels=levels;
            
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
        
        function set.observationPoints(obj,points)
            % function SET.OBSERVATIONPOINTS adds information to the
            % observation points based on receiver elements.
            %
            % SEE ALSO: unicycle.
            
            if ~iscell(points)
                points={points};
            end
            
            for k=1:length(indepoints)
                points{k}.xc=obj.xc(points{k}.index);
            end
            
            obj.observationPoints=points;
        end
        
    end % methods
    
end
