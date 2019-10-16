classdef gpsReceiver < unicycle.manifold.gps
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = gpsReceiver(network,evl,vecsize)
            % GPS is a class representing GPS data and Green's functions.
            %
            %   gps = manifold.gpsReceiver(network,evl,vecsize)
            %
            % creates a instance of GPS data connected to a simulation.
            %
            % INPUT:
            %
            % network  
            %
            % if network is a filename, for example 'sopac.dat', containing
            %
            %          # i NAME x1 x2 x3
            %            1 GPS1  0  0  0
            %            2 GPS2  1  1  0
            %
            % alternatively network may be a vector of size N x vecsize
            % containing locations for evaluating displacements
            %
            % evl      object of type ode.evolution containing a time
            %          series of fault displacement
            % vecsize  number of components of displacement vector
            %
            % SEE ALSO: unicycle
            
            import unicycle.greens.*
            
            obj.vecsize=vecsize;
            
            if isnumeric(network)
                obj.x=network;
                obj.D=size(network,1);
                obj.stationName=1:obj.D;
            else
                filename=network;
                [~,obj.stationName,x1,x2,x3]=...
                    textread(filename,'%d %s %f %f %f','commentstyle','shell');
                %obj.x=[x2,x1,0*x3];
                obj.x=[x2,x1,x3];
                obj.D=size(x2,1);
            end
            
           
            
           % source Green's functions (includes strike slip and dip slip)
            
            if isobject(evl.src)
                obj.knl.FO={'s', 'd'};
                if ~exist('./kernels/FO_s.grd','file')
                    if ~exist('./kernels/FO.mat','file')
                        [obj.FO{1,1},obj.FO{2,1}]=evl.src.displacementKernels(obj.x,vecsize);
                        for i = 1:numel(obj.FO)
                            fname=strcat('./kernels/FO_',obj.knl.FO{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],obj.FO{i},fname)
                        end
                    else
                        load('./kernels/FO.mat');
                        obj.FO=sFO;
                        clear sFO
                        for i = 1:numel(obj.FO)
                            fname=strcat('./kernels/FO_',obj.knl.FO{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],obj.FO{i},fname)
                        end
                    end
                else
                    for i = 1:numel(obj.FO)
                        fname=strcat('./kernels/FO_',obj.knl.FO{i},'.grd');
                        [~,~,obj.FO{i}]=unicycle.export.grdread(fname);
                    end
                end
            end
            
          
             % receiver Green's functions (includes strike slip and dip slip)
            
            if isobject(evl.flt)
                obj.knl.KO={'s', 'd'};
                if ~exist('./kernels/KO_s.grd','file')
                    if ~exist('./kernels/KO.mat','file')
                        [obj.KO{1,1},obj.KO{2,1}]=evl.flt.displacementKernels(obj.x,vecsize);
                        for i = 1:numel(obj.KO)
                            fname=strcat('./kernels/KO_',obj.knl.KO{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],obj.KO{i},fname)
                        end
                    else
                        load('./kernels/KO.mat');
                        obj.KO=sKO;
                        clear sKO
                        for i = 1:numel(obj.KO)
                            fname=strcat('./kernels/KO_',obj.knl.KO{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],obj.KO{i},fname)
                        end
                    end
                else
                    for i = 1:numel(obj.KO)
                        fname=strcat('./kernels/KO_',obj.knl.KO{i},'.grd');
                        [~,~,obj.KO{i}]=unicycle.export.grdread(fname);
                    end
                end
            end
                        
             % Strain Volume Green's functions (includes strike slip and dip slip)
             if isobject(evl.shz)
                 obj.knl.LO={'11','12','13','22','23','33'};
                if ~exist('./kernels/LO_11.grd','file')
                    if ~exist('./kernels/LO.mat','file')
                        [obj.LO{1,1},obj.LO{2,1},obj.LO{3,1}, ...
                            obj.LO{4,1},obj.LO{5,1},obj.LO{6,1}]=evl.shz.displacementKernels(obj.x,vecsize);
                        for i = 1:numel(obj.LO)
                            fname=strcat('./kernels/LO_',obj.knl.LO{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],obj.LO{i},fname)
                        end
                    else
                        load('./kernels/LO.mat');
                        obj.LO=sLO;
                        clear sLO
                        for i = 1:numel(obj.LO)
                            fname=strcat('./kernels/LO_',obj.knl.LO{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],obj.LO{i},fname)
                        end
                    end
                else
                    for i = 1:numel(obj.LO)
                        fname=strcat('./kernels/LO_',obj.knl.LO{i},'.grd');
                        [~,~,obj.LO{i}]=unicycle.export.grdread(fname);
                    end
                end
            end            
            
            % builds forward models of geodetic data if simulation exists
            if ~isempty(evl.y)
                obj.simulation(evl);
            end
            
        end % constructor
        
  
    end % methods
    
end
