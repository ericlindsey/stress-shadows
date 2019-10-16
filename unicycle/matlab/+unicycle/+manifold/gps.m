classdef gps < handle
    properties
        % number of GPS stations
        D;
        % station name
        stationName;
        % position
        x;
        % number of GPS components
        vecsize;
        % fault receiver displacement Green's functions
        KO=cell(2,1);
        % shear zone displacement Green's functions
        LO=cell(6,1);
        % source Green's functions
        FO=cell(1,1);
        % time series of displacement due to source
        us;
        % time series of displacement due to fault slip
        ua;
        % time series of displacement due to shear zones
        uv;
        % data vector
        d;
        % list of epochs
        epochs;
        % time coordinates
        t;
        %Kernel labels
        knl;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = gps(network,evl,vecsize)
            % GPS is a class representing GPS data.
            %
            %   gps = data.gps(network,evl,vecsize)
            %
            % creates a instance of GPS data.
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
            % alternatively network can be a vector of size N x vecsize
            % containing locations for evaluating displacements
            %
            % evl      object of type ode.evolution containing a time
            %          series of fault displacement
            % vecsize  number of components of displacement vector
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            if (0==nargin)
                return
            end
            
            obj.vecsize=vecsize;
                      
            if isnumeric(network)
                obj.x=network;
                obj.D=size(network,1);
                obj.stationName=1:obj.D;
            else
                filename=network;
                [~,obj.stationName,x1,x2,x3]=...
                    textread(filename,'%d %s %f %f %f','commentstyle','shell');
                obj.x=[x2,x1,0*x3];
                obj.D=size(x2,1);
            end
            
            
        end % constructor
        
        function [] = simulation(obj,evl)
            % method SIMULATION builds forward models of geodetic data if
            % simulation exists
            
            % save the time coordinates
            %obj.t=evl.t;
            
            % surface displacements due to steadily moving source faults
            if isobject(evl.src)
                if 0<evl.src.N
                    obj.us=zeros(obj.vecsize*obj.D,numel(evl.t));
                    for k=1:length(evl.t)
                        obj.us(:,k)= ...
                             obj.FO{1,1}*(evl.src.slip.*cosd(evl.src.rake)*(evl.t(k)-evl.t(1))) ...
                            +obj.FO{2,1}*(evl.src.slip.*sind(evl.src.rake)*(evl.t(k)-evl.t(1)));
                    end
                end
            end
            
            % surface displacements due to fault slip
            if isobject(evl.flt)
                obj.ua=zeros(obj.vecsize*obj.D,numel(evl.t));
                for k=1:length(evl.t)
                    % fault slip
                    obj.ua(:,k)=obj.KO{1,1}*evl.y(1:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),k) ...
                               +obj.KO{2,1}*evl.y(2:evl.flt.dgf:(evl.flt.N*evl.flt.dgf),k);
                end
            end
            
            % surface displacements due to strain on shear zones
            if isobject(evl.shz)
                obj.uv=zeros(obj.vecsize*obj.D,numel(evl.t));
                for k=1:length(evl.t)
                    obj.uv(:,k)=obj.LO{1,1}*evl.y(evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end,k) ...
                               +obj.LO{2,1}*evl.y(evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end,k) ...
                               +obj.LO{3,1}*evl.y(evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end,k) ...
                               +obj.LO{4,1}*evl.y(evl.flt.N*evl.flt.dgf+4:evl.shz.dgf:end,k) ...
                               +obj.LO{5,1}*evl.y(evl.flt.N*evl.flt.dgf+5:evl.shz.dgf:end,k) ...
                               +obj.LO{6,1}*evl.y(evl.flt.N*evl.flt.dgf+6:evl.shz.dgf:end,k);
                end
            end
        
            
        end % simulation
        
        function [] = exportTimeSeriesToNED(obj,t,wdir)
            
            for stationId=1:obj.D
                u1=0;
                u2=0;
                u3=0;
                if ~isempty(obj.ua)
                    u1=u1+obj.ua((stationId-1)*obj.vecsize+2,:)';
                    u2=u2+obj.ua((stationId-1)*obj.vecsize+1,:)';
                    u3=u3-obj.ua((stationId-1)*obj.vecsize+3,:)';
                end
                if ~isempty(obj.uv)
                    u1=u1+obj.uv((stationId-1)*obj.vecsize+2,:)';
                    u2=u2+obj.uv((stationId-1)*obj.vecsize+1,:)';
                    u3=u3-obj.uv((stationId-1)*obj.vecsize+3,:)';
                end
                
                fname=[wdir '/' obj.stationName{stationId} '.ned'];
                fid=fopen(fname,'wt');
                fprintf(fid,'# time    u1    u2    u3\n');
                fprintf(fid,'%f %f %f %f\n',[t,u1,u2,u3]');
                fclose(fid);
            end
        end
        function [dout] = exportTimeseriesArray(obj)
            dout = [];
            for stationId=1:obj.D
                u1=0;
                u2=0;
                u3=0;
                if ~isempty(obj.ua)
                    u1=u1+obj.ua((stationId-1)*obj.vecsize+2,:)';
                    u2=u2+obj.ua((stationId-1)*obj.vecsize+1,:)';
                    u3=u3-obj.ua((stationId-1)*obj.vecsize+3,:)';
                end
                if ~isempty(obj.uv)
                    u1=u1+obj.uv((stationId-1)*obj.vecsize+2,:)';
                    u2=u2+obj.uv((stationId-1)*obj.vecsize+1,:)';
                    u3=u3-obj.uv((stationId-1)*obj.vecsize+3,:)';
                end
                
                dout=[dout;u1];
                dout=[dout;u2];
                dout=[dout;u3];
            end
        end
        
        function [] = plotTimeSeries(obj,t,stationId,displacementComponent,varargin)
            
            if 4<nargin
                style=varargin{1};
            else
                style='b-';
            end

            disp=0;
            if ~isempty(obj.ua)
                disp=disp+obj.ua((stationId-1)*obj.vecsize+displacementComponent,:)';
            end
            if ~isempty(obj.uv)
                disp=disp+obj.uv((stationId-1)*obj.vecsize+displacementComponent,:)';
            end
            if ~isempty(obj.us)
                disp=disp+obj.us((stationId-1)*obj.vecsize+displacementComponent,:)';
            end
            
            %tG=[ones(size(t,1),1), t];
            %m=(tG'*tG)\(tG'*disp);
            %disp=disp-tG*m;
            plot(t,disp,style,'linewidth',5);
            
        end
        
        function d = dataVectorFromTimeSeries(obj,reduction,tlim)
            % method DATAVECTORFROMTIMESERIES constructs a data vector from
            % GPS time series
            %
            %    d = gps.dataVectorFromTimeSeries(reduction,tlim)
            %
            % and filters data out of the range [tlim(1) tlim(end)];
            
            assert(obj.D==length(obj.stationName),...
                'unicycle:manifold:gps:incorrect number of GPS stations');
            
            assert(3>=obj.vecsize,...
                'unicycle:manifold:gps:no implementation for vecsize > 3');
            
            obj.epochs=cell(obj.D,1);
            N=0;
            for k=1:obj.D
                try
                    filename=[reduction '/' obj.stationName{k}];
                    [gt]=textread(filename,'%f %*[^\n]','commentstyle','shell');
                catch ME
                    filename=[reduction '/' obj.stationName{k} '.ned'];
                    [gt]=textread(filename,'%f %*[^\n]','commentstyle','shell');
                end
                pos=find(gt>=tlim(1) & gt<= tlim(end));
                obj.epochs{k}=gt(pos);
                N=N+length(pos);
            end
            
            d=zeros(N*obj.vecsize,1);
            N=0;
            for k=1:obj.D
                if 1==obj.vecsize
                    try
                        filename=[reduction '/' obj.stationName{k}];
                        [gt,gx]=...
                            textread(filename,'%f %f %*[^\n]','commentstyle','shell');
                    catch ME
                        filename=[reduction '/' obj.stationName{k} '.ned'];
                        [gt,gx]=...
                            textread(filename,'%f %f %*[^\n]','commentstyle','shell');
                    end
                    pos=find(gt>=tlim(1) & gt<= tlim(end));
                    n=length(pos);
                    u=gx(pos);
                end
                if 2==obj.vecsize
                    try
                        filename=[reduction '/' obj.stationName{k}];
                        [gt,gy,gx]=...
                            textread(filename,'%f %f %f %*[^\n]','commentstyle','shell');
                    catch ME
                        filename=[reduction '/' obj.stationName{k} '.ned'];
                        [gt,gy,gx]=...
                            textread(filename,'%f %f %f %*[^\n]','commentstyle','shell');
                    end
                    pos=find(gt>=tlim(1) & gt<= tlim(end));
                    n=length(pos);
                    u=[gx(pos),gy(pos)];
                end
                if 3==obj.vecsize
                    try
                        filename=[reduction '/' obj.stationName{k}];
                        [gt,gy,gx,gd]=...
                            textread(filename,'%f %f %f %f %*[^\n]','commentstyle','shell');
                    catch ME
                        filename=[reduction '/' obj.stationName{k} '.ned'];
                        [gt,gy,gx,gd]=...
                            textread(filename,'%f %f %f %f %*[^\n]','commentstyle','shell');
                    end
                    pos=find(gt>=tlim(1) & gt<= tlim(end));
                    n=length(pos);
                    u=[gx(pos),gy(pos),-gd(pos)];
                end
                for j=1:obj.vecsize
                    d(N+1:N+n)=u(:,j);
                    N=N+n;
                end
            end
        end
        
        function d = dataVectorFromField(obj)
            d=[];
        end
        
        function g = dataVectorFromModel(obj,evl)
            % method DATAVECTORFROMMODEL creates a simulation data vector
            %
            %    g = gps.dataVectorFromModel(evl)
            %
            % for direct comparison with the data vector. Time series of
            % simulated data is interpolated at the time of GPS epochs.
            
            % builds forward models of geodetic data if simulation exists
            
            % this build obj.ur and obj.us
            obj.simulation(evl);
            
            g=zeros(size(obj.d));
            N=0;
            for k=1:obj.D
                for j=1:obj.vecsize
                    n=length(obj.epochs{k});
                    % interpolate unicycle simulation at time of GPS epochs
                    disp=0;
                    if ~isempty(obj.ua)
                        disp=disp+obj.ua((k-1)*obj.vecsize+j,:);
                    end
                    if ~isempty(obj.uv)
                        disp=disp+obj.uv((k-1)*obj.vecsize+j,:);
                    end
                    if ~isempty(obj.us)
                        disp=disp+obj.us((k-1)*obj.vecsize+j,:);
                    end
                    g(N+1:N+n)=interp1(evl.t,disp,obj.epochs{k});
                    N=N+n;
                end
            end
        end
        
        function d = dataVectorAtTime(obj,t,varargin)
            % method DATAVECTORATTIME constructs a data vector of static
            % displacement at a given time.
            %
            %    d = gps.dataVectorAtTime(t)
            %
            
            if 2<nargin
                interp=true;
            else
                interp=false;
            end
            
            d=zeros(obj.D*obj.vecsize,1);
            N=0;
            for k=1:obj.D
                dt=obj.epochs{k};
                n=length(dt);
                
                for j=1:obj.vecsize
                    if min(dt)>t
                        data=0;
                    else
                        if max(dt)<t
                            data=obj.d(N+n);
                        else
                            if interp
                                Gparabolic=[ones(length(dt),1) dt dt.^2 dt.^3];
                                Gp=[1 t t.^2 t.^3];
                                data=Gp*(pinv(Gparabolic'*Gparabolic)*(Gparabolic'*obj.d(N+1:N+n)));
                            else
                                data=interp1(obj.epochs{k},obj.d(N+1:N+n),t);
                            end
                        end 
                    end
                    % interpolate unicycle simulation at time of GPS epochs
                    d(1+(j-1)+(k-1)*obj.vecsize)=data;
                    N=N+n;
                end
            end
        end
        
        function g = dataVectorFromModelAtTime(obj,evl,t)
            % method DATAVECTORFROMMODELATTIME creates a simulation data vector
            %
            %    g = gps.dataVectorFromModelAtTime(evl,t)
            %
            % for direct comparison with the data vector. Time series of
            % simulated data is interpolated at the time of GPS epochs.
            
            % builds forward models of geodetic data if simulation exists
            
            % assumes obj.ur and obj.us created by obj.dataVectorFromModel
            
            if nargin<3
                error('unicycle::manifold:gps::NotEnoughArguments');
            end
            
            g=zeros(obj.D*obj.vecsize,1);
            for k=1:obj.D
                for j=1:obj.vecsize
                    % interpolate unicycle simulation at time specified
                    disp=0;
                    if ~isempty(obj.ua)
                        disp=disp+obj.ua((k-1)*obj.vecsize+j,:);
                    end
                    if ~isempty(obj.uv)
                        disp=disp+obj.uv((k-1)*obj.vecsize+j,:);
                    end
                    g(1+(j-1)+(k-1)*obj.vecsize)=interp1(evl.t,disp,t);
                    if ~isempty(obj.us)
                        g(1+(j-1)+(k-1)*obj.vecsize)=g(1+(j-1)+(k-1)*obj.vecsize)+interp1(evl.t,obj.us((k-1)*obj.vecsize+j,:),t);
                    end
                end
            end
        end
        
    end % methods
    
end
