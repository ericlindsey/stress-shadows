classdef evolution < handle
    properties
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                   %
        %        S T R E S S   K E R N E L   F R O M        %
        %      R E C E I V E R   T O   R E C E I V E R      %
        %                                                   %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        % fault self stress
        % KK{i,j} is traction s{j} due to fault slip d{i}, where
        % s={strike shear, dip shear, normal stress} and d={strike slip, dip slip}.
        KK=cell(2,3); % cell(causes, receivers)
        
        % stress on shear zones due to fault slip
        % KL{i,j} is stress s{j} due to slip d{i}, where d{i}={strike slip, dip slip}
        % and s=[s11,s12,s13,s22,s23,s33]
        KL=cell(2,6); % cell(causes, receivers)
        
        % traction on fault due to shear zones
        % LK{i,j} is traction s{j} due to strain e{i}, where 
        % e=[e11,e12,e13,e22,e23,e33] and s={strike shear, dip shear, normal stress}
        LK=cell(6,3); % cell(causes, receivers)

        % shear zones self stress
        % LL{i,j} is stress s{j} due to e{i} where s=[s11,s12,s13,s22,s23,s33]
        % and e=[e11,e12,e13,e22,e23,e33]
        LL=cell(6,6); % cell(causes, receivers)
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                   %
        %        S T R E S S   K E R N E L   F R O M        %
        %        S O U R C E   T O   R E C E I V E R        %
        %                                                   %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % fault stress due to continuous source
        % FK{i,j} is traction s{j} due to fault slip d{i}, where
        % s={strike shear, dip shear, normal stress} and d={strike slip, dip slip}.
        FK=cell(2,3); % cell(causes, receivers)
        
        % stress on shear zones due to continuous fault slip
        % FL{i,j} is stress s{j} due to slip d{i}, where d{i}={strike slip, dip slip}
        % and s=[s11,s12,s13,s22,s23,s33]
        FL=cell(2,6); % cell(causes, receivers)
        
        % coseismic events with slip distribution and stress kernels
        evt;
        % source geometry and slip distribution
        src;
        % fault receiver geometry
        flt;
        % fault receiver geometry
        shz;
        % simulation result (time and model)
        t;y;
        
        %kernel labels
        knl;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=evolution(src,flt,shz,evts)
            % EVOLUTION is a class encapsulating the geometry and physical
            % properties of fault patches for sources and receivers and
            % methods for solving integral equations for fault slip.
            %
            %   evl = ode.evolution(src,rcv);
            %
            % creates a instance from source (src), receiver (rcv).
            %
            % SEE ALSO: unicycle
            
            if (0==nargin)
                o.src=struct('N',0,'dgf',0);
                o.shz=struct('N',0,'dgf',0);
                o.flt=struct('N',0,'dgf',0);
                return
            end
            
            if isempty(src)
                o.src=struct('N',0,'dgf',0);
            end
            
            if isempty(shz)
                o.shz=struct('N',0,'dgf',0);
            end
            
            if isempty(flt)
                o.flt=struct('N',0,'dgf',0);
            end
            
            % receiver faults
            o.flt=flt;
            % receiver shear zones
            o.shz=shz;
            % source faults
            o.src=src;
            if ~exist('./kernels/','dir')
                mkdir('kernels');
            end
            
            o.knl.KK={'ss', 'sd', 'sn'; 'ds', 'dd', 'dn'};
            % stress interactions
            if isobject(flt)
                if ~exist('./kernels/KK_ss.grd','file')
                    if ~exist('./kernels/KK.mat','file')
                        [o.KK{:}]=flt.tractionKernels(flt);
                        for i = 1:numel(o.KK)
                            fname=strcat('./kernels/KK_',o.knl.KK{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],o.KK{i},fname)
                        end
                    else
                        load('./kernels/KK.mat');
                        o.KK=sKK;
                        clear sKK
                        for i = 1:numel(o.KK)
                            fname=strcat('./kernels/KK_',o.knl.KK{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],o.KK{i},fname)
                        end                    
                    end
                else
                    for i = 1:numel(o.KK)
                        fname=strcat('./kernels/KK_',o.knl.KK{i},'.grd');
                        [~,~,o.KK{i}]=unicycle.export.grdread(fname);
                    end
                end
            end
            if isobject(shz)
                if isobject(flt)
                    o.knl.KL={'s11', 's12', 's13','s22', 's23', 's33';...
                        'd11', 'd12', 'd13','d22', 'd23', 'd33'};
                    if ~exist('./kernels/KL_s11.grd','file')
                        if ~exist('./kernels/KL.mat','file')
                            [o.KL{:}]=flt.stressKernels(shz);
                            for i = 1:numel(o.KL)
                                fname=strcat('./kernels/KL_',o.knl.KL{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],o.KL{i},fname)
                            end
                        else
                            load('./kernels/KL.mat');
                            o.KL=sKL;
                            clear sKL
                            for i = 1:numel(o.KL)
                                fname=strcat('./kernels/KL_',o.knl.KL{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],o.KL{i},fname)
                            end
                        end
                    else
                        for i = 1:numel(o.KL)
                            fname=strcat('./kernels/KL_',o.knl.KL{i},'.grd');
                            [~,~,o.KL{i}]=unicycle.export.grdread(fname);
                        end
                    end
                    o.knl.LK={'s11', 'd11', 'n11'; 's12', 'd12', 'n12'; ...
                        's13', 'd13', 'n13'; 's22', 'd22', 'n22'; ...
                        's23', 'd23', 'n23'; 's33', 'd33', 'n33'};
                    if ~exist('./kernels/LK_s11.grd','file')
                        if ~exist('./kernels/LK.mat','file')
                            [o.LK{:}]=shz.tractionKernels(flt);
                            for i = 1:numel(o.LK)
                                fname=strcat('./kernels/LK_',o.knl.LK{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],o.LK{i},fname)
                            end
                        else
                            load('./kernels/LK.mat');
                            o.LK=sLK;
                            clear sLK
                            for i = 1:numel(o.LK)
                                fname=strcat('./kernels/LK_',o.knl.LK{i},'.grd');
                                unicycle.export.grdwrite([0 1], [0 1],o.LK{i},fname)
                            end
                        end
                    else
                        for i = 1:numel(o.LK)
                            fname=strcat('./kernels/LK_',o.knl.LK{i},'.grd');
                            [~,~,o.LK{i}]=unicycle.export.grdread(fname);
                        end
                    end
                end
                o.knl.LL={'1111','1211','1311','2211','2311','3311';...
                    '1112','1212','1312','2212','2312','3312';...
                    '1113','1213','1313','2213','2313','3313';...
                    '1122','1222','1322','2222','2322','3322';...
                    '1123','1223','1323','2223','2323','3323';...
                    '1133','1233','1333','2233','2333','3333'};
                if ~exist('./kernels/LL_1111.grd','file')
                    if ~exist('./kernels/LL.mat','file')
                        [o.LL{:}]=shz.stressKernels(shz);
                        for i = 1:numel(o.LL)
                            fname=strcat('./kernels/LL_',o.knl.LL{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],o.LL{i},fname)
                        end
                    else
                        load('./kernels/LL.mat');
                        o.LL=sLL;
                        clear sLL
                        for i = 1:numel(o.LL)
                            fname=strcat('./kernels/LL_',o.knl.LL{i},'.grd');
                            unicycle.export.grdwrite([0 1], [0 1],o.LL{i},fname)
                        end
                    end
                else
                    for i = 1:numel(o.LL)
                        fname=strcat('./kernels/LL_',o.knl.LL{i},'.grd');
                        [~,~,o.LL{i}]=unicycle.export.grdread(fname);
                    end
                end
            end
            
            
            % source/receiver stress interactions
            if isobject(src)
                if 1>o.src.N
                    o.FK=cell(2,3);
                    o.knl.FK={'ss', 'sd', 'sn'; 'ds', 'dd', 'dn'};
                    o.FL=cell(2,6);
                    o.knl.FL={'s11', 's12', 's13','s22', 's23', 's33';...
                        'd11', 'd12', 'd13','d22', 'd23', 'd33'};
                else
                    if isobject(flt)
                        if ~exist('./kernels/FK_ss.grd','file')
                            if ~exist('./kernels/FK.mat','file')
                                [o.FK{:}]=src.tractionKernels(flt);
                                for i = 1:numel(o.FK)
                                    fname=strcat('./kernels/FK_',o.knl.FK{i},'.grd');
                                    unicycle.export.grdwrite([0 1], [0 1],o.FK{i},fname)
                                end
                            else
                                load('./kernels/FK.mat');
                                o.FK=sFK;
                                clear sFK
                                for i = 1:numel(o.FK)
                                    fname=strcat('./kernels/FK_',o.knl.FK{i},'.grd');
                                    unicycle.export.grdwrite([0 1], [0 1],o.FK{i},fname)
                                end
                            end
                        else
                            for i = 1:numel(o.FK)
                                fname=strcat('./kernels/FK_',o.knl.FK{i},'.grd');
                                [~,~,o.FK{i}]=unicycle.export.grdread(fname);
                            end
                        end
                    end
                    if isobject(shz)
                        if ~exist('./kernels/FL_s11.grd','file')
                            if ~exist('./kernels/FL.mat','file')
                                [o.FL{:}]=src.stressKernels(shz);
                                for i = 1:numel(o.FL)
                                    fname=strcat('./kernels/FL_',o.knl.FL{i},'.grd');
                                    unicycle.export.grdwrite([0 1], [0 1],o.FL{i},fname)
                                end
                            else
                                load('./kernels/FL.mat');
                                o.FL=sFL;
                                clear sFL
                                for i = 1:numel(o.FL)
                                    fname=strcat('./kernels/FL_',o.knl.FL{i},'.grd');
                                    unicycle.export.grdwrite([0 1], [0 1],o.FL{i},fname)
                                end

                            end
                        else
                            for i = 1:numel(o.FL)
                                fname=strcat('./kernels/FL_',o.knl.FL{i},'.grd');
                                [~,~,o.FL{i}]=unicycle.export.grdread(fname);
                            end
                        end
                    end
                end
            end
            
            % coseismic events with prescribed slip distribution and timing
            o.evt=unicycle.ode.evolution.eventKernels(evts,flt,shz);
        end
        
        function []=ode(~,~,~)
            % ODE is a method that returns a column vector corresponding
            % to f(t,y) in order to solve the system of differential
            % equations y' = f(t,y).
            %
            % This is a meta method meant to be overloaded by a child
            % class.
        end
        
        function []=exportvtp(obj,t,y,scale,jump,wdir)
            % EXPORTVTK write a series of .vtp files to disk to visualize
            % the time series of fault slip in Paraview.
            %
            %   evl.exportvtp(evl.t,evl.y,1e-3,'xyz')
            %
            % input:
            % t         - time
            % y         - solution vector
            % scale     - spatial scale factor for fault position and dimension
            % jump      - increment between time steps
            % wdir      - export directory
            %
            % SEE ALSO: unicycle
            
            assert(length(t)==size(y,2));
            
            s2y=60*60*24*365;
            y2s=1./s2y;
            
            [xp,yp,zp,dim]=obj.rcv.computeVertexPosition();
            
            for k=1:jump:length(t)
                ss=y(1:obj.dgf:end,k);
                ds=y(2:obj.dgf:end,k);
                s=sqrt(ss.^2+ds.^2);
                if 1==k
                    if 2<=length(t)
                        vs=(y(1:obj.dgf:end,k+1)-y(1:obj.dgf:end,k+1))/(t(k+1)-t(k))*y2s;
                        vd=(y(2:obj.dgf:end,k+1)-y(2:obj.dgf:end,k+1))/(t(k+1)-t(k))*y2s;
                    else
                        vs=ss*0;
                        vd=ss*0;
                    end
                else
                    vs=(ss-y(1:obj.dgf:end,k-1))/(t(k)-t(k-1))*y2s;
                    vd=(ds-y(2:obj.dgf:end,k-1))/(t(k)-t(k-1))*y2s;
                end
                v=sqrt(vs.^2+vd.^2);

                fname=sprintf('%s/receiver-%05d.vtp',wdir,1+(k-1)/jump);
                unicycle.export.exportvtk_rfaults( ...
                    scale*xp, ...
                    scale*yp, ...
                    scale*zp, ...
                    dim, ...
                    fname,...
                    'strike slip',ss, ...
                    'dip slip',ds, ...
                    'slip',s,...
                    'strike velocity',vs, ...
                    'dip velocity',vd, ...
                    'velocity',v, ...
                    'log10 velocity',log10(max(v,1e-15)),...
                    'strike shear',y(3:obj.dgf:end,k), ...
                    'dip shear',y(4:obj.dgf:end,k))
            end
            
        end % end exportvtk
        
        function []=exportflt(obj,t,y,scale,wdir,varargin)
            % EXPORTFLT write a series of .flt files, one for each time step.
            %
            %   evl.exportflt(evl.t,evl.y,1e-3,'xyz')
            %
            % or
            %
            %   evl.exportflt(evl.t,evl.y,1e-3,'xyz',k)
            %
            % exports time step k.
            %
            % input:
            % t         - time
            % y         - solution vector
            % scale     - spatial scale factor for fault position and dimension
            % wdir      - export directory
            %
            % SEE ALSO: unicycle
            
            assert(length(t)==size(y,2));
            
            if nargin > 1
                k=varargin{1};
                ss=y(1:obj.dgf:end,k);
                ds=y(2:obj.dgf:end,k);
                s=sqrt(ss.^2+ds.^2);
                rake=atan2(ds,ss)*180/pi;
                
                fname=sprintf('%s/receiver-%05d.flt',wdir,k);
                unicycle.export.exportflt_rfaults(...
                    s, ...
                    +scale*obj.rcv.x(:,2), ...
                    +scale*obj.rcv.x(:,1), ...
                    -scale*obj.rcv.x(:,3), ...
                    scale*obj.rcv.L, ...
                    scale*obj.rcv.W, ...
                    obj.rcv.strike, ...
                    obj.rcv.dip, ...
                    rake, ...
                    fname, ...
                    'time (yr)',t(k))
            else
                for k=1:length(t)
                    ss=y(1:obj.dgf:end,k);
                    ds=y(2:obj.dgf:end,k);
                    s=sqrt(ss.^2+ds.^2);
                    rake=atan2(ds,ss)*180/pi;
                    
                    fname=sprintf('%s/receiver-%05d.flt',wdir,k);
                    unicycle.export.exportflt_rfaults(...
                        s, ...
                        +scale*obj.rcv.x(:,2), ...
                        +scale*obj.rcv.x(:,1), ...
                        -scale*obj.rcv.x(:,3), ...
                        scale*obj.rcv.L, ...
                        scale*obj.rcv.W, ...
                        obj.rcv.strike, ...
                        obj.rcv.dip, ...
                        rake, ...
                        fname, ...
                        'time (yr)',t(k));
                end
            end
            
        end % end exportflt
        
        function []=exportxyz(obj,t,y,scale,wdir)
            % EXPORTXYZ write a series of .xyz files, one for each time step.
            %
            %   evl.exportxyz(evl.t,evl.y,1e-3,'xyz')
            %
            % input:
            % t         - time
            % y         - solution vector
            % scale     - spatial scale factor for fault position and dimension
            % wdir      - export directory
            %
            % SEE ALSO: unicycle
            
            assert(length(t)==size(y,2));
            
            [xp,yp,zp,~]=unicycle.geometry.transform4patch_general(...
                obj.rcv.x(:,1),obj.rcv.x(:,2),-obj.rcv.x(:,3),...
                obj.rcv.L*0,obj.rcv.L,obj.rcv.W,obj.rcv.dip,obj.rcv.strike);
            
            for k=1:length(t)
                ss=y(1:obj.dgf:end,k);
                ds=y(2:obj.dgf:end,k);
                s=sqrt(ss.^2+ds.^2);
                
                fname=sprintf('%s/receiver-%05d.xyz',wdir,k);
                unicycle.export.exportxyz_rfaults(...
                    s, ...
                    scale*xp, ...
                    scale*yp, ...
                    -scale*zp, ...
                    fname, ...
                    'time (yr)',t(k))
            end
            
        end % end exportxyz
        
        function [] = plotHorizontalProfiles(obj,depth,threshold,seismicPeriod,aseismicPeriod,component)
            % PLOTHORIZONTALPROFILES plot cumulative fault slip across
            % a horizontal profile at a constant fault depth. Slip occurring
            % at seismic and aseismic velocities is plotting differently.
            %
            % INPUT:
            %
            % depth          depth of the horizontal profile
            % threshold      critical velocity for seismic slip (m/s)
            % seismicPeriod  sample period for seismic slip
            % aseismicPeriod sample period for aseismic slip
            % component      displacement component (1:strike slip, 2: dip slip)
            %
            % SEE ALSO: unicycle
            
            s2y=60*60*24*365;
            y2s=1./s2y;
            
            offset=0;
            for k=1:length(obj.rcv.segments)
                % patch index
                pos=(1:obj.rcv.segments{k}.nPatch)+obj.rcv.segments{k}.starti;
                wMin=min(obj.rcv.W(pos).*sind(obj.rcv.dip(pos)));
                pos=find(abs(obj.rcv.xc(pos,3)+depth-wMin/4)<=wMin/2);
                
                % along-strike coordinates
                x=(obj.rcv.x(pos,:)-repmat(obj.rcv.segments{k}.x,length(pos),1))*obj.rcv.segments{k}.sv'+offset;
                
                % model index (cumulative strike- or dip slip)
                index=component+(pos-1)*obj.dgf;
                
                % time steps
                Dt=diff(obj.t);
                % velocity
                v=diff(obj.y(index,:))./repmat(Dt,1,length(pos))*y2s;
                
                % seismic and aseismic index
                index1=find(max(v,[],2)>=threshold);
                if ~isempty(index1)
                    % one-second seismic index
                    pos1=getPeriodicIndex(obj.t(index1),seismicPeriod);
                    if ~isempty(pos1)
                        plot(x,obj.y(index,index1(pos1)),'r')
                    end
                else
                    fprintf('no velocity above %8.2e\n',threshold);
                end
                
                index2=find(max(v,[],2) <threshold);
                if ~isempty(index2)
                    % one-year aseismic index
                    pos2=getPeriodicIndex(obj.t(index2),aseismicPeriod);
                    if ~isempty(pos2)
                        plot(x,obj.y(index,index2(pos2)),'b--')
                    end
                else
                    fprintf('no velocity below %8.2e\n',threshold);
                end
                
                offset=obj.rcv.segments{k}.L;
            end
            
        end
        
        function [] = plotIndexProfiles(obj,pos,threshold,seismicPeriod,aseismicPeriod,component,transpose)
            % PLOTINDEXPROFILES plot cumulative fault slip across
            % a horizontal profile at a constant fault depth. Slip occurring
            % at seismic and aseismic velocities is plotting differently.
            %
            % INPUT:
            %
            % pos            index of patches to plot
            % threshold      critical velocity for seismic slip (m/s)
            % seismicPeriod  sample period for seismic slip
            % aseismicPeriod sample period for aseismic slip
            % component      displacement component (1:strike slip, 2: dip slip)
            % transpose      true for x,y, false for y,x
            
            s2y=60*60*24*365;
            y2s=1./s2y;
            
            % along-direction coordinates
            n=(obj.rcv.x(pos(end),:)-obj.rcv.x(pos(1),:))';
            n=n/sqrt(sum(n.^2));
            x=obj.rcv.x(pos,:)*n;
            
            % model index (cumulative strike- or dip slip)
            index=component+(pos-1)*obj.dgf;
            
            % time steps
            Dt=diff(obj.t);
            % velocity
            v=diff(obj.y(index,:))./repmat(Dt,1,length(pos))*y2s;
            
            % seismic and aseismic index
            index1=find(max(v,[],2)>=threshold);
            if ~isempty(index1)
                % one-second seismic index
                pos1=getPeriodicIndex(obj.t(index1),seismicPeriod);
                if ~isempty(pos1)
                    if transpose
                        plot(obj.y(index,index1(pos1)),x,'r')
                    else
                        plot(x,obj.y(index,index1(pos1)),'r')
                    end
                end
            else
                fprintf('no velocity above %8.2e\n',threshold);
            end
            
            index2=find(max(v,[],2) <threshold);
            if ~isempty(index2)
                % one-year aseismic index
                pos2=getPeriodicIndex(obj.t(index2),aseismicPeriod);
                if ~isempty(pos2)
                    if transpose
                        plot(obj.y(index,index2(pos2)),x,'b--')
                    else
                        plot(x,obj.y(index,index2(pos2)),'b--')
                    end
                end
            else
                fprintf('no velocity below %8.2e\n',threshold);
            end
            
        end
        
        function [] = plotVerticalProfiles(obj,lo,threshold,seismicPeriod,aseismicPeriod,component)
            % PLOTVERTICALPROFILES plot cumulative fault slip down
            % a vertical profile at a constant along-strike distance.
            % Slip occurring at seismic and aseismic velocities is plotting
            % differently.
            %
            % INPUT:
            %
            % lo             along-strike coordinate of vertical profile
            % threshold      critical velocity for seismic slip (m/s)
            % seismicPeriod  sample period for seismic slip
            % aseismicPeriod sample period for aseismic slip
            % component      displacement component (1:strike slip, 2: dip slip)
            %
            % SEE ALSO: unicycle
            
            s2y=60*60*24*365;
            y2s=1./s2y;
            
            lMin=min(obj.rcv.L);
            
            offset=0;
            for k=1:length(obj.rcv.segments)
                % along-strike coordinates
                x=(obj.rcv.x-repmat(obj.rcv.segments{k}.x,obj.rcv.segments{k}.nPatch,1))*...
                    obj.rcv.segments{k}.sv'+offset;
                
                % patch index
                pos=find(abs(x-lo)<=lMin/2);
                
                % down-dip coordinates
                y=(obj.rcv.x(pos,:)-repmat(obj.rcv.segments{k}.x,length(pos),1))*...
                    obj.rcv.segments{k}.dv'+offset;
                
                % model index (cumulative strike- or dip slip)
                index=component+(pos-1)*obj.dgf;
                
                % time steps
                Dt=diff(obj.t);
                % velocity
                v=diff(obj.y(:,index))./repmat(Dt,1,length(pos))*y2s;
                
                % seismic and aseismic index
                index1=find(max(v,[],2)>=threshold);
                index2=find(max(v,[],2) <threshold);
                
                if ~isempty(index1)
                    % one-second (seismicPeriod) seismic index
                    pos1=getPeriodicIndex(obj.t(index1),seismicPeriod);
                    if ~isempty(pos1)
                        plot(obj.y(index1(pos1),index),y,'r')
                    end
                else
                    fprintf('no velocity above %8.2e\n',threshold);
                end
                
                if ~isempty(index2)
                    % one-year (aseismicPeriod) aseismic index
                    pos2=getPeriodicIndex(obj.t(index2),aseismicPeriod);
                    if ~isempty(pos2)
                        plot(obj.y(index2(pos2),index),y,'b--')
                    end
                else
                    fprintf('no velocity below %8.2e\n',threshold);
                end
                
                offset=obj.rcv.segments{k}.L;
            end
            
        end
        
        function t = eventstress(obj,evtIndex)
            % function EVENTSTRESS computes the stress change due to an externally
            % imposed slip-distribution event.
            %
            % INPUT:
            %   evtIndex        - the number of the event from the object's event list
            %
            % OUTPUT:
            %   t               - a minimum size state vector containing stress components
            %
            % SEE ALSO: unicycle
            
            import hmmvp.*
            
            
            ss=obj.evt{evtIndex}.src.slip.*cosd(obj.evt{evtIndex}.src.rake);
            ds=obj.evt{evtIndex}.src.slip.*sind(obj.evt{evtIndex}.src.rake);

            if isobject(obj.flt)
                Ks=obj.evt{evtIndex}.KK{1,1}*ss+obj.evt{evtIndex}.KK{2,1}*ds;
                Kd=obj.evt{evtIndex}.KK{1,2}*ss+obj.evt{evtIndex}.KK{2,2}*ds;
                Kn=obj.evt{evtIndex}.KK{1,3}*ss+obj.evt{evtIndex}.KK{2,3}*ds;
            else
                Ks=[];
                Kd=[];
                Kn=[];
            end
            if isobject(obj.shz)
                L11=obj.evt{evtIndex}.KL{1,1}*ss+obj.evt{evtIndex}.KL{2,1}*ds;
                L12=obj.evt{evtIndex}.KL{1,2}*ss+obj.evt{evtIndex}.KL{2,2}*ds;
                L13=obj.evt{evtIndex}.KL{1,3}*ss+obj.evt{evtIndex}.KL{2,3}*ds;
                L22=obj.evt{evtIndex}.KL{1,4}*ss+obj.evt{evtIndex}.KL{2,4}*ds;
                L23=obj.evt{evtIndex}.KL{1,5}*ss+obj.evt{evtIndex}.KL{2,5}*ds;
                L33=obj.evt{evtIndex}.KL{1,6}*ss+obj.evt{evtIndex}.KL{2,6}*ds;
            else
                L11=[];
                L12=[];
                L13=[];
                L22=[];
                L23=[];
                L33=[];
            end
            
            k=[Ks,Kd,Kn]';
            l=[L11,L12,L13,L22,L23,L33]';
            t=[k(:);l(:)];
        end
        
        function migrateTo(obj,c)
            % method MIGRATETO converts obj to new evolution object.
            %
            % EXAMPLE:
            %
            %    evl2=unicycle.ode.rateStrengtheningBurgers();
            %    evl.migrateTo(evl2);
            %
            % SEE ALSO: unicycle.ode.evolution
            
            % fault self stress
            c.KK=obj.KK;
            
            % stress on shear zones due to fault slip
            c.KL=obj.KL;
            
            % traction on fault due to shear zones
            c.LK=obj.LK;
            
            % shear zones self stress
            c.LL=obj.LL;
            
            % fault stress due to continuous source
            c.FK=obj.FK;
            
            % stress on shear zones due to continuous fault slip
            c.FL=obj.FL;
            
            % coseismic events with slip distribution and stress kernels
            c.evt=obj.evt;
            
            % source geometry and slip distribution
            c.src=obj.src;
            
            % fault receiver geometry
            dgf=c.flt.dgf;
            c.flt=obj.flt;
            c.flt.dgf=dgf;
            
            % fault receiver geometry
            dgf=c.shz.dgf;
            c.shz=obj.shz;
            c.shz.dgf=dgf;
            
            % simulation result (time and model)
            c.t=obj.t;
            c.y=obj.y;
            
        end
        
    end % end methods
    
    methods(Static)
        function evt=eventKernels(evts,flt,shz)
            
            % coseismic events with prescribed slip distribution and timing
            evt=cell(length(evts),1);
            
            % order events by chronological order
            temp=zeros(length(evts),1);
            for k=1:length(evts)
                temp(k)=evts{k}.t0;
            end
            [~,pos]=sort(temp);
            
            for k=1:length(evts)
                evt{k}=struct('src',[],'KK',[],'KL',[]);
                evt{k}.src=evts{pos(k)};
 
                
                % OLD CODE
%                  fnameKK=['./kernels/evt_' num2str(k) '_KK.mat'];
%                  fnameKL=['./kernels/evt_' num2str(k) '_KL.mat'];
%                 if isobject(flt)
%                     if ~exist(fnameKK,'file')
%                         evt{k}.KK=cell(2,3);
%                         [evt{k}.KK{1,1},evt{k}.KK{2,1}, ...
%                             evt{k}.KK{1,2},evt{k}.KK{2,2}, ...
%                             evt{k}.KK{1,3},evt{k}.KK{2,3}]=...
%                             evt{pos(k)}.src.tractionKernels(flt);
%                         sKK=evt{k}.KK;
%                         save(fnameKK,'sKK');
%                         clear sKK
%                     else
%                         load(fnameKK);
%                         evt{k}.KK=sKK;
%                         clear sKK
%                     end
%                 end
                
                % NEW CODE
                            %o.knl.KK={'ss', 'sd', 'sn'; 'ds', 'dd', 'dn'};
                fnameKK={['./kernels/evt_' num2str(k) '_KK_ss.grd'],...
                    ['./kernels/evt_' num2str(k) '_KK_sd.grd'],...
                    ['./kernels/evt_' num2str(k) '_KK_sn.grd'];...
                    ['./kernels/evt_' num2str(k) '_KK_ds.grd'],...
                    ['./kernels/evt_' num2str(k) '_KK_dd.grd'],...
                    ['./kernels/evt_' num2str(k) '_KK_dn.grd']};
                fnameKKmat=['./kernels/evt_' num2str(k) '_KK.mat'];
                if isobject(flt)
                    evt{k}.KK=cell(2,3);
                    if ~exist(fnameKK{1},'file')
                        if ~exist(fnameKKmat,'file')
                            [evt{k}.KK{:}]=evt{pos(k)}.src.tractionKernels(flt);
                            for i = 1:numel(evt{k}.KK)
                                unicycle.export.grdwrite([0 1], [0 1],evt{k}.KK{i},fnameKK{i})
                            end
                        else
                            load(fnameKKmat);
                            evt{k}.KK=sKK;
                            clear sKK
                            for i = 1:numel(evt{k}.KK)
                                unicycle.export.grdwrite([0 1], [0 1],evt{k}.KK{i},fnameKK{i})
                            end
                        end
                    else
                        for i = 1:numel(evt{k}.KK)
                            [~,~,evt{k}.KK{i}]=unicycle.export.grdread(fnameKK{i});
                        end
                    end
                end
                                    
                fnameKL={['./kernels/evt_' num2str(k) '_KL_s11.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_s12.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_s13.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_s22.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_s23.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_s33.grd'];...
                    ['./kernels/evt_' num2str(k) '_KL_d11.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_d12.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_d13.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_d22.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_d23.grd'],...
                    ['./kernels/evt_' num2str(k) '_KL_d33.grd']};
                fnameKLmat=['./kernels/evt_' num2str(k) '_KL.mat'];
                %          new code
                if isobject(shz)
                    evt{k}.KL=cell(2,6);
                    if ~exist(fnameKL{1},'file')
                        if ~exist(fnameKLmat,'file')
                            [evt{k}.KL{:}]=evt{pos(k)}.src.stressKernels(shz);
                            for i = 1:numel(evt{k}.KL)
                                unicycle.export.grdwrite([0 1], [0 1],evt{k}.KL{i},fnameKL{i})
                            end
                        else
                            load(fnameKLmat);
                            evt{k}.KL=sKL;
                            clear sKL
                            for i = 1:numel(evt{k}.KL)
                                unicycle.export.grdwrite([0 1], [0 1],evt{k}.KL{i},fnameKL{i})
                            end
                        end
                    else
                        for i = 1:numel(evt{k}.KL)
                            [~,~,evt{k}.KL{i}]=unicycle.export.grdread(fnameKL{i});
                        end
                    end
                end
                %
               % fnameKL=['./kernels/evt_' num2str(k) '_KL.mat'];
                
                %                Old code
%                 if isobject(shz)
%                     
%                     if ~exist(fnameKL,'file')
%                         evt{k}.KL=cell(2,6);
%                         [evt{k}.KL{1,1},evt{k}.KL{2,1}, ...
%                             evt{k}.KL{1,2},evt{k}.KL{2,2}, ...
%                             evt{k}.KL{1,3},evt{k}.KL{2,3}, ...
%                             evt{k}.KL{1,4},evt{k}.KL{2,4}, ...
%                             evt{k}.KL{1,5},evt{k}.KL{2,5}, ...
%                             evt{k}.KL{1,6},evt{k}.KL{2,6}]=...
%                             evt{pos(k)}.src.stressKernels(shz);
%                         sKL=evt{k}.KL;
%                         save(fnameKL,'sKL');
%                         clear sKL
%                     else
%                         load(fnameKL);
%                         evt{k}.KL=sKL;
%                         clear sKL
%                     end
%                 end
            end
        end
    end % end static methods
end

function [pos] = getPeriodicIndex(t,period)
% function GETPERIODICINDEX finds the position in array t that
% are the closest to a periodic.

to=t(1);
pos=[];
for k=1:length(t)
    tc=t(k);
    if (tc-to>period)
        to=tc;
    end
    if tc>=to
        pos=[pos;k];
        to=tc+period;
    end
end
end

