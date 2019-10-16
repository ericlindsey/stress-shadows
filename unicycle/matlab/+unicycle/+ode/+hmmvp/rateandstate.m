classdef rateandstate < unicycle.ode.evolution
    properties
        % tolerance for hierarchical matrix
        tol;
        % number of threads for matrix init and matrix-vector multiplication
        nthreads;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=rateandstate(src,rcv,evts,G,nu)
            % RATEANDSTATE is a class encapsulating the geometry and physical
            % properties of fault patches for sources and receivers and
            % methods for solving integral equations for fault slip using
            % rate- and state-dependent friction with the Aging Law.
            %
            % The rate-and-state friction constitutive:
            %
            %       ds
            %   V = -- = 2*Vo*exp(-b/a*theta)*sinh(tau/(a*sigma))
            %       dt
            %
            % and
            %
            %   d theta
            %   ------- = (Vo*exp(-theta)-V)/L
            %      dt
            %
            % controls fault slip evolution. The governing equation is
            % static equilibrium (sum of forces equals zero).
            %
            %   evl = ode.rateandstate(src,rcv,G,nu);
            %
            % creates a instance from source (src), receiver (rcv) and
            % elastic properties (rigidity G, and Poisson's ratio nu).
            %
            % type methods(evl) and properties(evl) for a list of the class
            % methods and properties.
            
            import hmmvp.*

            obj.dgf=5;
            
            obj.G=G;
            obj.nu=nu;
            
            % receiver faults
            obj.rcv=rcv;
            % source faults
            obj.src=src;
            
            % initialize hierarchical matrices
            obj.tol=1e-5;
            obj.nthreads=12;
            
            base_fn='hmmvp/srcrcv/hm_';
            h=unihmmvp('build_hmatrices',src,rcv,G,nu,base_fn,obj.tol,obj.nthreads);
            
            % receiver stress interactions
            obj.Kss=h.Kss;
            obj.Ksd=h.Ksd;
            obj.Kds=h.Kds;
            obj.Kdd=h.Kdd;
            
            [obj.Kss.id,obj.Kss.nnz]=hmmvp('init',obj.Kss.hm_filename,obj.nthreads);
            [obj.Ksd.id,obj.Ksd.nnz]=hmmvp('init',obj.Ksd.hm_filename,obj.nthreads);
            [obj.Kds.id,obj.Kds.nnz]=hmmvp('init',obj.Kds.hm_filename,obj.nthreads);
            [obj.Kdd.id,obj.Kdd.nnz]=hmmvp('init',obj.Kdd.hm_filename,obj.nthreads);
            
            % source/receiver stress interactions
            if 1>obj.src.N
                obj.Fss=[];
                obj.Fsd=[];
                obj.Fds=[];
                obj.Fdd=[];
            else
                obj.Fss=h.Fss;
                obj.Fsd=h.Fsd;
                obj.Fds=h.Fds;
                obj.Fdd=h.Fdd;

                [obj.Fss.id,obj.Fss.nnz]=hmmvp('init',obj.Fss.hm_filename,obj.nthreads);
                [obj.Fsd.id,obj.Fsd.nnz]=hmmvp('init',obj.Fsd.hm_filename,obj.nthreads);
                [obj.Fds.id,obj.Fds.nnz]=hmmvp('init',obj.Fds.hm_filename,obj.nthreads);
                [obj.Fdd.id,obj.Fdd.nnz]=hmmvp('init',obj.Fdd.hm_filename,obj.nthreads);
            end
            
            % coseismic events with prescribed slip distribution and timing
            
            % order events by chronological order
            temp=zeros(length(evts),1);
            for k=1:length(evts)
                temp(k)=evts{k}.t0;
            end
            [~,pos]=sort(temp);
            
            for k=1:length(evts)
                obj.evt{k}=struct('src',[],'Fss',[],'Fsd',[],'Fds',[],'Fdd',[]);
                obj.evt{k}.src=evts{pos(k)};
                
                % initialize hierarchical matrices
                base_fn=sprintf('hmmvp/events/hm_%d_',k);
                h=unihmmvp('build_hmatrices',obj.evt{k}.src,rcv,G,nu,base_fn,obj.tol,obj.nthreads);
                
                % receiver stress interactions
                obj.evt{k}.Fss=h.Fss;
                obj.evt{k}.Fds=h.Fds;
                obj.evt{k}.Fsd=h.Fsd;
                obj.evt{k}.Fdd=h.Fdd;
                
                [obj.evt{k}.Fss.id,obj.evt{k}.Fss.nnz]=hmmvp('init',obj.evt{k}.Fss.hm_filename,obj.nthreads);
                [obj.evt{k}.Fds.id,obj.evt{k}.Fds.nnz]=hmmvp('init',obj.evt{k}.Fds.hm_filename,obj.nthreads);
                [obj.evt{k}.Fsd.id,obj.evt{k}.Fsd.nnz]=hmmvp('init',obj.evt{k}.Fsd.hm_filename,obj.nthreads);
                [obj.evt{k}.Fdd.id,obj.evt{k}.Fdd.nnz]=hmmvp('init',obj.evt{k}.Fdd.hm_filename,obj.nthreads);
                
            end
            
        end % constructor
        
        function velocity=v(obj,tau,theta)
            % function V computes the absolute value of velocity from
            % stress and state variable based on the governing equation.
            
            velocity=2*obj.rcv.Vo.*exp(-obj.rcv.b./obj.rcv.a.*theta).* ...
                sinh(tau./(obj.rcv.a.*obj.rcv.sigma));
        end
        
        function yp=ode(obj,t,y)
            % ODE is a method of class RATEANDSTATE that can be used to
            % solve the quasi-static equations for fault slip evolution
            % when slip is governed by rate-and-state friction:
            %
            %   [t,y]=ode.ode45(@evl.ratestate,[0 200],y0,options);
            %
            % where evl is an instance of class evolution.
            %
            % the rate-and-state friction consistutive equations are:
            %
            %       ds
            %   V = -- = 2*Vo*exp(-b/a*theta)*sinh(tau/(a*sigma))
            %       dt
            %
            % and
            %
            %   d theta
            %   ------- = (Vo*exp(-theta)-V)/L
            %      dt
            %
            % where a, b, Vo, and L are friction properties, and theta
            % is the state variable that describes the healing and weakening
            % of the fault contact.
            %
            
            import hmmvp.*
            
            % shear stress in strike direction
            tss=y(3:obj.dgf:end);
            
            % shear stress in dip direction
            tsd=y(4:obj.dgf:end);
            
            % state variable
            th=y(5:obj.dgf:end);
            
            % cumulative shear stress
            tau=sqrt(tss.^2+tsd.^2);
            
            % rake of shear stress and velocity
            rak=atan2(tsd,tss);
            
            % norm of velocity
            v=obj.v(tau,th);
            
            % enforce rake constraint
            if obj.rcv.isRakeConstraint
                mask=(tss.*cosd(obj.rcv.rake)+tsd.*sind(obj.rcv.rake))>=0;
            else
                mask=1;
            end
            
            % velocity in strike direction
            dss=v.*cos(rak).*mask;
            
            % velocity in dip direction
            dds=v.*sin(rak).*mask;
            
            % initialize yp
            yp=zeros(size(y));
            
            % strike-slip velocity
            yp(1:obj.dgf:end)=dss;
            
            % dip-slip velocity
            yp(2:obj.dgf:end)=dds;
            
            % shear stress rate of change in strike direction
            yp(3:obj.dgf:end)=hmmvp('mvp',obj.Kss.id,dss)+ ...
                              hmmvp('mvp',obj.Kds.id,dds);
            
            % shear stress rate of change in dip direction
            yp(4:obj.dgf:end)=hmmvp('mvp',obj.Ksd.id,dss)+ ...
                              hmmvp('mvp',obj.Kdd.id,dds);
            
            % constant loading
            if 0<obj.src.N
                % shear stress rate of change in strike direction
                yp(3:obj.dgf:end)=yp(3:obj.dgf:end)+...
                    hmmvp('mvp',obj.Fss.id,obj.src.slip.*cosd(obj.src.rake))+...
                    hmmvp('mvp',obj.Fds.id,obj.src.slip.*sind(obj.src.rake));
                
                % shear stress rate of change in dip direction
                yp(4:obj.dgf:end)=yp(4:obj.dgf:end)+...
                    hmmvp('mvp',obj.Fsd.id,obj.src.slip.*cosd(obj.src.rake))+...
                    hmmvp('mvp',obj.Fdd.id,obj.src.slip.*sind(obj.src.rake));
            end
            
            % rate of state
            yp(5:obj.dgf:end)=(obj.rcv.Vo.*exp(-th)-v)./obj.rcv.l;
            
        end % end ratestate
        
        function t = eventstress(obj,evtIndex,directionString)
            % function EVENTSTRESS computes the stress change due to an externally
            % imposed slip-distribution event.
            %
            % INPUT:
            %   evtIndex        - the number of the event from the object's event list
            %   directionString - the direction of coseismic stress change ('S' for
            %                     stress in strike-slip direction, 'D' for dip-slip
            %                     direction).
            
            import hmmvp.*
            
            if strcmp('S',directionString)
                t = hmmvp('mvp',obj.evt{evtIndex}.Fss.id,obj.evt{evtIndex}.src.slip.*cosd(obj.evt{evtIndex}.src.rake)) +...
                    hmmvp('mvp',obj.evt{evtIndex}.Fds.id,obj.evt{evtIndex}.src.slip.*sind(obj.evt{evtIndex}.src.rake));
            elseif strcmp('D',directionString)
                t = hmmvp('mvp',obj.evt{evtIndex}.Fsd.id,obj.evt{evtIndex}.src.slip.*cosd(obj.evt{evtIndex}.src.rake)) +...
                    hmmvp('mvp',obj.evt{evtIndex}.Fdd.id,obj.evt{evtIndex}.src.slip.*sind(obj.evt{evtIndex}.src.rake));
            end
        end
        
        function test_performance(obj)
            
            import hmmvp.*
            
            m=hmmvp('getm',obj.Kss.id);
            n=hmmvp('getn',obj.Kss.id);
            cf=m*n/obj.Kss.nnz;
            fprintf('compression factor: %f\n',cf);
            
            dss=rand(obj.rcv.N,1);
            tic
            for k=1:2000
                dum=hmmvp('mvp',obj.Kss.id,dss);
            end
            toc
            
        end
        
    end % end private methods
end


