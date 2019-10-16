classdef ratestrengthening_prestress < unicycle.ode.evolution
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=ratestrengthening_prestress(src,flt,evts)
            % RATESTRENGTHENING is a class encapsulating the geometry and physical
            % properties of fault patches for sources and receivers and
            % methods for solving integral equations for fault slip using
            % rate-dependent friction.
            %
            % The rate-and-state friction constitutive law
            %
            %       ds
            %   V = -- = 2*Vo*sinh(tau/((a-b)*sigma))
            %       dt
            %
            % controls fault slip evolution. The governing equation is
            % static equilibrium (sum of forces equals zero).
            %
            %   evl = ode.ratestrengthening_prestress(src,rcv,event,G,nu);
            %
            % creates a instance from source (src), receiver (rcv) and
            % elastic properties (rigidity G, and Poisson's ratio nu).
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            obj=obj@unicycle.ode.evolution(src,flt,[],evts);
            obj.flt.dgf=4;
            
        end % constructor
        
        function velocity=v(obj,tau)
            % function V computes the absolute value of velocity from
            % stress and state variable based on the governing equation.
            
            velocity=2*obj.flt.Vo.*sinh(tau./((obj.flt.a-obj.flt.b).*obj.flt.sigma));

        end
        
        function yp=ode(obj,t,y)
            % ODE is a method of class RATESTRENGTHENING that can be used to
            % solve the quasi-static equations for fault slip evolution
            % when slip is governed by rate-dependent friction:
            %
            %   [t,y]=ode.ode45(@evl.ratestate,[0 200],y0,options);
            %
            % where evl is an instance of class evolution.
            %
            % the rate-strengthening friction consistutive equation is:
            %
            %       ds
            %   V = -- = 2*Vo*sinh(tau/((a-b)*sigma))
            %       dt
            %
            % where a, b, and Vo are friction properties.
            %
            
            % afterslip shear stress in strike direction
            tss=y(3:obj.dgf:end);
            
            % afterslip shear stress in dip direction
            tsd=y(4:obj.dgf:end);
            
            % rake is fixed by user at the patches flagged by obj.flt.FixedRakePosition
            proj=tss(obj.flt.FixedRakePosition).*cosd(obj.flt.rake(obj.flt.FixedRakePosition))+...
                 tsd(obj.flt.FixedRakePosition).*sind(obj.flt.rake(obj.flt.FixedRakePosition));
            tss(obj.flt.FixedRakePosition)=proj.*cosd(obj.flt.rake(obj.flt.FixedRakePosition));
            tsd(obj.flt.FixedRakePosition)=proj.*sind(obj.flt.rake(obj.flt.FixedRakePosition));

            % enforce rake constraint for afterslip
            if obj.flt.isRakeConstraint
                ind=(tss.*cosd(obj.flt.rake)+tsd.*sind(obj.flt.rake))<0;
                tss(ind)=0;
                tsd(ind)=0;
            end

            % background stress
            cff0=obj.flt.tau-obj.flt.mu0.*obj.flt.sigma;
            tss=tss+cff0.*cosd(obj.flt.rake);
            tsd=tsd+cff0.*sind(obj.flt.rake);
            
            % cumulative shear stress
            tau=sqrt(tss.^2+tsd.^2);
            
            % rake of shear stress and velocity
            rak=atan2(tsd,tss);
            
            % norm of velocity
            v=obj.v(tau);
            
            % enforce rake constraint for afterslip + background stress
            if obj.flt.isRakeConstraint
                v((tss.*cosd(obj.flt.Vrake)+tsd.*sind(obj.flt.Vrake))<=0)=0;
            end

            % velocity in strike direction
            dss=v.*cos(rak);
            
            % velocity in dip direction
            dds=v.*sin(rak);
            
            % initialize yp
            yp=zeros(size(y));
            
            % strike-slip velocity
            yp(1:obj.dgf:end)=dss;
            
            % dip-slip velocity
            yp(2:obj.dgf:end)=dds;

            % shear stress rate of change in strike direction
            yp(3:obj.dgf:end)=obj.Kss*dss+obj.Kds*dds;
            
            % shear stress rate of change in dip direction
            yp(4:obj.dgf:end)=obj.Ksd*dss+obj.Kdd*dds;
            
            % subtract shear prestress rate of change in strike direction
            yp(3:obj.dgf:end)=yp(3:obj.dgf:end)-...
                    obj.Kss*(obj.flt.Vpl.*cosd(obj.flt.rake))-...
                    obj.Kds*(obj.flt.Vpl.*sind(obj.flt.rake));
                
            % subtract shear prestress rate of change in dip direction
            yp(4:obj.dgf:end)=yp(4:obj.dgf:end)-...
                    obj.Ksd*(obj.flt.Vpl.*cosd(obj.flt.rake))-...
                    obj.Kdd*(obj.flt.Vpl.*sind(obj.flt.rake));

            % constant loading
            if 0<obj.src.N
                % shear stress rate of change in strike direction
                yp(3:obj.dgf:end)=yp(3:obj.dgf:end)+...
                    obj.Fss*(obj.src.slip.*cosd(obj.src.rake))+...
                    obj.Fds*(obj.src.slip.*sind(obj.src.rake));
                
                % shear stress rate of change in dip direction
                yp(4:obj.dgf:end)=yp(4:obj.dgf:end)+...
                    obj.Fsd*(obj.src.slip.*cosd(obj.src.rake))+...
                    obj.Fdd*(obj.src.slip.*sind(obj.src.rake));
            end
            
        end % end ratestrengthening_prestress
        
    end % end methods
end


