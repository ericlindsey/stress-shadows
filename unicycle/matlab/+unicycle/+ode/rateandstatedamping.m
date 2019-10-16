classdef rateandstatedamping < unicycle.ode.evolution
    properties
        % damping Coefficient
        dampingCoefficient;
        % current velocity, used as initial guess for Newtown-Raphson
        previousVelocity;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=rateandstatedamping(src,rcv,evts)
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
            % The evolution law comes about from the classic evolution law
            %
            %    d psi
            %    ----- = 1 - V psi / L
            %      dt
            %
            % with the change of variable
            %
            %    theta = log(Vo psi / L)
            %
            % or
            %
            %    psi = exp(theta) / Vo * L
            %
            % EXAMPLE:
            %
            %   evl = ode.rateandstate(src,rcv,evts,dampingCoefficient)
            %
            % creates a instance from source (src), receiver (rcv), events
            % (evt) and radiation damping coefficient.
            %
            % type methods(evl) and properties(evl) for a list of the class
            % methods and properties.
            
            obj=obj@unicycle.ode.evolution(src,rcv,[],evts);
            obj.flt.dgf=5;
            
        end % constructor
        
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
            
            import unicycle.utils.*
            
            % normal stress
            dsigma=obj.KK{1,3}*y(1:obj.flt.dgf:end)+obj.KK{2,3}*y(2:obj.flt.dgf:end);
            
            % shear stress in strike direction
            tss=y(3:obj.flt.dgf:end);
            
            % shear stress in dip direction
            tsd=y(4:obj.flt.dgf:end);
            
            % state variable
            th=y(5:obj.flt.dgf:end);
            
            % cumulative shear stress
            tau=sqrt(tss.^2+tsd.^2);
            
            % rake of shear stress and velocity
            rak=atan2(tsd,tss);
            
            % norm of velocity
            v=(2.*obj.flt.Vs.*obj.flt.a.*obj.flt.sigma./obj.flt.earthModel.G).*...
                lambertW(obj.flt.earthModel.G*obj.flt.Vo./(2*obj.flt.Vs.*obj.flt.a.*obj.flt.sigma).*...
                exp((tau-obj.flt.mu0.*obj.flt.sigma-obj.flt.sigma.*obj.flt.b.*th)./(obj.flt.sigma.*obj.flt.a)));
            
            % enforce rake constraint
            if obj.flt.isRakeConstraint
                v((tss.*cosd(obj.flt.rake)+tsd.*sind(obj.flt.rake))<=0)=0;
            end
            
            % velocity in strike direction
            dss=v.*cos(rak);
            
            % velocity in dip direction
            dds=v.*sin(rak);
            
            % initialize yp
            yp=zeros(size(y));
            
            % strike-slip velocity
            yp(1:obj.flt.dgf:end)=dss;
            
            % dip-slip velocity
            yp(2:obj.flt.dgf:end)=dds;
            
            % forcing term from plate velocity
            dss=dss-obj.flt.Vpl.*cosd(obj.flt.Vrake);
            dds=dds-obj.flt.Vpl.*sind(obj.flt.Vrake);
            
            % shear stress rate of change in strike direction
            yp(3:obj.flt.dgf:end)=obj.KK{1,1}*dss+obj.KK{2,1}*dds;
            
            % shear stress rate of change in dip direction
            yp(4:obj.flt.dgf:end)=obj.KK{1,2}*dss+obj.KK{2,2}*dds;
            
            % constant loading
            if 0<obj.src.N
                % shear stress rate of change in strike direction
                yp(3:obj.flt.dgf:end)=yp(3:obj.flt.dgf:end)+...
                    obj.FK{1,1}*(obj.src.slip.*cosd(obj.src.rake))+...
                    obj.FK{2,1}*(obj.src.slip.*sind(obj.src.rake));
                
                % shear stress rate of change in dip direction
                yp(4:obj.flt.dgf:end)=yp(4:obj.flt.dgf:end)+...
                    obj.FK{1,2}*(obj.src.slip.*cosd(obj.src.rake))+...
                    obj.FK{2,2}*(obj.src.slip.*sind(obj.src.rake));
            end
            
            % rate of state
            yp(5:obj.flt.dgf:end)=(obj.flt.Vo.*exp(-th)-v)./obj.flt.l;
            
        end % end ratestate
        
    end % end methods
    
end
