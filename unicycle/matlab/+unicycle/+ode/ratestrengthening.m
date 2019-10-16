classdef ratestrengthening < unicycle.ode.evolution
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=ratestrengthening(varargin)
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
            %   evl = ode.ratestrengthening(src,rcv,event);
            %
            % creates a instance from source (src), receiver (rcv) and
            % events.
            %
            % SEE ALSO: unicycle
            
            obj=obj@unicycle.ode.evolution(varargin{:});
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
            
            % shear stress in strike direction
            tss=y(3:obj.flt.dgf:end);
            
            % shear stress in dip direction
            tsd=y(4:obj.flt.dgf:end);
            
            
            % rake is fixed by user at the patches flagged by obj.flt.FixedRakePosition
            proj=tss(obj.flt.FixedRakePosition).*cosd(obj.flt.rake(obj.flt.FixedRakePosition))+...
                 tsd(obj.flt.FixedRakePosition).*sind(obj.flt.rake(obj.flt.FixedRakePosition));
            tss(obj.flt.FixedRakePosition)=proj.*cosd(obj.flt.rake(obj.flt.FixedRakePosition));
            tsd(obj.flt.FixedRakePosition)=proj.*sind(obj.flt.rake(obj.flt.FixedRakePosition));
            
            % cumulative shear stress
            tau=sqrt(tss.^2+tsd.^2);
            
            % rake of shear stress and velocity
            rak=atan2(tsd,tss);
            
            % enforce rake constraint
            if obj.flt.isRakeConstraint
                tau((tss.*cosd(obj.flt.rake)+tsd.*sind(obj.flt.rake))<=0)=0;
            end
            
            % pin patches
            if 0~=numel(obj.flt.pinnedPosition)
                tau(obj.flt.pinnedPosition)=0;
            end
            
            % norm of velocity
            v=obj.v(tau);
            
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

            % shear stress rate of change in strike direction
            yp(3:obj.flt.dgf:end)=obj.KK{1,1}*dss+obj.KK{2,1}*dds;
            
            % shear stress rate of change in dip direction
            yp(4:obj.flt.dgf:end)=obj.KK{1,2}*dss+obj.KK{2,2}*dds;
            
            % constant loading
            if isobject(obj.src)
                if 0<obj.src.N
                    % shear stress rate of change in strike direction
                    yp(3:obj.flt.dgf:end)=yp(3:obj.flt.dgf:end)+...
                        obj.Fss*(obj.src.slip.*cosd(obj.src.rake))+...
                        obj.Fds*(obj.src.slip.*sind(obj.src.rake));

                    % shear stress rate of change in dip direction
                    yp(4:obj.flt.dgf:end)=yp(4:obj.flt.dgf:end)+...
                        obj.Fsd*(obj.src.slip.*cosd(obj.src.rake))+...
                        obj.Fdd*(obj.src.slip.*sind(obj.src.rake));
                end
            end 
            
        end % end ratestrengthening
        
    end % end methods
end


