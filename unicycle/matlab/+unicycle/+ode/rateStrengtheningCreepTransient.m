classdef rateStrengtheningCreepTransient < unicycle.ode.evolution
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=rateStrengtheningCreepTransient(varargin)
            % RATESTRENGTHENINGBurgers is a class encapsulating the geometry 
            % and physical properties of fault patches for fault and shear
            % zones under rate-strengthening friction and a linear
            % viscoelastic rheology.
            %
            % The rate-strengthening friction relation
            %
            %       ds
            %   V = -- = 2*Vo*sinh(tau/(a*sigma))
            %       dt
            %
            % controls fault slip evolution. And the distributed flow is
            % controlled by
            %
            %   d epsilon   tau
            %   --------- = ---
            %       dt      eta
            %
            % EXAMPLE:
            %
            %   evl = ode.rateStrengtheningPower(src,flt,shz,evt);
            %
            % creates a instance from source (src), receiver (flt) and
            % events (evt).
            %
            % type methods(evl) and properties(evl) for a list of the class
            % methods and properties.
            
            obj=obj@unicycle.ode.evolution(varargin{:});
            
            obj.flt.dgf=4;
            obj.shz.dgf=18; 
            
        end % constructor
        
        function yp=ode(o,t,y)
            % ODE is a method of class RATESTRENGTHENINGPOWER used to
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
            %   V = -- = 2*Vo*sinh(tau/(a*sigma))
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
            % The distributed flow is controlled by
            %
            %    dE    dEM   dEK
            %    -- =  --- + ---
            %    dt     dt    dt
            %
            % with
            %
            %    dEK   (||sigma' - 2*G EK||)^(m-1)*(sigma' - 2*G EK)
            %    --- = ---------------------------------------------
            %    dt                         etaK
            %
            %    dEM   ||sigma'||^(n-1)*sigma'
            %    --- = -----------------------
            %    dt             etaM
            %
            % where etaM is the linear Maxwell viscosity, etaK is the
            % linear Kelvin viscosity and (m,n) are the power exponents.
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            %
            %                     F A U L T S
            %
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            
            % initialize yp
            yp=zeros(size(y));
            
            % normal stress
            %dsigma=o.KK{1,3}*y(1:o.flt.dgf:(o.flt.N*o.flt.dgf))+o.KK{2,3}*y(2:o.flt.dgf:(o.flt.N*o.flt.dgf));
            
            % shear stress in strike direction
            tss=y(3:o.flt.dgf:(o.flt.N*o.flt.dgf));
            
            % shear stress in dip direction
            tsd=y(4:o.flt.dgf:(o.flt.N*o.flt.dgf));
            
            % cumulative shear stress
            tau=sqrt(tss.^2+tsd.^2);
            
            % rake of shear stress and velocity
            rak=atan2(tsd,tss);
            
            % enforce rake constraint
            if o.flt.isRakeConstraint
                tau((tss.*cosd(o.flt.Vrake)+tsd.*sind(o.flt.Vrake))<=0)=0;
            end
            
            % pin patches
            if 0~=numel(o.flt.pinnedPosition)
                tau(o.flt.pinnedPosition)=0;
            end
            
            % norm of velocity
            v=2*o.flt.Vo.*sinh(tau./(o.flt.sigma.*o.flt.a));
            
            % velocity in strike direction
            dss=v.*cos(rak);
            
            % velocity in dip direction
            dds=v.*sin(rak);
            
            % strike-slip velocity
            yp(1:o.flt.dgf:(o.flt.N*o.flt.dgf))=dss;
            
            % dip-slip velocity
            yp(2:o.flt.dgf:(o.flt.N*o.flt.dgf))=dds;
            
            % forcing term from plate velocity
            %dss=dss-o.flt.Vpl.*cosd(o.flt.Vrake);
            dds=dds-o.flt.Vpl.*sind(o.flt.Vrake);
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            %
            %                 S H E A R   Z O N E S
            %
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
            % isotropic stress
            p=(y(o.flt.N*o.flt.dgf+7:o.shz.dgf:end)+y(o.flt.N*o.flt.dgf+10:o.shz.dgf:end)+y(o.flt.N*o.flt.dgf+12:o.shz.dgf:end))/3;
            
            % deviatoric stress components
            s11p=y(o.flt.N*o.flt.dgf+ 7:o.shz.dgf:end)-p;
            s12p=y(o.flt.N*o.flt.dgf+ 8:o.shz.dgf:end);
            s13p=y(o.flt.N*o.flt.dgf+ 9:o.shz.dgf:end);
            s22p=y(o.flt.N*o.flt.dgf+10:o.shz.dgf:end)-p;
            s23p=y(o.flt.N*o.flt.dgf+11:o.shz.dgf:end);
            s33p=y(o.flt.N*o.flt.dgf+12:o.shz.dgf:end)-p;
            
            
            % Kelvin self-stress tensor
            Q11k=s11p-o.shz.Gk.*y(o.flt.N*o.flt.dgf+13:o.shz.dgf:end);
            Q12k=s12p-o.shz.Gk.*y(o.flt.N*o.flt.dgf+14:o.shz.dgf:end);
            Q13k=s13p-o.shz.Gk.*y(o.flt.N*o.flt.dgf+15:o.shz.dgf:end);
            Q22k=s22p-o.shz.Gk.*y(o.flt.N*o.flt.dgf+16:o.shz.dgf:end);
            Q23k=s23p-o.shz.Gk.*y(o.flt.N*o.flt.dgf+17:o.shz.dgf:end);
            Q33k=s33p-o.shz.Gk.*y(o.flt.N*o.flt.dgf+18:o.shz.dgf:end);
            
            % Frobenius norm of the self-stress tensor
            QL2=sqrt(Q11k.^2+2*Q12k.^2+2*Q13k.^2+Q22k.^2+2*Q23k.^2+Q33k.^2)/sqrt(2);
            
            % Kelvin strain rates d epsilonK / dt = Ak*(sigma' - 2G epsilonK)^m
            e11k=Q11k.*(QL2.^(o.shz.m-1)).*o.shz.Ak;
            e12k=Q12k.*(QL2.^(o.shz.m-1)).*o.shz.Ak;
            e13k=Q13k.*(QL2.^(o.shz.m-1)).*o.shz.Ak;
            e22k=Q22k.*(QL2.^(o.shz.m-1)).*o.shz.Ak;
            e23k=Q23k.*(QL2.^(o.shz.m-1)).*o.shz.Ak;
            e33k=Q33k.*(QL2.^(o.shz.m-1)).*o.shz.Ak;
            
            % Frobenius norm of the shear stress tensor
            tauL2=sqrt(s11p.^2+2*s12p.^2+2*s13p.^2+s22p.^2+2*s23p.^2+s33p.^2)/sqrt(2);

            % Total strain = Kelvin + Maxwell strain rates
            e11=e11k+s11p.*(tauL2.^(o.shz.n-1)).*o.shz.A;
            e12=e12k+s12p.*(tauL2.^(o.shz.n-1)).*o.shz.A;
            e13=e13k+s13p.*(tauL2.^(o.shz.n-1)).*o.shz.A;
            e22=e22k+s22p.*(tauL2.^(o.shz.n-1)).*o.shz.A;
            e23=e23k+s23p.*(tauL2.^(o.shz.n-1)).*o.shz.A;
            e33=e33k+s33p.*(tauL2.^(o.shz.n-1)).*o.shz.A;
            
            % rate of strain
            yp(o.flt.N*o.flt.dgf+1:o.shz.dgf:end)=e11;
            yp(o.flt.N*o.flt.dgf+2:o.shz.dgf:end)=e12;
            yp(o.flt.N*o.flt.dgf+3:o.shz.dgf:end)=e13;
            yp(o.flt.N*o.flt.dgf+4:o.shz.dgf:end)=e22;
            yp(o.flt.N*o.flt.dgf+5:o.shz.dgf:end)=e23;
            yp(o.flt.N*o.flt.dgf+6:o.shz.dgf:end)=e33;
            
            % rate of Kelvin strain
            yp(o.flt.N*o.flt.dgf+13:o.shz.dgf:end)=e11k;
            yp(o.flt.N*o.flt.dgf+14:o.shz.dgf:end)=e12k;
            yp(o.flt.N*o.flt.dgf+15:o.shz.dgf:end)=e13k;
            yp(o.flt.N*o.flt.dgf+16:o.shz.dgf:end)=e22k;
            yp(o.flt.N*o.flt.dgf+17:o.shz.dgf:end)=e23k;
            yp(o.flt.N*o.flt.dgf+18:o.shz.dgf:end)=e33k;
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %
            %
             %          S T R E S S   I N T E R A C T I O N S
            %
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            % shear stress rate of change in strike direction
            yp(3:o.flt.dgf:(o.flt.N*o.flt.dgf))= ...
                o.KK{1,1}*dss+o.KK{2,1}*dds ...
               +o.LK{1,1}*e11+o.LK{2,1}*e12+o.LK{3,1}*e13+o.LK{4,1}*e22+o.LK{5,1}*e23+o.LK{6,1}*e33;
            
            % shear stress rate of change in dip direction
            yp(4:o.flt.dgf:(o.flt.N*o.flt.dgf))= ...
                o.KK{1,2}*dss+o.KK{2,2}*dds ...
               +o.LK{1,2}*e11+o.LK{2,2}*e12+o.LK{3,2}*e13+o.LK{4,2}*e22+o.LK{5,2}*e23+o.LK{6,2}*e33; 
            
            % constant loading
            if isobject(o.src)
                if 0<o.src.N
                    % shear stress rate of change in strike direction
                    yp(3:o.dgf:(o.flt.N*o.flt.dgf))=yp(3:o.flt.dgf:(o.flt.N*o.flt.dgf))+...
                        o.FK{1,1}*(o.src.slip.*cosd(o.src.rake))+...
                        o.FK{2,1}*(o.src.slip.*sind(o.src.rake));
                    
                    % shear stress rate of change in dip direction
                    yp(4:o.dgf:(o.flt.N*o.flt.dgf))=yp(4:o.flt.dgf:(o.flt.N*o.flt.dgf))+...
                        o.FK{1,2}*(o.src.slip.*cosd(o.src.rake))+...
                        o.FK{2,2}*(o.src.slip.*sind(o.src.rake));
                end
            end
            
            % shear stress rate of change
            yp(o.flt.N*o.flt.dgf+ 7:o.shz.dgf:end)=o.LL{1,1}*e11+o.LL{2,1}*e12+o.LL{3,1}*e13+o.LL{4,1}*e22+o.LL{5,1}*e23+o.LL{6,1}*e33+o.KL{1,1}*dss+o.KL{2,1}*dds;
            yp(o.flt.N*o.flt.dgf+ 8:o.shz.dgf:end)=o.LL{1,2}*e11+o.LL{2,2}*e12+o.LL{3,2}*e13+o.LL{4,2}*e22+o.LL{5,2}*e23+o.LL{6,2}*e33+o.KL{1,2}*dss+o.KL{2,2}*dds;
            yp(o.flt.N*o.flt.dgf+ 9:o.shz.dgf:end)=o.LL{1,3}*e11+o.LL{2,3}*e12+o.LL{3,3}*e13+o.LL{4,3}*e22+o.LL{5,3}*e23+o.LL{6,3}*e33+o.KL{1,3}*dss+o.KL{2,3}*dds;
            yp(o.flt.N*o.flt.dgf+10:o.shz.dgf:end)=o.LL{1,4}*e11+o.LL{2,4}*e12+o.LL{3,4}*e13+o.LL{4,4}*e22+o.LL{5,4}*e23+o.LL{6,4}*e33+o.KL{1,4}*dss+o.KL{2,4}*dds;
            yp(o.flt.N*o.flt.dgf+11:o.shz.dgf:end)=o.LL{1,5}*e11+o.LL{2,5}*e12+o.LL{3,5}*e13+o.LL{4,5}*e22+o.LL{5,5}*e23+o.LL{6,5}*e33+o.KL{1,5}*dss+o.KL{2,5}*dds;
            yp(o.flt.N*o.flt.dgf+12:o.shz.dgf:end)=o.LL{1,6}*e11+o.LL{2,6}*e12+o.LL{3,6}*e13+o.LL{4,6}*e22+o.LL{5,6}*e23+o.LL{6,6}*e33+o.KL{1,6}*dss+o.KL{2,6}*dds;
            

            
        end % end ratestate
        

        
    end % end methods
    
end
