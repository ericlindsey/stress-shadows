classdef burgers < unicycle.ode.evolution
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=burgers(src,shz,evts)
            % BURGERS is a class encapsulating the geometry and physical 
            % properties of shear zones with a linear Burgers viscoelastic 
            % rheology.
            %
            % The distributed flow is controlled by
            %
            %    dE    dEM   dEK
            %    -- =  --- + ---
            %    dt     dt    dt
            %
            % with
            %
            %    dEK   sigma' - 2*G EK
            %    --- = ---------------
            %    dt         etaK
            %
            %    dEM   sigma'
            %    --- = ------
            %    dt     etaM
            %
            % where etaM is the linear Maxwell viscosity and etaK is the
            % linear Kelvin viscosity.
            %
            % EXAMPLE:
            %
            %   evl = ode.burgers(src,shz,evt);
            %
            % creates a instance from source (src), shear zone (shz) and
            % events (evt).
            %
            % type methods(evl) and properties(evl) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle
            
            obj=obj@unicycle.ode.evolution(src,[],shz,evts);
            obj.shz.dgf=18;
            
        end % constructor
        
        function yp=ode(o,t,y)
            % ODE is a method of class MAXWELL that can be used to
            % solve the quasi-static equations for viscoelastic flow
            % with a Maxwell (linear) rheology
            %
            %   ode.ode45(@evl.maxwell,[0 200],y0,options);
            %
            % where evl is an instance of class evolution, creates evl.t
            % and evl.y.
            %
            %    dE    dEM   dEK
            %    -- =  --- + ---
            %    dt     dt    dt
            %
            % with
            %
            %    dEK   sigma' - 2*G EK
            %    --- = ---------------
            %    dt         etaK
            %
            %    dEM   sigma'
            %    --- = ------
            %    dt     etaM
            %
            % where etaM is the linear Maxwell viscosity and etaK is the
            % linear Kelvin viscosity.
            
            % isotropic stress
            p=(y(7:o.shz.dgf:end)+y(10:o.shz.dgf:end)+y(12:o.shz.dgf:end))/3;
            
            % deviatoric stress components
            s11p=y( 7:o.shz.dgf:end)-p;
            s12p=y( 8:o.shz.dgf:end);
            s13p=y( 9:o.shz.dgf:end);
            s22p=y(10:o.shz.dgf:end)-p;
            s23p=y(11:o.shz.dgf:end);
            s33p=y(12:o.shz.dgf:end)-p;
            
            % initialize yp
            yp=zeros(size(y));
            
            % Kelvin strain rates d epsilonK / dt = ( sigma' - 2G epsilonK ) / etaK
            e11k=(s11p-o.shz.Gk.*y(13:o.shz.dgf:end))./o.shz.etaK;
            e12k=(s12p-o.shz.Gk.*y(14:o.shz.dgf:end))./o.shz.etaK;
            e13k=(s13p-o.shz.Gk.*y(15:o.shz.dgf:end))./o.shz.etaK;
            e22k=(s22p-o.shz.Gk.*y(16:o.shz.dgf:end))./o.shz.etaK;
            e23k=(s23p-o.shz.Gk.*y(17:o.shz.dgf:end))./o.shz.etaK;
            e33k=(s33p-o.shz.Gk.*y(18:o.shz.dgf:end))./o.shz.etaK;
            
            % total strain = Kelvin + Maxwell strain rates d epsilonM / dt = sigma' / etaM
            e11=e11k+s11p./o.shz.etaM;
            e12=e12k+s12p./o.shz.etaM;
            e13=e13k+s13p./o.shz.etaM;
            e22=e22k+s22p./o.shz.etaM;
            e23=e23k+s23p./o.shz.etaM;
            e33=e33k+s33p./o.shz.etaM;
            
            % rate of strain
            yp(1:o.shz.dgf:end)=e11;
            yp(2:o.shz.dgf:end)=e12;
            yp(3:o.shz.dgf:end)=e13;
            yp(4:o.shz.dgf:end)=e22;
            yp(5:o.shz.dgf:end)=e23;
            yp(6:o.shz.dgf:end)=e33;
            
            % rate of stress
            yp( 7:o.shz.dgf:end)=o.LL{1,1}*e11+o.LL{2,1}*e12+o.LL{3,1}*e13+o.LL{4,1}*e22+o.LL{5,1}*e23+o.LL{6,1}*e33;
            yp( 8:o.shz.dgf:end)=o.LL{1,2}*e11+o.LL{2,2}*e12+o.LL{3,2}*e13+o.LL{4,2}*e22+o.LL{5,2}*e23+o.LL{6,2}*e33;
            yp( 9:o.shz.dgf:end)=o.LL{1,3}*e11+o.LL{2,3}*e12+o.LL{3,3}*e13+o.LL{4,3}*e22+o.LL{5,3}*e23+o.LL{6,3}*e33;
            yp(10:o.shz.dgf:end)=o.LL{1,4}*e11+o.LL{2,4}*e12+o.LL{3,4}*e13+o.LL{4,4}*e22+o.LL{5,4}*e23+o.LL{6,4}*e33;
            yp(11:o.shz.dgf:end)=o.LL{1,5}*e11+o.LL{2,5}*e12+o.LL{3,5}*e13+o.LL{4,5}*e22+o.LL{5,5}*e23+o.LL{6,5}*e33;
            yp(12:o.shz.dgf:end)=o.LL{1,6}*e11+o.LL{2,6}*e12+o.LL{3,6}*e13+o.LL{4,6}*e22+o.LL{5,6}*e23+o.LL{6,6}*e33;
            
            % rate of Kelvin strain
            yp(13:o.shz.dgf:end)=e11k;
            yp(14:o.shz.dgf:end)=e12k;
            yp(15:o.shz.dgf:end)=e13k;
            yp(16:o.shz.dgf:end)=e22k;
            yp(17:o.shz.dgf:end)=e23k;
            yp(18:o.shz.dgf:end)=e33k;
            
        end % end ode
        
    end % end methods
    
end
