classdef maxwell < unicycle.ode.evolution
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj=maxwell(varargin)
            % MAXWELL is a class encapsulating the geometry and physical 
            % properties of shear zones with a linear Maxwell viscoelastic 
            % rheology.
            %
            % The distributed flow is controlled by
            %
            %   d epsilon   tau
            %   --------- = ---
            %       dt      eta
            %
            % EXAMPLE:
            %
            %   evl = ode.maxwell(src,shz,evt);
            %
            % creates a instance from source (src), shear zone (shz) and
            % events (evt).
            %
            % type methods(evl) and properties(evl) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle
            
            obj=obj@unicycle.ode.evolution(varargin{:});
            obj.shz.dgf=12;
            
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
            %    dE   sigma
            %    -- = -----
            %    dt    eta
            %
            % where eta is the linear Maxwell viscosity.
            
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
            
            % Maxwell strain rates d epsilon / dt = sigma' / etaM
            e11=s11p./o.shz.etaM;
            e12=s12p./o.shz.etaM;
            e13=s13p./o.shz.etaM;
            e22=s22p./o.shz.etaM;
            e23=s23p./o.shz.etaM;
            e33=s33p./o.shz.etaM;
            
            % rate of state
            yp(1:o.shz.dgf:end)=e11;
            yp(2:o.shz.dgf:end)=e12;
            yp(3:o.shz.dgf:end)=e13;
            yp(4:o.shz.dgf:end)=e22;
            yp(5:o.shz.dgf:end)=e23;
            yp(6:o.shz.dgf:end)=e33;
            
            % shear stress rate of change
            yp( 7:o.shz.dgf:end)=o.LL{1,1}*e11+o.LL{2,1}*e12+o.LL{3,1}*e13+o.LL{4,1}*e22+o.LL{5,1}*e23+o.LL{6,1}*e33;
            yp( 8:o.shz.dgf:end)=o.LL{1,2}*e11+o.LL{2,2}*e12+o.LL{3,2}*e13+o.LL{4,2}*e22+o.LL{5,2}*e23+o.LL{6,2}*e33;
            yp( 9:o.shz.dgf:end)=o.LL{1,3}*e11+o.LL{2,3}*e12+o.LL{3,3}*e13+o.LL{4,3}*e22+o.LL{5,3}*e23+o.LL{6,3}*e33;
            yp(10:o.shz.dgf:end)=o.LL{1,4}*e11+o.LL{2,4}*e12+o.LL{3,4}*e13+o.LL{4,4}*e22+o.LL{5,4}*e23+o.LL{6,4}*e33;
            yp(11:o.shz.dgf:end)=o.LL{1,5}*e11+o.LL{2,5}*e12+o.LL{3,5}*e13+o.LL{4,5}*e22+o.LL{5,5}*e23+o.LL{6,5}*e33;
            yp(12:o.shz.dgf:end)=o.LL{1,6}*e11+o.LL{2,6}*e12+o.LL{3,6}*e13+o.LL{4,6}*e22+o.LL{5,6}*e23+o.LL{6,6}*e33;
            
        end % end ode
        
    end % end methods
    
end
