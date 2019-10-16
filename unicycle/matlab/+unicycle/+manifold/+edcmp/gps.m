classdef gps < unicycle.manifold.gps
    properties
        
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = gps(filename,evl,vecsize,prefix)
            % GPS is a class representing GPS data and Green's functions.
            %
            %   gps = unicycle.manifold.edcmp.gps(network,evl,vecsize,prefix)
            %
            % creates a instance of GPS data.
            %
            % INPUT:
            %
            % network  filename, for example 'sopac.dat', containing
            %
            %          # i NAME x1 x2 x3
            %            1 GPS1  0  0  0
            %            2 GPS2  1  1  0
            %
            % evl      object of type ode.evolution containing a time
            %          series of fault displacement
            % vecsize  number of components of displacement vector
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            
            import unicycle.greens.edcmp
            
            obj.vecsize=vecsize;
            
            [~,obj.stationName,x1,x2,x3]=...
                textread(filename,'%d %s %f %f %f','commentstyle','shell');
            obj.x=[x2,x1,0*x3];
            obj.D=size(x2,1);
            
            % source Green's functions (includes strike slip and dip slip)
            obj.G=edcmp.G(evl.src,prefix,obj.x,vecsize);
            % receiver Green's functions (includes strike slip and dip slip)
            obj.H=edcmp.G(evl.rcv,prefix,obj.x,vecsize);
            
            % builds forward models of geodetic data if simulation exists
            if ~isempty(evl.y)
                obj.simulation(evl);
            end
            
        end % constructor
        
    end % methods
    
end
