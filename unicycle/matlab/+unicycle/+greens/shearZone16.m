classdef shearZone16 < unicycle.greens.earthModel
    properties
        % rigidity
        G;
        % Poisson's ratio
        nu;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=shearZone16(G,nu)
            % SHEARZONE16 is a class providing the necessary functions to
            % compute the stress interactions between source and receiver
            % shear zones.
            %
            %   earthModel = greens.shearZone16(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space.
            %
            % SEE ALSO: unicycle
            
            if (0==nargin)
                return
            end
            
            assert(0<=G,'unicycle.greens.shearZone16::rigidity must be positive.')
            assert(nu<=0.5,'unicycle.greens.shearZone16::Poisson''s ratio should be lower than 0.5.')
            assert(-1<=nu,'unicycle.greens.shearZone16::Poisson''s ratio should be greater than -1.')
            
            o.G=G;
            o.nu=nu;
        end
        
        function [varargout]=tractionKernels(obj,src,rcv)
            % TRACTIONKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % rcv - receiver fault
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeTractionKernelsVerticalShearZone(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=stressKernels(obj,src,rcv)
            % STRESSKERNELS computes the stress on receiver faults due to
            % motion of rectangular dislocations in a half space.
            %
            % rcv - receiver shear zone
            %
            % SEE ALSO: unicycle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeStressKernelsVerticalShearZone(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=displacementKernels(obj,src,x,vecsize)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of rectangular dislocations in a half space.
            %
            % src - source fault
            %
            % SEE ALSO: unicycle
            import unicycle.greens.*
            varargout = cell(1,nargout);
            [varargout{:}]=computeDisplacementKernelsVerticalShearZone(src,obj.nu,x,vecsize);
        end
    end % methods
end % class definition
