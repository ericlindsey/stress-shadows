classdef gimbutas12 < unicycle.greens.earthModel
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
        function o=gimbutas12(G,nu)
            % GIMBUTAS12 is a class providing the necessary functions to
            % compute the stress interactions between sources and receivers
            % and among the receivers.
            %
            %   earthModel = greens.gimbutas12(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space.
            %
            % SEE ALSO: unicycle, hmmvp
            
            if (0==nargin)
                return
            end
            
            o.G=G;
            o.nu=nu;
        end
        
        function [varargout]=stressKernels(obj,src,rcv)
            % STRESSKERNELS computes the stress on receiver faults due to
            % motion of dislocations.
            %
            % rcv - receiver fault
            %
            % SEE ALSO: unicycle, hmmvp
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeStressKernelsGimbutas12(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=displacementKernels(obj,src,x,vecsize)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of dislocations.
            %
            % src - source fault
            %
            % SEE ALSO: unicycle, unicycle.greens.computeDisplacementKernelsMeade07
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeDisplacementKernelsMeade07(src,obj.nu,x,vecsize);
        end
    end % methods
end % class definition
