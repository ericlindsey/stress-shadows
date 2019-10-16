classdef nikkhoo15 < unicycle.greens.earthModel
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
        function o=nikkhoo15(G,nu)
            % NIKKHOO15 is a class providing the necessary functions to
            % compute the stress interactions between sources and receivers
            % and among the receivers using triangular dislocations.
            %
            %   earthModel = greens.nikkhoo15(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space.
            %
            % SEE ALSO: unicycle, geometry.triangle, geometry.triangleReceiver
            
            if (0==nargin)
                return
            end
            
            assert(0<=G,'unicycle.greens.nikkhoo15::rigidity must be positive.')
            assert(nu<=0.5,'unicycle.greens.nikkhoo15::Poisson''s ratio should be lower than 0.5.')
            assert(-1<=nu,'unicycle.greens.nikkhoo15::Poisson''s ratio should be greater than -1.')
            
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
            [varargout{:}]=computeTractionKernelsNikkhoo15(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=stressKernels(obj,src,rcv)
            % STRESSKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % rcv - point
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeStressKernelsNikkhoo15(src,rcv,obj.G,obj.nu);
        end
        
        function [varargout]=displacementKernels(obj,src,x,vecsize)
            % DISPLACEMENTKERNELS computes the stress on receiver faults due to
            % motion of triangular dislocations.
            %
            % src - source fault
            %
            % SEE ALSO: unicycle, geometry.triangle
            
            import unicycle.greens.*
            
            varargout = cell(1,nargout);
            [varargout{:}]=computeDisplacementKernelsNikkhoo15(src,obj.nu,x,vecsize);
        end
    end % methods
end % class definition
