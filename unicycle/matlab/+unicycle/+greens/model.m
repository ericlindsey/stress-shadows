classdef model < handle
    % MODEL is a class containing functions to compute forward models
    % between fault slip and surface observation points.
    properties (Constant)
        % fault degrees of freedom (strike slip and dip slip)
        dgf=2;
    end
    methods (Static)
        function myG=G()
            % G=G(src,nu,x) computes inversion kernels (Green's 
            % functions) to connect fault slip to surface deformation.
            myG=[];
        end
    end % end Static method
end

