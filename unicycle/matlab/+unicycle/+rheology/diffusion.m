classdef diffusion < flowlaw
   methods
        function obj=diffusion(varargin)
            % -----------------------------------------------------------------
            % Diffusion : This is class to set the diffusion creep parameters
            % depending on the mineral. 
            % Default values of parameters.
            % d (Grain size) = 10 microns
            % T (temperature) = 1200 & 1350 C (quartz & olivine) 
            % 
            % These values are from following papers.
            % Quartz :
            % Rutter E.H. and K. H. Brodie., "Experimental grain size-sensitive 
            % flow of hot-pressed Brazilian quartz aggregates",Journal of structural
            % geology, 2004
            % 
            % Quartz flow law for dislocation creep is 
            %   .                    r  -m    |       |
            %   e = A (sigma)  (f   )  d  exp | - Q   |
            %                    H O          | ----- |
            %                     2           |   RT  |
            % Olivine :  
            % Hirth.G., D.Kohlstedt, "Rheology of the Upper Mantle and the Mantle 
            % Wedge: A View from the Experimentalists", 
            % Inside subduction zone Factory, 2003
            %
            % Olivine flow law for diffusion creep is 
            %   .                   r  -m     |          |
            %   e = A (sigma) (C   )  d   exp | - Q + pV |
            %                   OH            | -------- |
            %                                 |     RT   |
            % -----------------------------------------------------------------

            if (0==nargin)
                return;
            end

            wet=1;
            mineral=varargin{1}; 
            if (2==nargin)
                wet=varargin{2};
            end
            switch (mineral)
                % Quartz 
                case 1
                    if (1==wet)
                        obj.A=0.4;
                        obj.Q=220e3;
                        obj.n=1;
                        obj.r=1;
                        obj.V=0;
                        obj.d=10;
                        obj.m=2;
                        obj.coh=200;
                        obj.T=1200+273;
                    else
                        warning('Dry quartz flow law values are not well known');
                    end
                % Olivine
                case 2
                    % Wet diffusion constant C_OH case
                    if (1==wet)
                        obj.A=1e9;
                        obj.Q=335e3;
                        obj.n=1;
                        obj.r=1;
                        obj.V=4e-6;
                        obj.d=10;
                        obj.m=3;
                        obj.coh=1000;
                        obj.T=1350+273;
                    else
                    % Dry diffusion 
                        obj.A=1.5e9;
                        obj.Q=375e3;
                        obj.n=1;
                        obj.r=0;
                        obj.V=4e-6;
                        obj.d=10;
                        obj.m=3;
                        obj.coh=50;
                        obj.T=1350+273;
                    end 
                otherwise 
                    warning('Unexpected mineral. No values of have been set')
            end 
        end
    end
end
