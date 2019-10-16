classdef dislocation < flowlaw
   methods

        function obj=dislocation(varargin)
            % -----------------------------------------------------------------
            % Dislocation : This is class to set the dislocation creep parameters
            % depending on the mineral. 
            % Default values of parameters.
            % d (Grain size) = 10 microns
            % T (temperature) = 1200 & 1350 C (quartz & olivine) 
            % 
            % These values are from following papers.
            % Quartz :
            % Rutter E.H. and K. H. Brodie., "Experimental intracrystalline plastic 
            % flow in hot-pressed synthetic quartzite prepared from Brazilian 
            % quartz crystals", Journal of structural geology, 2004
            %
            % Quartz flow law for dislocation creep is 
            %   .            n       r  -m    |       |
            %   e = A (sigma)  (f   )  d  exp | - Q   |
            %                    H O          | ----- |
            %                     2           |   RT  |
            % Olivine :  
            % Hirth.G., D.Kohlstedt, "Rheology of the Upper Mantle and the Mantle 
            % Wedge: A View from the Experimentalists", 
            % Inside subduction zone Factory, 2003
            %
            % Olivine flow law for dislocation creep is 
            %   .            n       r      |          |
            %   e = A (sigma)  (C   )   exp | - Q + pV |
            %                    OH         | -------- |
            %                               |     RT   |
            %  sigma (deviatoric stress) is in MPa
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
                        obj.A=1.17e-5;
                        obj.Q=242e3;
                        obj.n=3;
                        obj.r=1;
                        obj.V=0;
                        obj.d=10;
                        obj.m=0;
                        obj.coh=1000;
                        obj.T=1200+273;
                    else
                        warning('Dry quartz flow law values are not well known');
                    end
                % Olivine
                case 2
                    % Wet dislocation constant C_OH case
                    if (1==wet)
                        obj.A=90;
                        obj.Q=480;
                        obj.n=3.5;
                        obj.r=1.2;
                        obj.V=11e-6;
                        obj.d=10;
                        obj.m=0;
                        obj.coh=1000;
                        obj.T=1350+273;
                    else
                        obj.A=1.1e5;
                        obj.Q=530e3;
                        obj.n=3.5;
                        obj.r=0;
                        obj.V=14e-6; % Karato & Jung, 2003 
                        obj.d=10;
                        obj.m=0;
                        obj.coh=50;
                        obj.T=1350+273;
                    end 
                otherwise 
                    warning('Unexpected mineral. No values of have been set')
            end 
        end
    end
end
