classdef flowlaw < handle
    properties(Constant)
        % Mantle density (kg/m^3)
        Rm=3330;
        % Specific heat (J/kg/K)
        Cp=1171;
        % Thermal conductivity (W/m/K)
        k=3.138;
        % Acceleration due to gravity (m/s^2)
        g=9.8;
        % Universal gas constant (J/mol/K)
        R=8.3144 ;
    end 
    
    properties
        % Pre-exponential factor (MPa^-n )
        A;
        % Activation energy (J/mol)
        Q;
        % Stress exponent
        n;
        % Water fagucity/content exponent
        r;
        % Activation volume (m^3/mol)
        V;
        % Grain-size (um)
        d;
        % Grain-size exponent
        m;  
        % Water content (H/10^6 Si)
        coh;
        % Temperature (K)
        T;             
   end

end
