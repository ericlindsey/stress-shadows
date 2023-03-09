function [mo,mag] = calc_moment_deficit_triangles(patchAreas,slipMag,shearModulus)
    % moment deficit for a triangular patch fault
    % Eric Lindsey, 2022
    
    % unicycle geometry object has an 'area' vector, units in square meters
    
    % we assume units of slipMag are also meters, or meters/year
    
    % shear modulus may be a vector or constant. Units are Pa (not GPa).
    % If not constant, careful about which side of the fault do you mean?
    
    % units are dyne-cm (extra factor of 1e7)
    mo = 1e7 * slipMag .* patchAreas .* shearModulus;
    mag = (2/3) * log10(mo) - 10.7;
    
end