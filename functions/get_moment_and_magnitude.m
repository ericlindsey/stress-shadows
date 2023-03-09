function  [Mo,Mw] = get_moment_and_magnitude(geom, slipmag)
    % return scalar moment and magnitude, for a unicycle geom object and
    % vector of slip magnitudes 
    assert(length(slipmag) == geom.N, 'Error: For moment calculation, slip magnitude must have length equal to geom.N')
    
    if isa( geom.earthModel, 'unicycle.greens.nikkhoo15')
        area = geom.area;
    elseif isa( geom.earthModel, 'unicycle.greens.okada92')
        area = geom.L .* geom.W ;
    else
        error('cannot compute magnitude, unknown source type')
    end
    
    % note, factors are 1e7 for conversion from N-m to dyne-cm, and 1e6 for
    % conversion of unicycle shear modulus G from MPa to Pa.
    Mo = 1e7 * 1e6 * geom.earthModel.G * sum(area .* slipmag);
    Mw = (2/3) * log10(Mo) - 10.7;

end
