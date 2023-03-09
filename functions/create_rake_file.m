function rake_out = create_rake_file(geom, azimuth, rate, fname)
    % function to create a 'rake file' with given rake angle relative to north. 
    % Computes the angle relative to the patch strike/dip vectors in the
    % unicycle 'geom' object, and scales by the 'rate' parameter. Inputs
    % azimuth and rate can be scalars, or vectors with length == geom.N.
    % Eric Lindsey, July 2019

    assert(length(azimuth) == length(rate) && (length(azimuth) == 1 || length(azimuth) == geom.N), 'Error: azimuth and rate inputs must be scalar or have length geom.N');
    
    rakeE = rate .* sind(azimuth);
    rakeN = rate .* cosd(azimuth);
    
    rake_ss = rakeE .* geom.sv(:,1) + rakeN .* geom.sv(:,2);
    rake_ds = -rakeE .* geom.sv(:,2) + rakeN .* geom.sv(:,1);

    rake_angle=rad2deg(atan2(rake_ds,rake_ss));

    if length(rakeE) == 1
        rakeE = rakeE*ones(size(rake_angle));
        rakeN = rakeN*ones(size(rake_angle));
    end
    
    rake_out=[rake_angle,rake_ss,rake_ds,rakeE,rakeN];
    save(fname,'rake_out','-ASCII')

end
