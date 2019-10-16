function save_geom_for_interp_GMT(output_filename,geom,values,lat0,lon0)
    % function to convert a geometry object to a GMT-plottable file
    
    % output file: contains 4 lines for each triangle:
    % >val1
    % lon1 lat1
    % lon2 lat2
    % lon3 lat3
    % >val2 
    %... etc
    assert(length(values)==geom.N,'Error: values vector must have the same length as the number of patches in Unicycle geometry object.')
    
    % convert patch centers (unicycle format is meters xyz) to lat/lon
    [plat,plon] = xy_to_latlon_polyconic(geom.xc(:,1)/1e3, geom.xc(:,2)/1e3, lon0, lat0);
    
    outdat=[plon plat values]; %#ok<NASGU>
    save(output_filename,'outdat','-ASCII');

    disp(['Wrote ' output_filename '. Use gmt surface to interpolate values in geographic coordinates.'])
    
end