function save_geom_for_GMT(output_filename,geom,values,lat0,lon0)
    % function to convert a geometry object to a GMT-plottable file
    
    % output file: contains 4 lines for each triangle:
    % >val1
    % lon1 lat1
    % lon2 lat2
    % lon3 lat3
    % >val2 
    %... etc
    assert(length(values)==geom.N,'Error: values vector must have the same length as the number of patches in Unicycle geometry object.')
    
    % convert vertices (unicycle format is meters xyz) to lat/lon
    [vlat,vlon] = xy_to_latlon_polyconic(geom.x(:,1)/1e3, geom.x(:,2)/1e3, lon0, lat0);
    vdep = abs(geom.x(:,3)/1e3);
    fileID = fopen(output_filename,'w');
    % for i=1:geom.N
    %     fprintf(fileID,'> -Z%.9f\n',values(i));
    %     fprintf(fileID,'%.9f %.9f\n',vlon(geom.vertices(i,1)),vlat(geom.vertices(i,1)));
    %     fprintf(fileID,'%.9f %.9f\n',vlon(geom.vertices(i,2)),vlat(geom.vertices(i,2)));
    %     fprintf(fileID,'%.9f %.9f\n',vlon(geom.vertices(i,3)),vlat(geom.vertices(i,3)));
    % end
    for i=1:geom.N
        fprintf(fileID,'> -Z%.9f\n',values(i));
        fprintf(fileID,'%.9f %.9f\n',vlat(geom.vertices(i,1)),vdep(geom.vertices(i,1)));
        fprintf(fileID,'%.9f %.9f\n',vlat(geom.vertices(i,2)),vdep(geom.vertices(i,2)));
        fprintf(fileID,'%.9f %.9f\n',vlat(geom.vertices(i,3)),vdep(geom.vertices(i,3)));
    end
    fclose(fileID);

    disp(['Wrote ' output_filename '. Use gmt psxy -L -C<colormap> to plot polygons in lat,depth coordinates.'])
    
end