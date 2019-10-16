function outgps = deuplicate_GSRM_gps_data(ingps)

    % find GPS data from the GSRM dataset that is repeated several times,
    % and keep only the version with smaller reported error.
    %
    % E. Lindsey, July 2019

    gpspts = [ingps(:,2), ingps(:,1)];
    gpsuncert = mean([ingps(:,5), ingps(:,6)],2);
    E=wgs84Ellipsoid('km');
    
    Ikeep=true(length(gpspts),1);
    
    for i=1:length(gpspts)
        if Ikeep(i)==1
            ptdist = distance(gpspts(i,:),gpspts,E);
            Idup = find(ptdist<1); %threshold for duplicate: 1km.
            ndup=length(Idup);
            if ndup > 1 
                Ibest = find(gpsuncert(Idup) == min(gpsuncert(Idup)),1);
                Idup(Ibest)=[];
                Ikeep(Idup)=0;
            end
        end
    end
    outgps=ingps(Ikeep,:);
    
end