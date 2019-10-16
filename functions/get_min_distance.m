function mindist = get_min_distance(xi,yi,X,Y)
    mindist=inf;
    for i=1:length(X)
        disti = sqrt((xi-X(i)).^2 + (yi-Y(i)).^2);
        if disti<mindist
            mindist=disti;
        end
    end
end

