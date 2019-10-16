function Rx = rot_matrix(lat,lon)
% Given a list of lat,lon, create a rotation matrix.
% The rotation matrix exists as a 2x3 block for each site:
% [[ -sin(lat)cos(lon) -sin(lat)sin(lon) cos(lat)]
%      [  sin(lon)         -cos(lon)          0]]
%     Thus, the final result has dimensions 2nx3
Rx=[];
for i =1:length(lat)
    Ri=[[-sind(lat(i))*cosd(lon(i)), -sind(lat(i))*sind(lon(i)), cosd(lat(i))];...
        [sind(lon(i)), -cosd(lon(i)), 0]];
    Rx=[Rx;Ri];
end

end
        
   
        
