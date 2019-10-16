function check_triangleReceiver_signs(triangleReceiver)
    % due to an undiagnosed bug inside Unicycle, the stresses/displacements
    % generated by a triangleReceiver object will be inconsistent depending
    % on the sense of orientation of the triangle vertices. To enforce
    % consistency with the rest of jointinv, we require that the normal 
    % vector formed by the cross product of the first and second sides of 
    % the triangle points downward.
    %
    % Eric Lindsey, June 2019
    
    p1 = triangleReceiver.x(triangleReceiver.vertices(:,1),:);
    p2 = triangleReceiver.x(triangleReceiver.vertices(:,2),:);
    p3 = triangleReceiver.x(triangleReceiver.vertices(:,3),:);

    v1 = p2 - p1;
    v2 = p3 - p2;

    cp=cross(v1,v2);

    assert(max(cp(:,3))<=0,'Error: For use with Jointinv, vertices for all triangles in triangleReceiver object must be ordered in a clockwise sense.')

end
