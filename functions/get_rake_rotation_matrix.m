function [Rmat,Vpl] = get_rake_rotation_matrix(rakedata)
    % compute a rotation matrix and slip rate for each patch from a dataset
    % containing: [rakeangle,strike_rate,dip_rate] in the first 3 columns
    %
    % Eric Lindsey, August 2018
    % Modified E. Lindsey, June 2019
    
    Vpl = (rakedata(:,2).^2+rakedata(:,3).^2).^0.5;
    rakeangle = rakedata(:,1);
    Rmat = [diag(cosd(rakeangle)), diag(-sind(rakeangle)); diag(sind(rakeangle)), diag(cosd(rakeangle))];

end
