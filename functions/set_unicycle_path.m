function set_unicycle_path(unicyclepath)
    % Eric Lindsey, August 2018
    
    if (~contains(path,unicyclepath))
        disp(['Setting unicycle path: ', unicyclepath])
        addpath(genpath(unicyclepath));
        import unicycle.*
    end
end