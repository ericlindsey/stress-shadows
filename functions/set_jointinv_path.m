function set_jointinv_path(varargin)
    %
    %SET_JOINTINV_PATH(folder)
    %
    % add specified folder to path (optional). If no argument is supplied,
    % script locates the directory one above the current one. (this script 
    % resides in the functions/ subfolder of Jointinv).
    %
    % Eric Lindsey, Sept 2018
    %
    if nargin == 1
        jointinvpath = varargin{1};
    else
        thispath = pwd;
        idcs = strfind(thispath,filesep);
        jointinvpath = thispath(1:idcs(end)-1); 
    end
    % set the path only if it is not already set.
    if (~contains(path,jointinvpath))
        disp(['Setting jointinv path: ', jointinvpath])
        addpath(genpath(jointinvpath));
    end
end