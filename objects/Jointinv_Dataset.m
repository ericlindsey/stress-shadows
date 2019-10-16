classdef Jointinv_Dataset < handle
    % DATASET contains the general properties required by jointinv for any type of data
    
    % It is inherited by:
    
    %   Static_GPS_Dataset
    
    properties
        
        fileName         % name of the file used to load the dataset
        coordinates      % n x [2, 3, or 4] - independent variables (locations, times) for each data point
        
        dataVector       % N x 1 - column vector containing the dependent variables (data)
                         % the data should be unraveled by point, not by component.
                         % eg. (pt1_E, pt1_N, pt1_U, pt2_E, pt2_N, pt2_U, ...)
                         % note that N = n*k
        
        predVector       % predicted value of data from the inversion result
        
        covarianceMatrix % N x N - covariance matrix containing the data uncertainties/weights
        
        numComponents    % dimension k of the data (data is unraveled as a 
                         % vector; this value is needed to determine how
                         % the coordinates match to the data)
                         
    end
    
    methods
        % no default methods for this class
    end
    
end
