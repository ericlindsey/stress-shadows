function [E,N,U] = predict_enu_model(source, strikeslip, dipslip, coords, kernelFolder)
    % Output predicted E,N,U slip from a slip model and Unicycle source
    % object
    %
    % Eric Lindsey, 2018
    
    G=unicycle_displacement_kernel(source, coords, [1,2], kernelFolder);
    m=[strikeslip(:); dipslip(:)];
    pred=G*m;
    E=pred(1:3:end);
    N=pred(2:3:end);
    U=pred(3:3:end);

end
