function optvals = find_best_alpha_abic(scenario)
    % use ABIC to optimize value of alpha for jointinv

    % set starting point, bounds, and options for the optimization
    % note, these are log_10(hyperparameter)
    x0=-2;
    xmin=-4;
    xmax=0;
    
    lsqnonlin_options = optimoptions(@lsqnonlin,'FunctionTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 500);
    % pass these lsqnonlin parameters and options as a struct
    inv_options.x0=x0;
    inv_options.xmin=xmin;
    inv_options.xmax=xmax;
    inv_options.lsqnonlin_options=lsqnonlin_options;
    
    [optvals,~]=lsqnonlin( @(x) run_jointinv_abic(x, scenario) , ...
        inv_options.x0, inv_options.xmin, inv_options.xmax, inv_options.lsqnonlin_options );  
    
end
