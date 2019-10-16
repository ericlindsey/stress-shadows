function optvals = find_best_alphabeta_abic(scenario)
    % use ABIC to optimize value of alpha and beta for jointinv

    % set starting point for the optimization
    % note, these are log_10(hyperparameter)
    x0=[-2,-3];

    [optvals,~] = fminunc( @(x) run_jointinv_abic(x, scenario), x0);  
    
end
