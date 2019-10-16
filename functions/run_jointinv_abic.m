function resid = run_jointinv_abic(logweights,scenario)
    % run jointinv inversion on training dataset and report the ABIC value
    % 
    % weights are passed in log form to enforce nonnegativity,
    % for jointinv we use 10^weights(i).
    %
    % Eric Lindsey, May 2019
    
    % assign weight values
    for i=1:length(logweights)
        scenario.userParams.smoothingWeights{1}{i}=10^logweights(i);
    end
    
    % run the inversion
    scenario.run_inversion()

    % old: get (W^0.5)*(m-d)
    %resid = scenario_test.datasets{1}.covarianceMatrix^0.5*(scenario_test.predVector - scenario_test.dataVector);
    
    % new: get ABIC value
    resid = abic_alphabeta(scenario);

end
