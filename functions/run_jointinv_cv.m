function resid = run_jointinv_cv(logweights,scenario_test,scenario_train)
    % run jointinv inversion on training dataset and report the chi2 value
    % on test dataset.
    % weights are passed in log form to enforce nonnegativity,
    % for jointinv we use 10^weights(i).
    %
    % Eric Lindsey, May 2019
    
    % assign weight values
    for i=1:length(logweights)
        scenario_train.userParams.smoothingWeights{1}{i}=10^logweights(i);
    end
    
    % run the inversion on the training dataset
    scenario_train.run_inversion()

    % compute the weighted residuals on the test dataset
    scenario_test.modelVector = scenario_train.modelVector;
    scenario_test.predVector = scenario_test.designMatrix*scenario_test.modelVector;
    %resid = scenario_test.datasets{1}.covarianceMatrix^0.5*(scenario_test.predVector - scenario_test.dataVector);
    resid = abic_alphabeta(scenario_test.dataVector, scenario_test.modelVector, scenario_test.datasets{1}.covarianceMatrix, scenario_test.designMatrix, scenario_test.smoothingMatrix(1:end-1,:), scenario_train.userParams.smoothingWeights{1}{1}, scenario_train.userParams.smoothingWeights{1}{2});

end
