function [cv_results,cv_models]=kfold_cv_jointinv(scenario, Nk, inv_options)
    % k-fold cross validation on a jointinv object, repeated Ncv times
    % Eric Lindsey, May 2019
    
    % create matrices to store results
    cv_results = zeros(Nk,1+length(inv_options.x0));
    cv_models = zeros(Nk,length(scenario.modelVector));
    
    % create k disjoint lists ('k-folds'). Stored as a cell list, due to
    % possibly different lengths
    klists = create_kfold_ind(length(scenario.dataVector),Nk);
    
    % loop: use each disjoint list as the test set once
    for i=1:Nk
        % split dataset into training and test datasets
        Itest=klists{i};
        [scenario_train,scenario_test]=split_jointinv_traintest(scenario,Itest);

        % run the inversion to find best hyperparameters to fit the train/test datasets
        [optvals,resnorm]=lsqnonlin( @(x) run_jointinv_cv(x, scenario_train, scenario_test) , ...
            inv_options.x0, inv_options.xmin, inv_options.xmax, inv_options.lsqnonlin_options );        

        % store the results from each k-fold test
        cv_results(i,:)=[optvals,resnorm];
        cv_models(i,:)=scenario_train.modelVector;

    end
    
end