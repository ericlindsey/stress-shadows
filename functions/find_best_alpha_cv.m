function [best_alpha,sig_alpha] = find_best_alpha_cv(scenario,Ncv,Nk)
    % use k-fold CV to optimize value of alpha.
    % (beta should be set in the scenario object separately)

    % set starting point, bounds, and options for the optimization
    % note, these are log_10(alpha)
    x0=-2;
    xmin=-3;
    xmax=0;
    lsqnonlin_options = optimoptions(@lsqnonlin,'FunctionTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 500);
    % pass these lsqnonlin parameters and options as a struct
    inv_options.x0=x0;
    inv_options.xmin=xmin;
    inv_options.xmax=xmax;
    inv_options.lsqnonlin_options=lsqnonlin_options;

    % total number of models is Ncv*Nk. Allocate result matrices:
    fullcv_results = zeros(Ncv*Nk,1+length(inv_options.x0));
    
    for i=1:Ncv
        % k-fold cross validation. Will run Nk optimizations for (alpha, [beta]).
        [kcv_results,~]=kfold_cv_jointinv(scenario,Nk,inv_options);

        % save results into larger matrix
        fullcv_results((i-1)*Nk+1:i*Nk,:)=kcv_results;

    end

    best_alpha=mean(fullcv_results(:,1));
    sig_alpha=std(fullcv_results(:,1));
    
end
