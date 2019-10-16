
% choose the input file and create blank object
expNumber=12; 
scenario = Jointinv(expNumber);

% read the input file and set up matrices
scenario.run_setup();

% run the inversion
scenario.run_inversion()

% compute desired products
results = calc_coupling_result_components(scenario);

% display model misfit
disp(['chi2: ' num2str(results.chi2)])

% generate custom plots
plot_coupling_inversion(scenario,results,expNumber);

% output data for plotting in GMT
save_coupling_inversion(scenario,results);

