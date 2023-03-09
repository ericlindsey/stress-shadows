
% choose the input file and create blank object
expNumber=101; 
scenario = Jointinv(expNumber);

% read the input file and set up matrices
scenario.run_setup();

% run the inversion
scenario.run_inversion()

% generate custom plots
figNum = expNumber;
figure(figNum),clf,hold on
plot_jointinv_slipmodel(scenario,figNum,100);

% output data for plotting in GMT
% save_eq_slipmodel(scenario);
