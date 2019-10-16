%% run many realizations of data noise to find the range of acceptable models

fix_alpha=0.029; %abic-determined best alpha for this dataset

colors={'r','b'};

% search over a set of beta values
nbeta=20;
betas=[1e-8,logspace(-3,-2.1,nbeta-1)];

% try several realizations of data noise
Nnoise=100000;
% timing info: 
% 0.47 sec for 1
%  4.4 sec for 10
%   44 sec for 100
%43257 sec for 100000  (~12 hrs)

abics=nan*zeros(2,Nnoise,nbeta);

no_noise_file = load('test_data_2d/zero_noise_data.dat');
zero_noise_data = no_noise_file(:,4);


maxfname='test_data_2d/test_abic_maxmodel.dat';
minfname='test_data_2d/test_abic_minmodel.dat';

%initialize min/max models
minmodel=1*ones(104,2);
maxmodel=0*ones(104,2);
save(maxfname,'maxmodel','-ASCII');
save(minfname,'minmodel','-ASCII');

maxmodel=load(maxfname);
minmodel=load(minfname);

%maxmodel=zeros(length(truemodel),2);
%minmodel=1 + maxmodel;


tic
for expNumber=0:1
    
    for inoise=1:Nnoise
        
        min_abic=1e9;
        
        scenario = Jointinv(expNumber);
        scenario.run_setup();

        % add data noise
        sigma=1;
        scenario.dataVector = zero_noise_data + sigma*randn(size(scenario.dataVector));

        %fix the value of alpha
        scenario.userParams.smoothingWeights{1}{1} = fix_alpha;

        for ibeta=1:nbeta
            %set beta
            scenario.userParams.smoothingWeights{1}{2} = betas(ibeta);
            %run inversion for this parameter set
            scenario.run_inversion();
            %save results
            %chi2red = (scenario.predVector - scenario.dataVector)'*scenario.datasets{1}.covarianceMatrix*(scenario.predVector - scenario.dataVector)/length(scenario.dataVector);
            abic = abic_alphabeta(scenario.dataVector, scenario.modelVector, scenario.datasets{1}.covarianceMatrix, scenario.designMatrix, scenario.smoothingMatrix(1:end-1,:), scenario.userParams.smoothingWeights{1}{1}, scenario.userParams.smoothingWeights{1}{2});

            min_abic = min(abic, min_abic);
            
            if abic<1.1*min_abic
 
                % save results
                maxmodel(:,expNumber+1)=max(maxmodel(:,expNumber+1),1-scenario.modelVector/(100/cosd(10)));
                minmodel(:,expNumber+1)=min(minmodel(:,expNumber+1),1-scenario.modelVector/(100/cosd(10)));
                
                %models(expNumber+1,inoise,ibeta,:) = 1-scenario.modelVector/(100/cosd(10));
                abics(expNumber+1,inoise,ibeta)=abic;             
            end
            
        end
        
    end
    
end
toc

save(maxfname,'maxmodel','-ASCII');
save(minfname,'minmodel','-ASCII');

%% make plots

figure(9),clf

d = scenario.dataVector;
xdata=scenario.datasets{1}.coordinates(:,1)/1e3;
xfault=scenario.sources{1}.geom.xc(:,1)/1e3;
truemodelfname='test_data_2d/test_true_modelVector.dat';
truemodel = load(truemodelfname);
truemodelx=[0;truemodel(:,1)/1e3];
truemodely=[truemodel(1,2)/100*cosd(10); truemodel(:,2)/100*cosd(10)];


%subplot(2,1,2)
plot(truemodel(:,1)/1e3, 1-truemodel(:,2)/(100/cosd(10)),'-k'), hold on

%     subplot(2,1,2)
%     plot(xfault, 1-scenario.modelVector/(100/cosd(10)),'-','linewidth',1,'color',colors{expNumber+1}), hold on
%     subplot(2,1,1)
%     plot(log10(betas),squeeze(chi2reds(expNumber+1,inoise,:)),'-x','color',colors{expNumber+1}),hold on

poly_x = [0; xfault; flipud(xfault); 0; 0];
poly_y_nosc=[minmodel(1,1); minmodel(:,1); flipud(maxmodel(:,1)); maxmodel(1,1); minmodel(1,1)];
poly_y_sc  =[minmodel(1,2); minmodel(:,2); flipud(maxmodel(:,2)); maxmodel(1,2); minmodel(1,2)];

patch(poly_x,1-poly_y_nosc,'r')
patch(poly_x,1-poly_y_sc,'b')
alpha(0.3)

out_poly_nosc = [poly_x, poly_y_nosc];
out_poly_sc = [poly_x, poly_y_sc];
out_truemodel = [truemodelx, truemodely];

nosc_fname='test_data_2d/polygon_noconstraint.dat';
sc_fname='test_data_2d/polygon_stress_constraint.dat';
truemodel_fname = 'test_data_2d/true_model.dat';
save(nosc_fname,'out_poly_nosc','-ASCII');
save(sc_fname,'out_poly_sc','-ASCII');
save(truemodel_fname,'out_truemodel','-ASCII');

    
    
    
      %
    % save the results to a file
    %betafname=['test_data_2d/test_beta_fixalpha_exp' num2str(expNumber) '.dat'];
    %betatest=[betas', chi2red];
    %save(betafname,'betatest','-ASCII')

%     subplot(2,1,2)
%     legend('Unconstrained inversion','Stress-constrained inversion','location','NorthWest')
%     xlabel('log_{10}(beta)')
%     ylabel('reduced \chi^2 with fixed alpha')
%     set(gca,'fontsize',14)
%     ylim([0.7,1.1])
%     
    