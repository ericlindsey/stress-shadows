function abic = abic_smoothingonly(scenario)
% function to calculate A-BIC value given 
% d - data 
% m - model
% G - design matrix for a linear problem
% W - weighting function (generally inverse of data covariance matrix)
% alpha - smoothing hyper-parameter
% L - smoothing matrix (generally finite-difference Laplacian operator) 
% OUTPUT
% abic - information criteria
% Rishav Mallick, EOS, 2019
% Modified by E. Lindsey, June 2019
% reference Fukuda and Johnson (2008; 2010), Funning et al., (2014)
%%
d = scenario.dataVector;
m = scenario.modelVector;
%W = inv(scenario.datasets{1}.covarianceMatrix);

G = scenario.designMatrix;
smoothingLength = scenario.smoothingLengths{1}{1};
L = scenario.smoothingMatrix(1:smoothingLength,:);
alpha = scenario.userParams.smoothingWeights{1}{1};

Wdgm = scenario.datasets{1}.covarianceMatrix \ (d - G*m);
WG = scenario.datasets{1}.covarianceMatrix \ G;

% compute bet-fitting model:
% modified by Eric to make the model vector "m" a user input instead of computed.
%m = ( G'*W*G + (alpha^2)*(L'*L) )\(G'*W*d);

%abic1 = length(d)*log( (d-G*m)'*W*(d-G*m) + (alpha^2)*(L*m)'*(L*m) );
abic1 = length(d)*log( (d-G*m)'* Wdgm + (alpha^2)*(L*m)'*(L*m) );

eigprior = (alpha^2)*abs(eig(L'*L));
eigprior(eigprior<=0) = [];
abic2 = sum(log(eigprior));

%eigdes = abs(eig(G'*W*G + LL));
eigdes = eig(G'*WG + alpha^2*(L'*L));
eigdes(eigdes<=0) = [];
abic3 = sum(log(eigdes));

% return ABIC value
abic = abic1 - abic2 + abic3;

figure(4),hold on
%plot(log10(alpha),abic1,'b.')
plot(log10(alpha),(d-G*m)'* Wdgm,'rx')
plot(log10(alpha),(alpha^2)*(L*m)'*(L*m),'bx')
title('abic1')
% 
% figure(5),hold on
% plot(log10(alpha),abic2-500*log10(alpha),'rx')
% title('abic2 - 500*log10(alpha)')
% figure(6),hold on
% plot(log10(alpha),abic3-500*log10(alpha),'gx')
% title('abic3 - 500*log10(alpha)')
 figure(7),hold on
 plot(log10(alpha),abic,'k.')
 title('abic')
% figure(8),hold on
% plot(log10(alpha),abic3-abic2,'mx')
% title('abic3-abic2')


end