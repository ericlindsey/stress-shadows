function abic = abic_alphabeta(scenario)
% function to calculate A-BIC value given 
% d - data 
% m - model
% G - design matrix for a linear problem
% W - weighting function (generally inverse of data covariance matrix)
% alpha - smoothing hyper-parameter
% L - smoothing matrix (generally finite-difference Laplacian operator) 
% beta - length hyper-parameter
% OUTPUT
% abic - information criteria
% Rishav Mallick, EOS, 2019
% Modified by E. Lindsey, June 2019
% reference Fukuda and Johnson (2008; 2010), Funning et al., (2014)

d = scenario.dataVector;
m = scenario.modelVector;
W = inv(scenario.datasets{1}.covarianceMatrix);
G = scenario.designMatrix;
smoothingLength = scenario.smoothingLengths{1}{1};
L = scenario.smoothingMatrix(1:smoothingLength,:);
alpha = scenario.userParams.smoothingWeights{1}{1};
beta = scenario.userParams.smoothingWeights{1}{2};

% note, this method assumes you are applying beta to all model parameters equally.
npatch = length(scenario.modelVector);

% compute bet-fitting model:
% modified by Eric to make the model vector "m" a user input instead of computed.
%m = (G'*W*G + (alpha^2)*(L'*L) + (beta^2)*eye(npatch))\(G'*W*d);

abic1 = length(d)*log((d-G*m)'*W*(d-G*m) + (alpha^2)*(L*m)'*(L*m) + (beta^2)*(m'*m));

LL = (alpha^2)*(L'*L) + (beta^2)*eye(npatch);
eigprior = abs(eig(LL));
eigprior(eigprior<=0) = [];
abic2 = sum(log(eigprior));

eigdes = abs(eig(G'*W*G + LL));
eigdes(eigdes<=0) = [];
abic3 = sum(log(eigdes));

% return ABIC value
abic = abic1 - abic2 + abic3;

end