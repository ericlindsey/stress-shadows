
% set up a simple model to test the difference between L = [ 1 1 ...] and L=eye(N).
N=length(scenario.modelVector);

G=scenario.designMatrix;
d=scenario.dataVector;

Lsmooth = scenario.smoothingMatrix(1:N,1:N);

% case 1: 
Lsum = ones(1,N);
% case 2: 
Lval = eye(N);

betas = logspace(-5,1,30);
nbeta=length(betas);

figure(1),clf
title('Sum constraint'),hold on
figure(2),clf
title('Value constraint'),hold on

cmap = parula(nbeta);

lb=scenario.lowerBounds;
ub=scenario.upperBounds;

%plot(d,'k.'),hold on

alpha = 10^-1.65;

for ibeta=1:nbeta
    msum = lsqlin([G; alpha*Lsmooth; betas(ibeta)/(sqrt(N/2))*Lsum],[d; zeros(N,1); 0],[],[],[],[],lb,ub);
    mval = lsqlin([G; alpha*Lsmooth; betas(ibeta)*Lval],[d; zeros(N,1); zeros(N,1)],[],[],[],[],lb,ub);
    figure(1)
    plot(msum,'color',cmap(ibeta,:)),hold on
    figure(2)
    plot(mval,'color',cmap(ibeta,:)), hold on
    
    misfitsum(ibeta)=(G*msum - d)'*(G*msum - d);
    misfitval(ibeta)=(G*mval - d)'*(G*mval - d);

    
end

figure(3),clf,hold on
plot(misfitsum,'b')
plot(misfitval,'r')
legend('sum','val')
%xlabel('log_{10}(\beta)')





%%


% set up a simple model to test the difference between L = [ 1 1 ...] and L=eye(N).
N = 1;

G=[ones(N,1), (1:N)']

d=2*(1:N)' + ones(N,1)+5*randn(N,1);

% case 1: 
Lsum = [1 1];
% case 2: 
Lval = eye(2);

betas = logspace(-3,5,100);

figure(1),clf

lb=zeros(2,1);
ub=inf+zeros(2,1);

%plot(d,'k.'),hold on

for ibeta=1:nbeta
    msum = lsqlin([G;betas(ibeta)./sqrt(2)*Lsum],[d; 0],[],[],[],[],lb,ub);
    mval = lsqlin([G;betas(ibeta)*Lval],[d; zeros(2,1)],[],[],[],[],lb,ub);
    plot(msum,'bx'),hold on
    plot(mval,'Rs')
end
%xlabel('log_{10}(\beta)')
title(['N = ' num2str(N)])
xlim([-1,3])