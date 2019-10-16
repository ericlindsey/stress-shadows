function samples = mh(varargin)
% res = MH(y,x0,@genfunc,[LB,UB,params])
%
% INPUT:
%   y        = data (Nx1)
%   x0       = initial parameters (Px1)
%   @genfunc = generative model
%
%   Output is nsamples*P where nsamples=njumps/sampleevery
%
% Example:
%   % define a forward model (here y=a*exp(-bx))
%   myfun=@(m,x)(exp(-m(1)*x)+m(2));
%   % generate some noisy data
%   true_m = [1;2];
%   x=linspace(1,10,100);
%   y=myfun(true_m,x) + .05*randn(1,100);
%   % estimate parameters
%   mt=[1;2]; % you can get mt using nonlinear opt
%   samples=mh(y,mt,@(m)(myfun(m,x)));
%   figure,plot(samples)
%
% OPTIONAL INPUT PARAMETERS:
%   LB = lower bounds on m (default=[-inf]*ones(size(x0)))
%   UB = upper bounds on m (default=[+inf]*ones(size(x0)))
%   params.burnin = #burnin iterations (default=1000)
%   params.njumps = #jumps (default=5000)
%   params.sampleevery = keep every n jumps (default=10)
%   params.update = update proposal every n jumps (default=20)
%
% AUTHOR:
%   S. Jbabdi 01/12


[y,x0,genfunc,LB,UB,params]=parseinputargs(varargin,nargin);
N=length(y);

p=x0;
s=genfunc(x0);
e=N/2*log(norm(y-s));

acc=zeros(1,length(p));
rej=zeros(1,length(p));
keepp=zeros(params.njumps,length(p));
prop=ones(1,length(p));
for i=1:params.njumps+params.burnin
    for k=1:length(p)
        oldp=p;
        p(k)=p(k)+randn*prop(k);
        if(p(k)<LB(k) || p(k)>UB(k))
            p(k)=oldp(k);
            rej(k)=rej(k)+1;
        else
            s=genfunc(p);
            olde=e;
            e=N/2*log(norm(y-s));
            if(exp(olde-e)>rand)
                acc(k)=acc(k)+1;
            else
                p(k)=oldp(k);
                rej(k)=rej(k)+1;
                e=olde;
            end
        end
    end
    keepp(i,:)=p;
    if(rem(i,params.update)==0)
        prop=prop.*sqrt((1+acc)./(1+rej));
        acc=0*acc;
        rej=0*rej;
    end
end
samples=keepp(params.burnin+1:params.sampleevery:end,:);

function [y,x0,genfunc,LB,UB,params]=parseinputargs(varargin,nargin)
y=varargin{1};
x0=varargin{2};
genfunc=varargin{3};
if(nargin>3);LB=varargin{4};else LB=-inf*ones(size(x0));end
if(nargin>4);UB=varargin{5};else UB=+inf*ones(size(x0));end
if(nargin>5);
    params=varargin{6};
else
    params=struct('burnin',1000,'njumps',5000,'sampleevery',10,'update',20);
end
if(~isfield(params,'burnin'));params.burnin=1000;end
if(~isfield(params,'njumps'));params.njumps=5000;end
if(~isfield(params,'sampleevery'));params.sampleevery=10;end
if(~isfield(params,'update'));params.update=20;end




