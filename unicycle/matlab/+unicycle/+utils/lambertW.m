function w = lambertW(x,branch)
% lambertW  Functional inverse of x = w*exp(w).
% w = lambertW(x), same as lambertW(x,0)
% w = lambertW(x,0)  Primary or upper branch, W_0(x)
% w = lambertW(x,-1)  Lower branch, W_{-1}(x)
%
% See: http://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-function/

% Copyright 2013 The MathWorks, Inc.

% Effective starting guess
if nargin < 2 || branch ~= -1
   w = ones(size(x));  % Start above -1
else  
   w = -2*ones(size(x));  % Start below -1
end
v = inf*w;

% Halley's method
while any(abs(w - v)./abs(w) > 1.e-8)
   v = w;
   e = exp(w);
   f = w.*e - x;  % Iterate to make this quantity zero
   w = w - f./((e.*(w+1) - (w+2).*f./(2*w+2)));
end
