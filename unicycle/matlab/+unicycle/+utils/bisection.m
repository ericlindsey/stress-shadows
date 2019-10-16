function x=bisection(lb,ub,fun,accuracy)
% Function BISECTION(lb,ub,fun,accuracy) finds the root of the 
% function @fun within the bounds [lb,ub] using the bisection 
% method. Accuracy sets the precision of the result.
% The function fun is assumed monotonically increasing.

import unicycle.utils.bisection

x=(lb+ub)/2;
y=fun(x);

if abs(y)<accuracy
    return;
else
    n=y>0;
    p=~n;
    x=bisection(x.*p+lb.*n,ub.*p+x.*n,fun,accuracy);
end

end