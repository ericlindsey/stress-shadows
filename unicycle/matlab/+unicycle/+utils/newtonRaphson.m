function x = newtonRaphson(f,dfdx,x0,tolerance)
% function NEWTONRAPHSON performs a gradient descent to find the zero
% of a target function given the function itself and its derivative.

err=Inf;
x=x0;
while max(abs(err))>tolerance
    xPrevious=x;
    x=xPrevious-f(xPrevious)./dfdx(xPrevious);
    err=x-xPrevious;
end

end