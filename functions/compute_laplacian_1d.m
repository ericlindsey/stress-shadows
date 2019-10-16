function L = compute_laplacian_1d(N)
    % for a 1-d fault with uniform spacing, give a simple 2nd-order finite
    % difference operator for the laplacian
    
    e = ones(N,1);
    L = eye(N)*spdiags([e -2*e e], -1:1, N, N);
    % fix the edges - 2nd order finite difference operator
    L(1,:)=L(2,:);
    L(end,:)=L(end-1,:);
    
end