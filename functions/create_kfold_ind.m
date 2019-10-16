function klists = create_kfold_ind(N,k)

    % create a k by N/k matrix with rows containing randomly shuffled
    % disjoint sets of indices, for k-fold cross-validation
    
    p=randperm(N);
    resid=mod(N,k);
    klists={}; %use cells because the lists may not all be the same length
    
    %loop method is ugly but it works
    lastI=1;
    for i=1:k
        if i<=resid
            newI = lastI + ceil(N/k);
        else
            newI = lastI + floor(N/k);
        end
        klists{i}=p(lastI:newI-1);
        lastI=newI;
    end

end