function p = f_test_menke(chi2A,nA,chi2B,nB)
    % return the P value for the F-test
    % reference: Menke p. 112
    %
    % E. Lindsey and R. Mallick, June 2019

    Fobs = (chi2A/nA) / (chi2B/nB);
    if( Fobs<1)
        Fobs = 1/Fobs;
    end
    p = 1 - (fcdf(Fobs,nA,nB)-fcdf(1/Fobs,nA,nB));
    
end