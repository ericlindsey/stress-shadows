function [train,test]=split_jointinv_traintest(scenario,Itest)
    % split the dataset from scenario0 into two new objects for training and
    % testing

    % first we have to create independent objects, because matlab copy uses
    % pointers
    train=Jointinv(scenario.expNumber);
    train.run_setup();
    test =Jointinv(scenario.expNumber);
    test.run_setup();
    
    %generate two disjoint lists of IDs - length Ntest and length-Ntest
    Itrain=1:length(scenario.dataVector);
    Itrain(Itest)=[];
    delete_jointinv_data(train,Itest);
    delete_jointinv_data(test,Itrain);

end