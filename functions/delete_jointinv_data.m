function delete_jointinv_data(scenario0,Idel)
    % delete data from the jointinv object, and design matrix

    %scenario0.datasets{1}.name(Idel)=[]; % name deletion is broken
    %scenario0.datasets{1}.coordinates(Idel,:)=[];
    %scenario0.datasets{1}.lat(Idel)=[];
    %scenario0.datasets{1}.lon(Idel)=[];
    
    % warning! The current implementation drops individual lines from the
    % datavector, but this has little physical meaning
    
    %if length(scenario0.userParams.dataComponents)==1
        scenario0.dataVector(Idel)=[];
        scenario0.designMatrix(Idel,:)=[];
        scenario0.datasets{1}.dataVector(Idel)=[];
        scenario0.datasets{1}.covarianceMatrix(Idel,:)=[];
        scenario0.datasets{1}.covarianceMatrix(:,Idel)=[];
        
%     elseif  length(scenario0.userParams.dataComponents)==2
%         warning('2d case not implemented')
%     
%     else %length(scenario0.userParams.dataComponents)==3 assumed
%         warning('3D case not implemented')
%     
%     end
    
end
