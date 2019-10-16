% Copula approach for optimized calculation of posterior distribution of model
% parameters. In particular, here we use Gaussian copulas because of their 
% tractability. For more details, refer UNICYCLE documentation. 
%
% Note : 
% For details on how to use this check ../../examples/tutorials/copula_tutorial 
% directory.
%
% AUTHOR 
% Original : Yu Hang, NTU, Jun. 2016
% Modified : Sagar Masuti, Dec 2016
% Modified : Dongju Peng , Mar 2017

classdef copula 
    properties
        % Number of forward models computed/to be computed.
        Ns;                      
        % Number of model parameters.
        M;                         
        % Bounds on the model parameters.
        bounds;    
        % Type of prior distribution. (check constant Properties below) 
        prior_dist_type;           
        % Marginal estimation method. (check constant Properties below) 
        marginal_est_method;      
        % Whether to compute the forward model or its already computed.
        bforward;
        % Whether the noise is constant or variable data dependent
        stdnoise
        % constant std noise (if not given will be calculated internally) 
        v; 
        % model location
        location;
        % gps name
        gps;
        % Number of cores to use for parallel computing.
        nproc;
    end
    properties (Constant)
        % Prior distribution type.
        uniform=1;
        lognormal=2;
        normal=3;

        % Marginal estimation method.
        kernel_density_method=1;
        gaussian_mixture_method=2;
        empirical_method=3;
    end
    
    methods (Static)
        function s=getline(fid)
            nc=0;
            while (0 == nc)
                t=fgetl(fid);
                if ('#' == t(1))
                    continue;
                else
                    nc=1;
                    s=t;
                end
            end
        end
    end
 
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        % Two options to create the instance of this class. 
        % 1) copobj=copula();
        %    copobj.Ns=100;
        %    copobj.M=2;
        %    ...
        %    ...
        %    copobj.nproc=2;
        % 
        % 2) copobj=copula('/path/to/config/file', 'config filename');
        % ---------------------------------------------------------------------
        function obj=copula(varargin)

            obj.v=0;
            if isempty(varargin)
                return
            end
            
            if (2>nargin)
                error('unicycle.optim.copula: Not enough input argument');
            end

            filename=strcat(varargin{1}, varargin{2});
            fid=fopen(filename, 'r');
            
            obj.Ns=str2num(obj.getline(fid));            
            fprintf('# Number of forward models: \n%d\n',obj.Ns)
            
            obj.M=str2num(obj.getline(fid));
            fprintf('# Number of model parameters: \n%d\n',obj.M);
            
            fprintf('# Bounds on the model parameters:\n');
            for i=1:obj.M
                bds=str2num(obj.getline(fid));
                obj.bounds=[obj.bounds;bds];
                fprintf('%.5f\t%.5f\n',bds(1),bds(2));
            end
            
            obj.prior_dist_type=str2num(obj.getline(fid));
            fprintf('# Type of prior distribution. (1 -> Uniform, 2 -> Lognormal, 3 -> Normal): \n%d\n',...
                    obj.prior_dist_type);
                
            obj.marginal_est_method=str2num(obj.getline(fid));
            fprintf('# Marginal estimation method. (1-> Kernel density method, 2-> Gaussian mixture method): \n%d\n',...
                    obj.marginal_est_method);

            obj.bforward=str2num(obj.getline(fid));
            fprintf('# Whether to compute the forward model or not: \n%d\n',...
                     obj.bforward);
                 
            obj.stdnoise=str2num(obj.getline(fid));
            fprintf('# whether the noise is constant or data dependent if it is considerd: \n%d\n',...
                     obj.stdnoise);

            obj.location=obj.getline(fid);
            fprintf('# Model location: \n%s\n', obj.location);
            
            obj.gps=obj.getline(fid);
            fprintf('# GPS network file location: \n%s\n', obj.gps);         
            
            obj.nproc=str2num(obj.getline(fid));
            fprintf('# Number of cores to use for parallel computing : \n%d\n',...
                     obj.nproc);

            fclose(fid);
        end

        function cdf = empcdf_con(o,Data)
            % Xu Shiyan, Yu Hang, NTU
            % convert pdf to its corresponding cdf
            [sortedData oldIndex] = sort(Data);%sort in ascending order
            cdf = ones(size(Data));
            for i=1:size(cdf,2);%1~25
                m = 1;
                n = 1;
                for j=2:size(cdf,1)%2~250
                      if sortedData(j,i)==sortedData(j-1,i)
                          n = n + 1;
                          continue;
                      else
                          for k=m:n
                              cdf(oldIndex(k,i),i) = j-1;
                          end
                          n = n + 1;
                          m = n;
                      end
                end
                for j=m:n
                    cdf(oldIndex(j,i),i) = size(Data,1);
                end
            end
            cdf = (cdf-0.5)/size(sortedData,1);
        end

        
        function [mp,Sm,ds]=posterior(obj,d)
            % Computes the posterior mean and covariance of conditional 
            % distribution p(m|dobs). This method uses the quasi Monte Carlo 
            % sampling to sample the prior distribution and calls a function 
            % called "forward". The forward function is the place to provide 
            % the physical model. Its the users responsibility to have this 
            % forward function in the path so this method can invoke it.
            % 
            % INPUT:
            % d - Observed data
            % 
            % OUTPUT:
            % mp - Posterior mean of the conditional distribution p(m|d)
            % Sm - Posterior covariance of the conditional distribution p(m|d)
            % 
            % Invoking:
            % copobj=copula('config/file/path','configfile');
            % [mp,sm]=copobj.posterior();
            % Also refer to example in "../../examples/tutorials/copula_tutorial"
            % ------------------------------------------------------------------


            %% ------------------------- Configurations ------------------------
            Ns=obj.Ns;                      
            M=obj.M;                        
            bounds=obj.bounds;              
            marginal_est_method=obj.marginal_est_method;      

            mu_prior=1*ones(1,M); 
            sgm_prior=3*ones(1,M);

            %% ------------------ Compute/load synthetic data ------------------
            N=size(d,1);               % Number of data points in observed data.            
            d(isnan(d))=0;
            if (obj.bforward==1)
                % Quasi Monte Carlo Sampling
                urnd=haltonset(M,'Skip',1e3,'Leap',1e2);
                urnd=scramble(urnd,'RR2');
                urnd=qrandstream(urnd);
                msCDF=qrand(urnd,Ns);

                ms=zeros(Ns,M);
                switch obj.prior_dist_type
                case obj.uniform
                    for i=1:M
                        ms(:,i)=bounds(i,1)+(bounds(i,2)-bounds(i,1))*msCDF(:,i);
                    end
                case obj.lognormal
                    for i=1:M
                        p_bounds=logncdf(bounds(i,:),0,3);
                        ms(:,i)=logninv(p_bounds(i,1) + msCDF(:,i).*diff(p_bounds(i,:)),0,3);
                    end
                case obj.normal
                    for i=1:M
                        p_bounds(i,:)=normcdf(bounds(i,:),mu_prior(i),sgm_prior(i));
                        ms(:,i)=norminv(p_bounds(i,1)+msCDF(:,i).*diff(p_bounds(i,:)),mu_prior(i),sgm_prior(i));
                    end
                end

                % learn conditional density p(d_i|m)
                % compute the values of d_i given samples of m
                ds=[];
                for i=1:Ns
                    tic
                    minput=ms(i,:);
                    disp(sprintf ('running forward model: %d',i));
                    a=forward(minput).';
                    ds=[ds;a];
                    toc
                end
            else
                load ds;
                load ms;
                msCDF=zeros(Ns,M);
                
                switch obj.prior_dist_type
                    case obj.uniform
                        for i=1:M
                            msCDF(:,i)=(ms(:,i)-bounds(i,1))/(bounds(i,2)-bounds(i,1));
                        end
                    case obj.lognormal
                        for i=1:M
                            p_bounds(i,:)=logncdf(bounds(i,:),0,3);
                            msCDF(:,i)=(logncdf(ms(:,i),0,3)-p_bounds(i,1))/diff(p_bounds(i,:));
                        end
                    case obj.normal
                        for i=1:M
                            mu_prior(1,i)=sum(bounds(i,:))/2;
                            sgm_prior(1,i)=1;
                            p_bounds(i,:)=normcdf(bounds(i,:),mu_prior(i),sgm_prior(i));
                            msCDF(:,i)=(normcdf(ms(:,i),mu_prior(i),sgm_prior(i))-p_bounds(i,1))/diff(p_bounds(i,:));
                        end
                end
            end
            %%
            ds(isnan(ds))=0;
            lb_ds=min(ds,[],1); % lower bound
            ub_ds=max(ds,[],1); % upper bound


            id1=find(d.'>=lb_ds & d.'<=ub_ds);    % d(id1) are within the bound of ds(:,id1)
            id2=find(d.'==lb_ds & d.'==ub_ds);    % ds(:,id2) are deterministic rather than random
            id3=setdiff(1:N,id1);                 % d(id3) are outside the bound of ds(id3);
            id1=setdiff(id1,id2);

            n1=length(id1);
            n2=length(id2);
            n3=length(id3);

            % reorder ds and d to reduce necessary communication during parallel computing
            ds=[ds(:,id1),ds(:,id2),ds(:,id3)];
            d=[d(id1);d(id2);d(id3)];

            pool=parpool(obj.nproc);

            dCDF=zeros(1,N);

            Nss=100*Ns;
            msg=norminv(msCDF);           % convert non-Gaussian samples to Gaussia
            msgs=kron(ones(100,1),msg);

            Sc=zeros(1,1,N);              % conditional covariance
            Sr=zeros(1,M,N);              % Sdm/Smm

            Smm=msg.'*msg/Ns; %eye(M);   % independent prior for each component of m
            maxiter=1000;

            % --------------------- Marginal density estimation ----------------------
            % determine whether the noise is constant or data dependent
            if obj.stdnoise==1  
                % A single noise standard deviation for all data points
                if (0==obj.v)
                    if any(d.'<lb_ds | d.'>ub_ds)
                       v=max(min(abs(repmat(d.',Ns,1)-ds)))/3;
                    else
                       v=min(std(repmat(d.',Ns,1)-ds,[],2));
                    end
                else 
                    v=obj.v;
                end
                
                parfor i=1:N
                    dsi=ds(:,i);
                    dss=kron(ones(100,1),dsi)+v*randn(Nss,1);

                    %if (obj.kernel_density_method==marginal_est_method)
                    if (1==marginal_est_method)
                        % kernel density estimation method.
                        dsCDF=0;
                        for j=1:Ns
                            dsCDF=dsCDF+normcdf(dss,dsi(j),v);
                        end
                        dsCDF=dsCDF/Ns;
                        dCDF(i)=mean(normcdf(d(i),dsi,v));
                    elseif (2==marginal_est_method)
                        % Gaussian mixture model
                        k=min(Ns/20,32); % no. of mixtures
                        options=statset('MaxIter',maxiter);
                        o=gmdistribution.fit(dss,k,'Regularize',1e-9,'Start','plus','Options',options); %,'Replicates',1
                        dsCDF=cdf(o,dss);
                        dCDF(i)=cdf(o,d(i));
                        if dCDF(i)==1
                            dCDF(i);
                        end
                    else
                        % empirical CDF
                        CDF=obj.empcdf_con([dss;d(i)]);
                        dsCDF=CDF(1:end-1);
                        dCDF(i)=CDF(end);
                    end
                    dsg=norminv(dsCDF);
                    Sdm=dsg.'*msgs/Nss;
                    Sr(:,:,i)=Sdm/Smm;
                    Sc(:,:,i)=dsg.'*dsg/Nss - Sdm/Smm*Sdm.';
                end
            else
                % different noise standard deviation for different data points
                v0=min(std(repmat(d(1:n1).',Ns,1)-ds(:,1:n1),[],2));
                parfor i=1:N
                    if i > n1 && i <=n1+n2
                        dCDF(i)=0.5;
                        dsg=randn(Nss,1);
                        Sdm=dsg.'*msgs/Nss;
                        Sr(:,:,i)=Sdm/Smm;
                        Sc(:,:,i)=dsg.'*dsg/Nss - Sdm/Smm*Sdm.';
                    else
                        dsi=ds(:,i);
                        if i <=n1
                            v=v0;
                        else
                            v=min(abs(dsi-d(i)))/3;
                        end
                        
                        dss=kron(ones(100,1),dsi)+v*randn(Nss,1);
                        %if (obj.kernel_density_method==marginal_est_method)
                        if (1==marginal_est_method)
                            % kernel density estimation method.
                            [~,~,bw]=ksdensity(dsi);
                            bw=max(bw,v);
                            dsCDF=0;
                            for j=1:Ns
                                dsCDF=dsCDF+normcdf(dss,dsi(j),bw);
                            end
                            dsCDF=dsCDF/Ns;
                            dCDF(i)=mean(normcdf(d(i),dsi,bw));
                        %elseif (obj.gaussian_mixture_method==marginal_est_method)
                        elseif (2==marginal_est_method)
                            % Gaussian mixture model
                            if i <=n1
                                k=min(Ns/20,32);
                            else
                                k=min(Nss/20,32); % no. of mixtures
                            end
                            options=statset('MaxIter',maxiter);
                            o=gmdistribution.fit(dss,k,'Regularize',1e-9,'Start','plus',...
                                'Options',options); %,'Replicates',1
                            dsCDF=cdf(o,dss);
                            dCDF(i)=cdf(o,d(i));
                            if dCDF(i)==1
                                dCDF(i);
                            end
                        else
                            % empirical CDF
                            CDF=obj.empcdf_con([dss;d(i)]);
                            dsCDF=CDF(1:end-1);
                            dCDF(i)=CDF(end);
                        end
                        dsg=norminv(dsCDF);
                        Sdm=dsg.'*msgs/Nss;
                        Sr(:,:,i)=Sdm/Smm;
                        Sc(:,:,i)=dsg.'*dsg/Nss - Sdm/Smm*Sdm.';
                    end
                end
            end

            delete(gcp);

        % --------------------- Posterior distribution p(m|d) -------------------------
            % compute the CDF of data d
            dg=norminv(dCDF);

            % compute the posterior Gaussian copula
            Km=eye(M);        % prior preicision matrix  (assume all components of m are independent in p(m))
            hm=0;             % prior potential vector
            for i=1:N
                Km=Km+Sr(:,:,i).'/Sc(:,:,i)*Sr(:,:,i);
                hm=hm+Sc(:,:,i)\Sr(:,:,i).'*dg(i);
            end
            mp=Km\hm;         % posterior mean
            Sm=inv(Km);       % posterior covariance
        end % end posterior
    end % end methods
end % end class
