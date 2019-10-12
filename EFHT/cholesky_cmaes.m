function cholesky_cmaes(strfitnessfct,stopfitness,psigma,VRmin,VRmax)
% CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
% nonlinear function minimization.
%
% This code is an excerpt from cmaes.m and implements the key parts
% of the algorithm. It is intendend to be used for READING and
% UNDERSTANDING the basic flow and all details of the CMA *algorithm*.
% Computational efficiency is sometimes disregarded.

% -------------------- Initialization --------------------------------

% User defined input parameters (need to be edited)
%strfitnessfct = 'rosenbrock'; % name of objective/fitness function
sample_number=100;
interval=50;
%N = 20; % number of objective variables/problem dimension

for N=5:5:30
    %****** initialize container ******
    vals=[];
    average_gains=[];
    arz=[];
    arx=[];
    
    sigma = psigma*(VRmax-VRmin); % coordinate wise standard deviation (step-size)
    %stopfitness = 1e-10; % stop if fitness < stopfitness (minimization)
    stopeval = 50000; % stop after stopeval number of function evaluations

    % Strategy parameter setting: Selection
    lambda=12;  % population size, offspring number
    mu=6;    % number of parents/points for recombination
    mat_name=['original_population_' strfitnessfct '_' num2str(N) '.mat'];
    load(mat_name);
    %xmean = rand(N,1)*(VRmax-VRmin)+VRmin;%xmean = rand(N,1); % objective variables initial point
    %save(mat_name,'xmean');
    weights = log(mu+1/2)-log(1:mu)'; % muXone recombination weights
    weights = weights/sum(weights); % normalize recombination weights array

    A=eye(N);
    p_c=zeros(N,1); 
    p_s=zeros(N,1); % evolution paths for C and sigma

    mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu
    cs = (mueff+2)/(N+mueff+5); % t-const for cumulation for sigma control
    damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
    cc = (4+mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
    c1 = 2 / ((N+1.3)^2+mueff); % learning rate for rank-one update of C
    cmu = 2 * (mueff-2+1/mueff) / ((N+2)^2+2*mueff/2); % and for rank-mu update
    chiN=N^0.5*(1-1/(4*N)+1/(21*N^2)); % expectation of

    % -------------------- Generation Loop --------------------------------
    best_val=feval(strfitnessfct,xmean);
    counteval = 1; % the next 40 lines contain the 20 lines of interesting code
    while counteval < stopeval
       sentinel=(mod(counteval,interval*lambda)==1);
       gains=[];      
       for s=1:sample_number
           % Generate and evaluate lambda offspring
           for k=1:lambda,
               arz(:,k) = randn(N,1); % standard normally distributed vector
               arx(:,k) = xmean + sigma * (A * arz(:,k)); % add mutation % Eq. 40
               % fix x into bounds
               idl=arx(:,k)<VRmin;
               arx(idl,k)=VRmin;
               idu=arx(:,k)>VRmax;
               arx(idu,k)=VRmax;
               
               arfitness(k) = feval(strfitnessfct, arx(:,k)); % objective function call          
           end
           [arfitness, arindex] = sort(arfitness); % minimization

           if sentinel==0
               break;
           end
           
           %****** collect gains ******
            if arfitness(1)<best_val
                gains=[gains,best_val-arfitness(1)];
            else
                gains=[gains,0];
            end
       end%****** sampling ******
       
       %****** record history best val ******
       if sentinel==1
           vals=[vals,best_val];
           average_gains=[average_gains,mean(gains)];      
       end
       
       if arfitness(1)<best_val
           best_val=arfitness(1);
       end
       
       counteval = counteval+lambda;
       % Sort by fitness and compute weighted mean into xmean
       
       xmean_pre=xmean;  %record previous xmean
       xmean = arx(:,arindex(1:mu))*weights; % recombination % Eq. 42
       p_c = (1-cc)*p_c + sqrt(cc*(2-cc)*mueff) * (xmean-xmean_pre)/sigma; % Eq. 45

       A=sqrt(1-c1-cmu)*A;
       A=rank_one_update(A,c1,p_c);

       for i=1:mu
           A=rank_one_update(A,cmu*weights(i),(arx(:,i)-xmean_pre)/sigma);
       end

       % Cumulation: Update evolution paths
       p_s = (1-cs)*p_s + (sqrt(cs*(2-cs)*mueff)) *inv(A)*(xmean-xmean_pre)/sigma; % Eq. 43

       % Adapt step-size sigma
       sigma = sigma * exp((cs/damps)*(norm(p_s)/chiN - 1)); % Eq. 44

       % Break, if fitness is good enough
       if arfitness(1) <= stopfitness
           break;
       end

       disp([num2str(counteval) ': ' num2str(arfitness(1))]);

    end % while, end generation loop

        % -------------------- Final Message ---------------------------------

        disp([num2str(counteval) ':' num2str(arfitness(1))]);
        % Return best point of last generation.
        
        %****** output ******
        dir_name=['collect_' strfitnessfct];
        mkdir(dir_name);
        val_filename=[dir_name '/val_N' num2str(N) '.txt'];
        fid=fopen(val_filename,'wt');
        fprintf(fid,'%g\n',vals);
        fclose(fid);
        average_gains_filename=[dir_name '/average_gains_N' num2str(N) '.txt'];
        fid=fopen(average_gains_filename,'wt');
        fprintf(fid,'%g\n',average_gains);
        fclose(fid);
        %process data
        data_manager(dir_name,vals,average_gains,N);
end