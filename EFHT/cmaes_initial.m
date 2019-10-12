%CMA-ES for non-linear function minization
function xmin=cmaes_initial(strfitnessfct,stopfitness,psigma,x0)
    %****** collect parameter setting ******
    sample_number=100;
    interval=50;

    for N = 5:5:30
        %****** initialize collecting container ******
        val=[];
        average_gains=[];
        arz=[];
        arx=[];
        
        %Set dimension, fitness fct, stop criteria, start value¡­
        %strfitnessfct='fdiffpow';
        maxeval=100000;
        %stopfitness=1e-10;% stop criteria
        x=ones(N,1)*x0; % object parameter start poing (weighted mean)
        sigma=psigma;
        minsigma=1e-10;% step size, minimal step size

        %parameter setting: selection
        lambda=10;
        mu=1;
        arweights=log((lambda+1)/2)-log(1:mu)'; 
        % parameter setting: adaptation
        cc=4/(N+4);
        ccov=2/(N+2^0.5)^2;
        cs=4/(N+4);
        damp=1/cs+1;

        % Initialize dynamic strategy parameters and constants
        B=eye(N);
        D=eye(N);
        BD=B*D;
        C=BD*transpose(BD);
        pc=zeros(N,1);
        ps=zeros(N,1);
        cw=sum(arweights)/norm(arweights);
        chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));

        % Generation loop
        arfitness(1)=feval(strfitnessfct,x);
        counteval=1;
        best_val=arfitness(1);
        while arfitness(1) > stopfitness && counteval < maxeval
            sentinel=(mod(counteval,interval*lambda)==1);
            gains=[];   %****** container for gains ******
            
            for i=1:sample_number
                % Generate and evaluate lambda offspring
                for k=1:lambda
                    %repeat the next two lines until arx(:,k) is feasible
                    arz(:,k)=randn(N,1);
                    arx(:,k)=x+sigma*(BD*arz(:,k));     % Eq. (13)
                    arfitness(k)=feval(strfitnessfct,arx(:,k));                   
                end

                % Sort by fitness and compute weighted mean
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
            
            end %****** end sn ******
            
            %****** record history best val ******
            if sentinel==1
                val=[val,best_val];
                average_gains=[average_gains,mean(gains)];
            end
            if arfitness(1)<best_val
                best_val=arfitness(1);
            end               
            
            counteval=counteval+lambda;
            x=arx(:,arindex(1));
            z=arz(:,arindex(1));

            % Adapt covariance matrix
            pc = (1-cc)*pc + ( sqrt(cc*(2-cc)) *cw) * (BD*z);      % Eq. (14)
            C = (1-ccov)*C+ccov*pc*transpose(pc);                              % Eq. (15)
            %adapt sigma
            ps = (1-cs)*ps + ( sqrt(cs*(2-cs)) *cw) * (B*z);          % Eq. (16)
            sigma = sigma * exp((norm(ps)-chiN)/chiN/damp);            % Eq.(17)

            % Update B and D from C
            if mod(counteval/lambda, N/10) <1
                C=triu(C)+transpose(triu(C,1));   % enforce sysmmetry
                [B,D] = eig(C);
                %limit condition of C to 1e14 + 1
                if max(diag(D)) > 1e14*min(diag(D))
                    tmp=max(diag(D)/1e14 - min(diag(D)));
                    C=C +tmp*eye(N);
                    D=D +tmp*eye(N);
                end
                D=diag(sqrt(diag(D)));  % D contains standard deviations now
                BD=B*D; % for speed up only
            end % if mod

            % Adjust minimal step size
            if (sigma*min(diag(D)) < minsigma) | (arfitness(1) == arfitness(min(mu+1,lambda))) | (x==x+0.2*sigma*BD(:,1+floor(mod(counteval/lambda,N))))
                sigma = 1.4*sigma;
            end
        end % while,, end generation loop 

        disp([num2str(counteval) ':' num2str(arfitness(1))]);
        xmin = arx(:,arindex(1)); %return best point of last generation
    
        %****** output ******
        dir_name=['collect_' strfitnessfct];
        mkdir(dir_name);
        val_filename=[dir_name '/val_N' num2str(N) '.txt'];
        fid=fopen(val_filename,'wt');
        fprintf(fid,'%g\n',val);
        fclose(fid);
        average_gains_filename=[dir_name '/average_gains_N' num2str(N) '.txt'];
        fid=fopen(average_gains_filename,'wt');
        fprintf(fid,'%g\n',average_gains);
        fclose(fid);
        %process data
        data_manager(dir_name,val,average_gains,N);
        
    end %****** end N******
    
        




