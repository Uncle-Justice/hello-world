function xmin=ES(strfitnessfct,stopfitness,psigma,x0)
% x0 is the fitness difference of the initial solution

%------------- Collection setting -----------------
sample_num=100;
interval=50;

%------------- Parameter setting ---------------
lambda=10;   %the number of offsprings
maxeval=100000;
%epsilon=1e-10;   %terminate criterion
%strfitnessfct='ar_oriented_hyperellipsoid';    % name of fitness benchmark function

for N=5:5:30     % dimension
    c=1/sqrt(N);
    cu=sqrt(c/(2-c));
    beta=1/sqrt(N);
    beta_scal=1/N;

    %------------- Initialization --------------
    sigma=psigma;     % step side
    sigma_scal=ones(N,1);
    z=zeros(N,1);
    x=ones(N,1)*x0;   %first generation
    %f=feval(strfitnessfct,x,N);

    %****** initialize collecting container ******
    val=[];
    average_gains=[];
    arz=[];
    arx=[];
        
    arf(1)=feval(strfitnessfct,x);
    counteval=1;
    best_f=arf(1);
    while (arf(1)>stopfitness)&& (counteval<maxeval)
       gains=[];         % initialize the container for gains
       
       % sentinel是一个标志，在每一个interval=50代的时候取ave gain
       sentinel=(mod(counteval,lambda*interval)==1);

       % 这个for循环相当于一代
        for s=1:sample_num
            % Generate and evaluate lambda offspring 遍历该代种群个体
            for l=1:lambda
                arz(:,l)=randn(N,1);
                arx(:,l)=x+sigma*sigma_scal.*arz(:,l);
                arf(l)=feval(strfitnessfct,arx(:,l));
            end   %end lambda
            
             % Sort by fitness
            [arf,arindex]=sort(arf);
            
            % 如果这个interval已经结束，就要跳出arf的计算，开始计算这一整个interval的average gain，
            if sentinel==0
                break; 
            end
            
            % 如果这个interval还没结束，就要继续收集gain
            %****** collect gain ******
            if arf(1)<best_f
                gains=[gains,best_f-arf(1)];
            else
                gains=[gains,0];
            end
            
        end %****** end sn ******

        %collect best fitness value and average gain
        % 如果这一代结束的时候正好这个interval也结束了，就计算average gain
        if sentinel==1
            val=[val,best_f];
            average_gains=[average_gains,mean(gains)];
        end
        %record best fitness value
        % 每一代结束之后，更新best_f
        if arf(1)<best_f
            best_f=arf(1);
        end
        
        counteval=counteval+lambda;
        x=arx(:,arindex(1));

        % randomly select a arz
        z=(1-c)*z+c*arz(:,arindex(1));

        b=norm(z)/sqrt(N)/cu;
        sigma=sigma*( exp(norm(z)/sqrt(N)/cu -1+1/(5*N)) )^beta;
        a=norm(z)/cu;
        sigma_scal=sigma_scal*( norm(z)*cu+0.85 )^beta_scal;

        disp([num2str(counteval) ':' num2str(best_f)]);
    end   %end iterative process 

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
    
    xmin=x;
end % end N



