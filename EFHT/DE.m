function xmin=DE(strfitnessfct,stopfitness,psigma,Lbound,Ubound,Cr,F)
% x0 is the fitness difference of the initial solution

%------------- Collection setting -----------------
sample_num=100;
interval=50;

%------------- Parameter setting ---------------
lambda=10;   %the number of offsprings
NP = 10;
maxeval=100000;
%epsilon=1e-10;   %terminate criterion
%strfitnessfct='ar_oriented_hyperellipsoid';    % name of fitness benchmark function

for N=5:5:30     % dimension


    %------------- Initialization --------------
%     x=ones(N,1)*x0;   %first generation
    x = []
    for i = 1:NP
        x= [x,Lbound + (Ubound-Lbound).*rand(N,1)];
    end
        %f=feval(strfitnessfct,x,N);

    %****** initialize collecting container ******
    fitness_difference=[];
    min_fitness_difference=1e9;
    val=[];
    average_gains=[];

        
    arf(1)=feval(strfitnessfct,x(:,1));
    
    counteval=1;
    best_f=arf(1);
    while (arf(1)>stopfitness)&& (counteval<maxeval)
       gains=[];         % initialize the container for gains
       
       % sentinel是一个标志，在每一个interval=50代的时候取ave gain
       sentinel=(mod(counteval,lambda*interval)==1);

       % 这个for循环相当于一代
        for s=1:sample_num
            % Generate and evaluate lambda offspring 遍历该代种群个体
            u = [];
            for l=1:NP
                r1 = randi(NP);
                while(r1 == l)
                    r1 = randi(NP);
                end
                r2 = randi(NP);
                while(r2 == l | r2 == r1)
                    r2 = randi(NP);
                end
                r3 = randi(NP);
                while(r3 == l | r3 == r1 | r3 == r2)
                    r3 = randi(NP);
                end
                ui = x(:,r1) + F .* (x(:,r2) - x(:,r3));
                for j = 1:size(ui,1)
                    if(rand()<Cr | j == randi(size(ui,1)))                   
                    else
                        ui(j) = x(j,l);
                    end
                end
                
                if(feval(strfitnessfct,x(:,l))>feval(strfitnessfct,ui))
                    x(:,l) = ui;
                end
                arf(l)=feval(strfitnessfct,x(:,l));
%                 min_fitness_difference=min(min_fitness_difference,max(0,arf(l)-0));
            end   %end lambda
            
             % Sort by fitness
            [arf,arindex]=sort(arf);
            
            % 如果这个interval已经结束，就要跳出arf的计算，开始计算这一整个interval的average gain，
            if sentinel==0
                break; 
            end
            
            % 如果这个interval还没结束，就要继续收集gain
            %****** collect gain and fitness_difference ******
            if arf(1)<best_f
                
                gains=[gains,best_f-arf(1)];
            else
                gains=[gains,0];
            end
            
        end %****** end sn ******

        %collect best fitness value and average gain
        % 如果这一代结束的时候正好这个interval也结束了，就计算average gain
        if sentinel==1
%             if(best_f>5000)
%                 val=[val,5000];
%             else
            fitness_difference = [fitness_difference,min_fitness_difference];
            val=[val,best_f];
%             end

            average_gains=[average_gains,mean(gains)];

        end
        %record best fitness value
        % 每一代结束之后，更新best_f
        if arf(1)<best_f
            best_f=arf(1);
        end
        
        counteval=counteval+lambda;


        disp([num2str(counteval) ':' num2str(best_f)]);
    end   %end iterative process 

    %****** output ******
    dir_name=['collect_' strfitnessfct];
    mkdir(dir_name);
    val_filename=[dir_name '/val_N' num2str(N) '.txt'];
    fid=fopen(val_filename,'wt');
    fprintf(fid,'%g\n',val);
    fclose(fid);
    val_filename=[dir_name '/fitness_difference_N' num2str(N) '.txt'];
    fid=fopen(val_filename,'wt');
    fprintf(fid,'%g\n',fitness_difference);
    fclose(fid);
    average_gains_filename=[dir_name '/average_gains_N' num2str(N) '.txt'];
    fid=fopen(average_gains_filename,'wt');
    fprintf(fid,'%g\n',average_gains);
    fclose(fid);
    %process data
    
    %这句不知道为什么就是报错
%     data_manager(dir_name,val,average_gains,N);
    
    xmin=x;
end % end N



