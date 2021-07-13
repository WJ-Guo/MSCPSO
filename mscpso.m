%PSO >> function for the PSO ALGORITHM
%
% USAGES:   1.) [fxmin, xmin, Swarm, history] = MSCPSO(psoOptions);
%           2.) [fxmin, xmin, Swarm, history] = MSCPSO;
%           3.) fxmin = PSO(psoOptions);
%           3.) PSO 
%           etc.
%
% Arguments     : psoOptions--> A Matlab stucture containing all PSO related options. (see also: get_psoOptions)
% Return Values : [fxmin, xmin, Swarm, history]
%                     |     |       |       |_The history of the algorithm. Depicts how the function value of GBest changes over the run.
%                     |     |       |_The final Swarm. (A matrix containing co-ordinates of all particles)
%                     |     |_The co-ordinates of Best (ever) particle found during the PSO's run.
%                     |__The objective value of Best (^xmin) particle.
%
%  History        :   Author      :   JAG (Jagatpreet Singh)
%                     Created on  :   05022003 (Friday. 2nd May, 2003)
%                     Comments    :   The basic PSO algorithm.
%                     Modified on :   0710003 (Thursday. 10th July, 2003)
%                     Comments    :   It uses psoOptions structure now. More organized.
%  
%                     see also: get_psoOptions
function [fxmin, xmin, Swarm, history] = MSCPSO(psoOptions,initSWARM,initStep)


%Globals
global psoFlags;
global psoVars;
global psoSParameters;
global notifications;


upbnd = 600; % Upper bound for init. of the swarm
lwbnd = 300; % Lower bound for init. of the swarm
GM = 0; % Global minimum (used in the stopping criterion)
ErrGoal = 1e-10; % Desired accuracy
% 

%Initializations
if nargin == 0
    psoOptions = get_psoOptions;
end


%For Displaying 
if psoOptions.Flags.ShowViz
    global vizAxes; %Use the specified axes if using GUI or create a new global if called from command window
    vizAxes = plot(0,0, '.');
    axis([-1000 1000 -1000 1000 -1000 1000]);   %Initially set to a cube of this size
    axis square;
    grid off;
    set(vizAxes,'EraseMode','xor','MarkerSize',15); %Set it to show particles.
    pause(1);
end
%End Display initialization

% Initializing variables
success = 0; % Success Flag
iter = 0;   % Iterations' counter
fevals = 0; % Function evaluations' counter

% Using params---
% Determine the value of weight change
w_start = psoOptions.SParams.w_start;   %Initial inertia weight's value
w_end = psoOptions.SParams.w_end;       %Final inertia weight
w_varyfor = floor(psoOptions.SParams.w_varyfor*psoOptions.Vars.Iterations); 
%Weight change step. Defines total number of iterations for which weight is changed.
w_now = w_start;
inertdec = (w_start-w_end)/w_varyfor; %Inertia weight's change per iteration

% Initialize Swarm and Velocity初始化的PSO的粒子数量
SwarmSize = psoOptions.Vars.SwarmSize;
if (exist('initSWARM'))

		Swarm = initSWARM*(psoOptions.Obj.ub-psoOptions.Obj.lb) + psoOptions.Obj.lb;
	else
		Swarm = rand(SwarmSize, psoOptions.Vars.Dim)*(psoOptions.Obj.ub-psoOptions.Obj.lb) + psoOptions.Obj.lb;
	end


	if (exist('initStep'))
	   VStep = initStep;
	else
	   VStep = rand(SwarmSize, psoOptions.Vars.Dim);
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%	%%%%%%%%%%%%%%%%%%%%%function parameter define.
	[SwarmSize, Dim] = size(Swarm);
   [M1,M2,shifto,lambda10,lambda100]=functionparameter(SwarmSize, Dim);

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

f2eval = psoOptions.Obj.f2eval; %The objective function to optimize.需优化的目标函数

%Find initial function values.
fSwarm = feval(f2eval, Swarm,M1,M2,shifto,lambda10,lambda100);%最小值最好
fevals = fevals + SwarmSize;

% Initializing the Best positions matrix and
% the corresponding function values
PBest = Swarm;
fPBest = fSwarm;

% Finding best particle in initial population
[fGBest, g] = min(fSwarm);
lastbpf = fGBest;
Best = Swarm(g,:); %Used to keep track of the Best particle ever
fBest = fGBest;
history = [0, fGBest];

if psoOptions.Flags.Neighbor
    % Define social neighborhoods for all the particles
    for i = 1:SwarmSize
        lo = mod(i-psoOptions.SParam.Nhood:i+psoOptions.SParam.Nhood, SwarmSize);
        nhood(i,:) = [lo];
    end
    nhood(find(nhood==0)) = SwarmSize; %Replace zeros with the index of last particle.
end

if psoOptions.Disp.Interval & (rem(iter, psoOptions.Disp.Interval) == 0)
    disp(sprintf('Iterations\t\tfGBest\t\t\tfevals'));
end


%初始阈值 0.0005(0,1)  0.005(-5.12,5.12),0.05(-50,50),0.1(-100,100),0.3(-300,300)
T_d = repmat(0.1,1,psoOptions.Vars.Dim); %形成一个行向量[0.1 0.1 0.1 0.1]表示为V的阈值
F_d = repmat(0,1,psoOptions.Vars.Dim);   %%
K_1 = 5;
K_2 = 10;

Initmutestd = (psoOptions.Obj.ub-psoOptions.Obj.lb);%初始化方差

M = 2;                %%%scale
P = SwarmSize/M;
	%初始化标准方差
	InitmutestdQ = repmat(Initmutestd,[1,M])
	SubMutesize = 5;
cloneQSize = M*SubMutesize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  THE  PSO  LOOP                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while( (success == 0) & (iter <= psoOptions.Vars.Iterations) )
    iter = iter+1;
    
    % Update the value of the inertia weight w
    if (iter<=w_varyfor) & (iter > 1)
        w_now = w_now - inertdec; %Change inertia weight
    end
    
    
    %%%%%%%%%%%%%%%%%
    % The PLAIN PSO %
    
    % Set GBest
    A = repmat(Swarm(g,:), SwarmSize, 1); %A = GBest. repmat(X, m, n) repeats the matrix X in m rows by n columns.
    B = A; %B wil be nBest (best neighbor) matrix
    
            
    % Generate Random Numbers
    R1 = rand(SwarmSize, psoOptions.Vars.Dim);
    R2 = rand(SwarmSize, psoOptions.Vars.Dim);
    
    % Calculate Velocity
    if ~psoOptions.Flags.Neighbor   %Normal
        VStep = psoOptions.SParams.c1*R1.*(PBest-Swarm) + psoOptions.SParams.c2*R2.*(A-Swarm);   %%without Wi
    else %With neighborhood
        R3 = rand(SwarmSize, psoOptions.Vars.Dim); %random nos for neighborhood
        VStep = psoOptions.SParams.c1*R1.*(PBest-Swarm) + psoOptionsSParams.c2*R2.*(A-Swarm) + psoOptionsSParams.c3*R3.*(B-Swarm);
    end
    
    			
%CMEPSO给出了一种多尺度，定向、定时的变异操作，即当速度小于一定变异操作，
    %自动分家并寻找新集聚地繁衍后代的特征，将会通过逃逸运动寻找新的生存地。这时忘记
    %自身的历史最佳位置而仅记住种群的最佳位置。			
   %  d=[1 2 3;2 1 5]
 % d =
  %     1     2     3
 %      2     1     5
 % >> d==2
 % ans =
  %     0     1     0
  %     1     0     0
 % >> find(d==2)
 % ans =
 %      2
  %     3
 % >> [a b]=find(d==2)
 % a =
 %      2
  %     1
 % b =
 %      1
 %      2 
    % Apply Vmax Operator for v > Vmax
    changeRows = VStep > psoOptions.SParams.Vmax;
    VStep(find(changeRows)) = psoOptions.SParams.Vmax;
    % Apply Vmax Operator for v < -Vmax
    changeRows = VStep < -psoOptions.SParams.Vmax;
    VStep(find(changeRows)) = -psoOptions.SParams.Vmax;   
    
    
    % ::UPDATE POSITIONS OF PARTICLES::
    Swarm = Swarm + psoOptions.SParams.Chi * VStep;    % Evaluate new Swarm
  
    
    
		   	
    %为了防止早熟收敛，搜索能力随代数增加而逐渐下降的特点，ESPSO引入随机变异操作，以便逃出局部极小点
    EsacepeRows = (abs(VStep) <= repmat(T_d,SwarmSize,1));                   %%%particles satisfy mutation condition      
     %EsacepeRows = (VStep == 0);
    [rowE,col,val]=find(EsacepeRows);                           %%%save position
		
		rowOE=unique(rowE);
		
		for i=1:length(rowOE)
				 curSwarm = Swarm(rowOE(i),:);              %%提取种群可变异的粒子
				 
				
			
				 colE = col(find(rowE==rowOE(i)));          %发现所有列中出现上述小与自己维度阈值的列
				 
				 VStepBak = VStep(rowOE(i),:);              
			
				 VStepBak1 = VStepBak;
				 if (M~=1)
				 		VStepBak1remp = repmat(VStepBak1,M+2,1);
				 else
				  	VStepBak1remp = repmat(VStepBak1,M+1,1);
				 end
				 
		     VStepBak1remp(1:M,colE) = randn(M,length(colE)).*repmat(InitmutestdQ',1,length(colE));
		     %InitmutestdQ每一个维度的方差都相同
		     
   			 VStepBak1remp(M+1,colE) = rand(1,length(colE)).*psoOptions.SParams.Vmax;	
   			 %含两种变异策略，最后一行是Vmax，均匀变异
   			 
   			 
   			 %实现各种变异并产生相应的结果
   			 if (M~=1)
   			 		VStepBak1remp(M+2,colE) = randn(1,length(colE)).*repmat(Initmutestd,1,length(colE));
		        SwarmBak = repmat(curSwarm,M+2,1) + psoOptions.SParams.Chi .* VStepBak1remp;    % Evaluate new Swarm
         else
           SwarmBak = repmat(curSwarm,M+1,1) + psoOptions.SParams.Chi .* VStepBak1remp;   % Evaluate new Swarm
         end	
   
         %检查变异后的可能出现的出界情况
         changeRows = SwarmBak > psoOptions.Obj.ub;
    		 SwarmBak(find(changeRows)) = rand*(psoOptions.Obj.ub-psoOptions.Obj.lb) + psoOptions.Obj.lb;;
         % Apply Vmax Operator for v < -Vmax
         changeRows = SwarmBak < psoOptions.Obj.lb;;
       	 SwarmBak(find(changeRows)) = rand*(psoOptions.Obj.ub-psoOptions.Obj.lb) + psoOptions.Obj.lb;;
          
          
           %SwarmBak 
           %pause
		   [swarmsize1,dim1]=size(SwarmBak);
         shiftoBak=shifto(1:swarmsize1,1:dim1);
		 fSwarmBak = feval(f2eval,SwarmBak,M1,M2,shiftoBak,lambda10,lambda100);
		 
         if (M~=1)
         		fevals = fevals + M+2;
         else
         		fevals = fevals + M+1;%记录计算目标函数的次数
         end
         
  	  	[sortfSwarm,indxfSwarm]=sort(fSwarmBak);
  	  	%最后的变异后的结果代替
  	  	Swarm(rowOE(i),:) = SwarmBak(indxfSwarm(1));
  	 end 	
  	
  	
    
    fSwarm = feval(f2eval, Swarm,M1,M2,shifto,lambda10,lambda100);
    fevals = fevals + SwarmSize;
    [sortfSwarm,findex] = sort(fSwarm);	          %%升序排列
    
    	
    	%计算方差的变化大小,将种群划分为M份，每份含SWarmSize/M个个体
				for i =1:M,
					FixXi(i)=sum(sortfSwarm((i-1)*P+1:i*P));	
				end;
				
				FixXimax = max(FixXi);
				FixXimin = min(FixXi);
				
				if (FixXimax~=FixXimin)
				    InitmutestdQ = InitmutestdQ.*exp((M*FixXi-sum(FixXi))/(M*(FixXimax-FixXimin)) );
			 	end
			 	
				%防止变换后的方差出界	
				InitmutestdQ=StdQstandard(InitmutestdQ,psoOptions.Obj.ub-psoOptions.Obj.lb);
			 
        
        
        
        B_id = zeros(SwarmSize,psoOptions.Vars.Dim);
		    B_id(find(EsacepeRows))=1;
		    %将那些v不好的所对应的swarm和其对应的列标志为1
		    
		    
		    b_SUM = sum(B_id,1);%形成行向量
		    
		    F_d = F_d + b_SUM;%记录行向量中的每一个Dim所出现小与阈值的次数，记录的是累积数
		      
		   	T_d(find(F_d > K_1))=T_d(find(F_d > K_1))./K_2;
		   	%表明你设得阈值太大，需要减少一些，
		   	F_d(find(F_d > K_1))=0; %由于已经做了调整，所以清零，重新计算累积
        
        
    
    % Updating the best position for each particle更新局部最优
    changeRows = fSwarm < fPBest;
    fPBest(find(changeRows)) = fSwarm(find(changeRows));
    PBest(find(changeRows), :) = Swarm(find(changeRows), :);
    
    lastbpart = PBest(g, :);%记录原来是最优位置现经过一轮后现在变异成的粒子值
    
    %AEPSO给出了一种定向、定时的变异操作，即当速度小于一定变异操作，
    %自动分家并寻找新集聚地繁衍后代的特征，将会通过逃逸运动寻找新的生存地。这时忘记
    
    
   % fPBest(rowOE) = fSwarm(rowOE);
   % PBest(rowOE,:) = Swarm(rowOE,:);
    
    %自身的历史最佳位置而仅记住种群的最佳位置。
    
 
    
   
    
    
    % Updating index g
    [fGBest, g] = min(fPBest);

    %Update Best. Only if fitness has improved.
    if fGBest < lastbpf
        [fBest, b] = min(fPBest);
        Best = PBest(b,:);
    end
   
   
    
    %%OUTPUT%%
    if psoOptions.Save.Interval & (rem(iter, psoOptions.Save.Interval) == 0)
        history((size(history,1)+1), :) = [iter, fBest];
    end
    
    if psoOptions.Disp.Interval & (rem(iter, psoOptions.Disp.Interval) == 0)
        disp(sprintf('%4d\t\t\t%.5g\t\t\t%5d', iter, fGBest, fevals));
    end

    if psoOptions.Flags.ShowViz
        [fworst, worst] = max(fGBest);
        DrawSwarm(Swarm, SwarmSize, iter, psoOptions.Vars.Dim, Swarm(g,:), vizAxes);
    end
    
    %%TERMINATION%%
    %if abs(fGBest-psoOptions.Obj.GM) <= psoOptions.Vars.ErrGoal     %GBest
    %    success = 1;
    %elseif abs(fBest-psoOptions.Obj.GM)<=psoOptions.Vars.ErrGoal    %Best
    %    success = 1
    %else
    %    lastbpf = fGBest; %To be used to find Best
    %end

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  END  OF PSO  LOOP                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[fxmin, b] = min(fPBest);
xmin = PBest(b, :);

history = history(:,2);
%Comment below line to Return Swarm. Uncomment to return previous best positions.
Swarm = PBest; %Return PBest


function [stdq]=StdQstandard(q,W);	
		while(find(abs(q)>W/2)>0) 
			q(find(abs(q)>W/2))=abs(W/2-q(find(abs(q)>W/2)));
		end
		stdq=q;
	
	function [stdx]=xstandard(x,a,b);	
		while(find(x>b)>0) 
			x(find(x>b))=2*b-x(find(x>b));
		end
		while(find(x<a)>0) 
			x(find(x<a))=2*a-x(find(x<a));
		end;
		stdx =x;
		
	
