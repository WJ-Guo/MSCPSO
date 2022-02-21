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

% Initialize Swarm and Velocity��ʼ����PSO����������
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

f2eval = psoOptions.Obj.f2eval; %The objective function to optimize.���Ż���Ŀ�꺯��

%Find initial function values.
fSwarm = feval(f2eval, Swarm,M1,M2,shifto,lambda10,lambda100);%��Сֵ���
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


%��ʼ��ֵ 0.0005(0,1)  0.005(-5.12,5.12),0.05(-50,50),0.1(-100,100),0.3(-300,300)
T_d = repmat(0.1,1,psoOptions.Vars.Dim); %�γ�һ��������[0.1 0.1 0.1 0.1]��ʾΪV����ֵ
F_d = repmat(0,1,psoOptions.Vars.Dim);   %%
K_1 = 5;
K_2 = 10;

Initmutestd = (psoOptions.Obj.ub-psoOptions.Obj.lb);%��ʼ������

M = 2;                %%%scale
P = SwarmSize/M;
	%��ʼ����׼����
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
    
    			
%CMEPSO������һ�ֶ�߶ȣ����򡢶�ʱ�ı�������������ٶ�С��һ�����������
    %�Զ��ּҲ�Ѱ���¼��۵ط��ܺ��������������ͨ�������˶�Ѱ���µ�����ء���ʱ����
    %�������ʷ���λ�ö�����ס��Ⱥ�����λ�á�			
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
  
    
    
		   	
    %Ϊ�˷�ֹ��������������������������Ӷ����½����ص㣬ESPSO�����������������Ա��ӳ��ֲ���С��
    EsacepeRows = (abs(VStep) <= repmat(T_d,SwarmSize,1));                   %%%particles satisfy mutation condition      
     %EsacepeRows = (VStep == 0);
    [rowE,col,val]=find(EsacepeRows);                           %%%save position
		
		rowOE=unique(rowE);
		
		for i=1:length(rowOE)
				 curSwarm = Swarm(rowOE(i),:);              %%��ȡ��Ⱥ�ɱ��������
				 
				
			
				 colE = col(find(rowE==rowOE(i)));          %�����������г�������С���Լ�ά����ֵ����
				 
				 VStepBak = VStep(rowOE(i),:);              
			
				 VStepBak1 = VStepBak;
				 if (M~=1)
				 		VStepBak1remp = repmat(VStepBak1,M+2,1);
				 else
				  	VStepBak1remp = repmat(VStepBak1,M+1,1);
				 end
				 
		     VStepBak1remp(1:M,colE) = randn(M,length(colE)).*repmat(InitmutestdQ',1,length(colE));
		     %InitmutestdQÿһ��ά�ȵķ����ͬ
		     
   			 VStepBak1remp(M+1,colE) = rand(1,length(colE)).*psoOptions.SParams.Vmax;	
   			 %�����ֱ�����ԣ����һ����Vmax�����ȱ���
   			 
   			 
   			 %ʵ�ָ��ֱ��첢������Ӧ�Ľ��
   			 if (M~=1)
   			 		VStepBak1remp(M+2,colE) = randn(1,length(colE)).*repmat(Initmutestd,1,length(colE));
		        SwarmBak = repmat(curSwarm,M+2,1) + psoOptions.SParams.Chi .* VStepBak1remp;    % Evaluate new Swarm
         else
           SwarmBak = repmat(curSwarm,M+1,1) + psoOptions.SParams.Chi .* VStepBak1remp;   % Evaluate new Swarm
         end	
   
         %�������Ŀ��ܳ��ֵĳ������
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
         		fevals = fevals + M+1;%��¼����Ŀ�꺯���Ĵ���
         end
         
  	  	[sortfSwarm,indxfSwarm]=sort(fSwarmBak);
  	  	%���ı����Ľ������
  	  	Swarm(rowOE(i),:) = SwarmBak(indxfSwarm(1));
  	 end 	
  	
  	
    
    fSwarm = feval(f2eval, Swarm,M1,M2,shifto,lambda10,lambda100);
    fevals = fevals + SwarmSize;
    [sortfSwarm,findex] = sort(fSwarm);	          %%��������
    
    	
    	%���㷽��ı仯��С,����Ⱥ����ΪM�ݣ�ÿ�ݺ�SWarmSize/M������
				for i =1:M,
					FixXi(i)=sum(sortfSwarm((i-1)*P+1:i*P));	
				end;
				
				FixXimax = max(FixXi);
				FixXimin = min(FixXi);
				
				if (FixXimax~=FixXimin)
				    InitmutestdQ = InitmutestdQ.*exp((M*FixXi-sum(FixXi))/(M*(FixXimax-FixXimin)) );
			 	end
			 	
				%��ֹ�任��ķ������	
				InitmutestdQ=StdQstandard(InitmutestdQ,psoOptions.Obj.ub-psoOptions.Obj.lb);
			 
        
        
        
        B_id = zeros(SwarmSize,psoOptions.Vars.Dim);
		    B_id(find(EsacepeRows))=1;
		    %����Щv���õ�����Ӧ��swarm�����Ӧ���б�־Ϊ1
		    
		    
		    b_SUM = sum(B_id,1);%�γ�������
		    
		    F_d = F_d + b_SUM;%��¼�������е�ÿһ��Dim������С����ֵ�Ĵ�������¼�����ۻ���
		      
		   	T_d(find(F_d > K_1))=T_d(find(F_d > K_1))./K_2;
		   	%�����������ֵ̫����Ҫ����һЩ��
		   	F_d(find(F_d > K_1))=0; %�����Ѿ����˵������������㣬���¼����ۻ�
        
        
    
    % Updating the best position for each particle���¾ֲ�����
    changeRows = fSwarm < fPBest;
    fPBest(find(changeRows)) = fSwarm(find(changeRows));
    PBest(find(changeRows), :) = Swarm(find(changeRows), :);
    
    lastbpart = PBest(g, :);%��¼ԭ��������λ���־���һ�ֺ����ڱ���ɵ�����ֵ
    
    %AEPSO������һ�ֶ��򡢶�ʱ�ı�������������ٶ�С��һ�����������
    %�Զ��ּҲ�Ѱ���¼��۵ط��ܺ��������������ͨ�������˶�Ѱ���µ�����ء���ʱ����
    
    
   % fPBest(rowOE) = fSwarm(rowOE);
   % PBest(rowOE,:) = Swarm(rowOE,:);
    
    %�������ʷ���λ�ö�����ס��Ⱥ�����λ�á�
    
 
    
   
    
    
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
		
	