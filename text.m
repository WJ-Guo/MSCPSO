clc;clear all;
psoOptions= get_psoOptionsRosenbrock;
	 SwarmSize = psoOptions.Vars.SwarmSize;
	 Dim = psoOptions.Vars.Dim;
	 initrandSWARM = rand(SwarmSize,Dim);
	 initrandDim = rand(SwarmSize,Dim);
	   format  shortE; 
sumfxmin1 = [];sumhistory1=[];sumfxmin2 = [];sumhistory2=[];
% sumfxmin3 = [];sumhistory3=[];sumfxmin4 = [];sumhistory4=[];
% sumfxmin5 = [];sumhistory5=[];sumfxmin6 = [];sumhistory6=[];
% sumfxmin7 = [];sumhistory7=[];
MaxCrossValid =20;
for i=1:MaxCrossValid


		[fxmin1, xmin1, Swarm1, history1] = MSCPSO(get_psoOptionsRosenbrock);
		% [fxmin2, xmin2, Swarm2, history2] = DEPSO(get_psoOptionsRosenbrock);
		% [fxmin3, xmin3, Swarm3, history3] = R2HMOPSO(get_psoOptionsRosenbrock);
		% [fxmin4, xmin4, Swarm4, history4] = GPAMPSO(get_psoOptionsRosenbrock);
		% [fxmin5, xmin5, Swarm5, history5] = AESPSO(get_psoOptionsRosenbrock);
		% [fxmin6, xmin6, Swarm6, history6] = ALPSO(get_psoOptionsRosenbrock);
		% [fxmin7, xmin7, Swarm7, history7] = MSCPSO(get_psoOptionsRosenbrock);
	

	  
	
		sumfxmin1 = [sumfxmin1;fxmin1];sumhistory1=[sumhistory1,history1];
		% sumfxmin2 = [sumfxmin2;fxmin2];sumhistory2=[sumhistory2,history2];
		% sumfxmin3 = [sumfxmin3;fxmin3];sumhistory3=[sumhistory3,history3];sumfxmin4 = [sumfxmin4;fxmin4];sumhistory4=[sumhistory4,history4];
		% sumfxmin5 = [sumfxmin5;fxmin5];sumhistory5=[sumhistory5,history5];sumfxmin6 = [sumfxmin6;fxmin6];sumhistory6=[sumhistory6,history6];
		% sumfxmin7 = [sumfxmin7;fxmin7];sumhistory7=[sumhistory7,history7];
		
end;

    max(sumfxmin1)
    min(sumfxmin1)
	mean(sumfxmin1)	
	std(sumfxmin1)
	sumhistory1 = mean(sumhistory1,2);
	
    
	% max(sumfxmin2)
    % min(sumfxmin2)
	% mean(sumfxmin2)	
	% std(sumfxmin2)
	% sumhistory2 = mean(sumhistory2,2);
	
	
	% max(sumfxmin3)
	% min(sumfxmin3)
	% mean(sumfxmin3)	
	% std(sumfxmin3)
	% sumhistory3 = mean(sumhistory3,2);
	
	% max(sumfxmin4)
	% min(sumfxmin4)
	% mean(sumfxmin4)	
	% std(sumfxmin4)
    % sumhistory4 = mean(sumhistory4,2);	
	
	% max(sumfxmin5)
	% min(sumfxmin5)	
	% mean(sumfxmin5)	
	% std(sumfxmin5)
	% sumhistory5 = mean(sumhistory5,2);		
	

	% max(sumfxmin6)	
	% min(sumfxmin6)	
	% mean(sumfxmin6)	
	% std(sumfxmin6)
    % sumhistory6 = mean(sumhistory6,2);
	
    
	% max(sumfxmin7)
	% min(sumfxmin7)
	% mean(sumfxmin7)	
	% std(sumfxmin7)
    % sumhistory7 = mean(sumhistory7,2);	
	
		
		
		
		save sumhistory1MPBPSORosenbrock100.dat  sumhistory1 -ascii;
		% save sumhistory2DEPSORosenbrock100.dat  sumhistory2 -ascii;
		% save sumhistory3R2HMOPSORosenbrock100.dat  sumhistory3 -ascii;
		% save sumhistory4GPAMPSORosenbrock100.dat  sumhistory4 -ascii;
		% save sumhistory5AESPSORosenbrock100.dat  sumhistory5 -ascii;	
		% save sumhistory6ALPSORosenbrock100.dat   sumhistory6 -ascii;
		% save sumhistory7MSCPSORosenbrock100.dat sumhistory7 -ascii;
		
		load sumhistory1MPBPSORosenbrock100.dat;
		sumhistory1=sumhistory1MPBPSORosenbrock100;
		
		% load sumhistory2DEPSORosenbrock100.dat;
		% sumhistory2=sumhistory2DEPSORosenbrock100;
		
		% load sumhistory3R2HMOPSORosenbrock100.dat;
		% sumhistory3=sumhistory3R2HMOPSORosenbrock100;
		
		% load sumhistory4GPAMPSORosenbrock100.dat;
		 % sumhistory4=sumhistory4GPAMPSORosenbrock100;
		
		% load sumhistory5AESPSORosenbrock100.dat;
		% sumhistory5=sumhistory5AESPSORosenbrock100;
		
		% load sumhistory6ALPSORosenbrock100.dat;
		% sumhistory6=sumhistory6ALPSORosenbrock100;
		
		% load sumhistory7MSCPSORosenbrock100.dat;
		% sumhistory7=sumhistory7MSCPSORosenbrock100;
		
		 slen = 100;
		start =1;
		interval = 5;
		figure
		plot(start:interval:slen,(sumhistory1(start:interval:slen)),'bv--');hold on;
		% plot(start:interval:slen,(sumhistory2(start:interval:slen)),'bs-.');
		% plot(start:interval:slen,(sumhistory3(start:interval:slen)),'b*-');
		% plot(start:interval:slen,(sumhistory4(start:interval:slen)),'bo--'); 
		% plot(start:interval:slen,(sumhistory5(start:interval:slen)),'b+:'); 
		% plot(start:interval:slen,(sumhistory6(start:interval:slen)),'bx-.');
		% plot(start:interval:slen,(sumhistory7(start:interval:slen)),'b.-'); 
        %axis([start,slen,0,max([log(sumhistory1(start)),log(sumhistory2(start)),log(sumhistory3(start)),log(sumhistory4(start)),log(sumhistory5(start)),log(sumhistory6(start)),log(sumhistory7(start))])]);
        xlabel('Iteration');ylabel('Fitness value');title('Rosenbrock')
		
		legend('MPBPSO');