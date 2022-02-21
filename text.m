clc;clear all;
psoOptions= get_psoOptionsRosenbrock;
	 SwarmSize = psoOptions.Vars.SwarmSize;
	 Dim = psoOptions.Vars.Dim;
	 initrandSWARM = rand(SwarmSize,Dim);
	 initrandDim = rand(SwarmSize,Dim);
	   format  shortE; 
sumfxmin1 = [];sumhistory1=[];

MaxCrossValid =20;
for i=1:MaxCrossValid


		[fxmin1, xmin1, Swarm1, history1] = MSCPSO(get_psoOptionsRosenbrock);
		sumfxmin1 = [sumfxmin1;fxmin1];sumhistory1=[sumhistory1,history1];		
end;

    max(sumfxmin1)
    min(sumfxmin1)
	mean(sumfxmin1)	
	std(sumfxmin1)
	sumhistory1 = mean(sumhistory1,2);
	
		save sumhistory1MPBPSORosenbrock100.dat  sumhistory1 -ascii;
		
		load sumhistory1MPBPSORosenbrock100.dat;
		sumhistory1=sumhistory1MPBPSORosenbrock100;
			
		 slen = 100;
		start =1;
		interval = 5;
		figure
		plot(start:interval:slen,(sumhistory1(start:interval:slen)),'bv--');hold on;

        %axis([start,slen,0,max([log(sumhistory1(start)),log(sumhistory2(start)),log(sumhistory3(start)),log(sumhistory4(start)),log(sumhistory5(start)),log(sumhistory6(start)),log(sumhistory7(start))])]);
        xlabel('Iteration');ylabel('Fitness value');title('Rosenbrock')
		
		legend('MPBPSO');
