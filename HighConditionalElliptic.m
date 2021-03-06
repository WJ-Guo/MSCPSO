%The HighConditionalElliptic function for use with the psotoolbox
%
% Function Description:
% Equation -> sum ( 1000000.^(i-1/D-i)*x(i)^2 )
%      xmin  = [0, 0, 0.....0]  (all zeroes)
%      fxmin = 0                  (zero)
%      -100 <= x(i) <= 100
function HighConditionalEllipticed = HighConditionalElliptic(Swarm,M1,M2,shifto,lambda10,lambda100)
[SwarmSize, Dim] = size(Swarm);
Swarm1 = Swarm(:, 1:(Dim-1));
         % shifto1=shifto(1:SwarmSize,:);
		 % Swarm=Swarm-shifto1;
HighConditionalEllipticed =zeros(SwarmSize,1);
for i=1:Dim
   HighConditionalEllipticed = HighConditionalEllipticed+1000000.^((i-1)/(Dim-1))*Swarm(:,i).^2;
end   

