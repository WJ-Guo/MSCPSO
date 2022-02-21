function Sphered = Sphere(Swarm,M1,M2,shifto,lambda10,lambda100)
   [SwarmSize, Dim] = size(Swarm);
    % shifto1=shifto(1:SwarmSize,:);
    % Swarm=Swarm-shifto1;
   
    Sphered =  sum(Swarm.^2,2);