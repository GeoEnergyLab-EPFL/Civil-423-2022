function pressure = line_source_solution(xy,time,Q,k,S)
% ANALYTICAL SOLUTION FOR RADIAL FLOW OF LINE SOURCE AT CONSTANT INJECTION 
% RATE IN INFINITE MEDIUM 

% OUTPUT
% pressure: pressure(x,t) --> Matrix (time, spatial coordinate)

% INPUT
% xy: (x,y) spatial coordinates of a set of points in the domain (size = number of points x 2)
% time: time vector (size = number of time steps)
% Q: Injection flow rate
% k: Permeability coefficient
% S: Specific storage

pressure=zeros(length(time),length(xy));
for i = 1:length(time)
    r = sqrt((xy(:,1).^2+xy(:,2).^2)); % radius
    pressure(i,:) = (Q/(4*pi*k))*expint(r.^2/(4*(k/S)*time(i)));   
end