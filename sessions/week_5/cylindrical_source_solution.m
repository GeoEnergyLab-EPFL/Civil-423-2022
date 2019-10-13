function pressure = cylindrical_source_solution(xy,time,Q,k,S,rw,R)
% ANALYTICAL SOLUTION FOR RADIAL FLOW OF CYLINDRICAL SOURCE AT CONSTANT 
% INJECTION RATE IN CYLINDRICAL AND FINITE MEDIUM WITH NO-FLOW OUTER 
% BOUNDARY CONDITION

% OUTPUT
% pressure: pressure(x,t) --> Matrix (time, spatial coordinate)

% INPUT
% xy: (x,y) spatial coordinates of a set of points in the domain (size = number of points x 2)
% time: time vector (size = number of time steps)
% Q: Injection flow rate
% k: Permeability coefficient
% S: Specific storage
% rw: Wellbore radius
% R: Reservoir radius

Rd = R/rw; % dimensionless radius of reservoir
rd = sqrt((xy(:,1).^2+xy(:,2).^2))/rw; % dimensionless radius
pressure=zeros(length(time),length(xy));
for i = 1:length(time)
    td = time(i)*k/S/rw^2; % dimensionless time
    pressure(i,:) = ((Q/(2*pi*k)))*((2/(Rd^2-1))*(rd.^2/4+td)-...
        Rd^2*log(rd)/(Rd^2-1)-...
        (3*Rd^4-4*Rd^4*log(Rd)-2*Rd^2-1)/4/(Rd^2-1)^2);
end