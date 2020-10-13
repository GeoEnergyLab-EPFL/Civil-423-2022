%% Exercise 4: Transient flow - Terzaghi's Consolidation Problem 
 
% The goal of this exercise is to understand how to solve a one dimensional 
% problem of transient flow by the finite volume / finite difference method 
% and using the theta-method for the time integration scheme.

% The problem to solve is the pore pressure diffusion of the Terzaghi's
% consolidation problem.

%% Mesh

ne =100; % Number of elements / nodes
h = 1./ne; % Size of elements
p_o = 1; % Initial pressure
c = 1; % Consolidation coefficient
l = 1; % Size of the domain (thickness of soil layer)
xend =linspace(0,1.,ne+1); % Coordinates at the ends of each element
%(note that "ghost cell" is not considered here, it is "located" at -h/2)
xs=xend(1:ne)+h/2; % Coordinates at the center of each element

%% L Matrix

L=sparse(ne,ne); % Initialization of L Matrix

% Now you have to fill the L Matrix

%% Theta-method and numerical solution

theta=.8; % theta parameter - time integration scheme choice (theta in [0,1])
time_step=30*h^2./2.; % Size of the time step as times * CFL condition   
tMax=1.; % Maximum time up to which we seek the solution
iter_max=2000; % Maximum number of iterations 

po =p_o*ones(ne,1);  % Initial pore pressure
pressure=[ ]; % Initialization of pore pressure solution matrix (with the 
% solution of each time step (1 step = 1 row)
pressure(1,:)=po; % Initial condition

tn=0.; % Initial time = 0
time=[tn]; % Vector where all the time steps will be stored

Analytical =[]; % Matrix for the corresponding analytical solution results 
Analytical(1,:)=po; % Initial condition   

Id =speye(ne,ne); % Identity matrix necessary to solve by theta-method

j=0;
while tn<tMax && j<=iter_max
    j=j+1;
    tn=--------; % Increase the time by delta t
    dp=--------; % Compute the increment of pressure dp
    po=--------; % Adding to the previous pressure
    pressure--------; % Store the current pressure results in the corresponding pressure matrix
    time--------; % Store the current time in the corresponding vector
    Analytical---------; % Compute the corresponding analytical
    %solution at the center of each element / cell for the given time tn
end

%% Graphics - Numerical versus analytical solution

% Plot pressure profile at different times

figure(1)
title('Pressure profile at different times'); hold on;
h1=plot(xs,pressure(1:round(length(time)/10):length(time),:),'.-k'); hold on;
h2=plot(xs,Analytical(1:round(length(time)/10):length(time),:),'-r'); hold on;
xlabel('x');
ylabel('Pressure');
legend([h1(1),h2(1)],'Numerical','Analytical');

% Plot pressure versus time at a given element / cell

figure(2)
cell_id = 20; % element / cell number to be plotted
title(strcat('Pressure versus time at  cell #',string(cell_id))); hold on;
plot(time,pressure(:,cell_id),'.-k') ; hold on;
plot(time,Analytical(:,cell_id),'-r');
xlabel(' time');
ylabel(strcat('Pressure at cell #', string(cell_id)));
legend('Numerical','Analytical');

% Plot relative error at a given element / cell
rel_error=abs(Analytical(:,10)-pressure(:,10))./Analytical(:,10);
figure(3)
plot(time,rel_error); hold on;
xlabel('time');
ylabel(strcat('relative error on Pressure at cell #', string(cell_id)));