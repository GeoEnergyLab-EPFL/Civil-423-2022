%% Exercise 5: Transient flow - Fluid injection from a wellbore into a 
% cylindrical reservoir
 
% The goal of this exercise is to understand how to solve a two dimensional 
% problem of transient flow by the finite element method using the 
% theta-method for the time integration scheme.

% The problem to solve is the pore pressure diffusion due to fluid 
% injection at constant injection rate

%% MESH

% Outer circle - Reservoir boundary 
% Here you have to define the outer boundary of the mesh 

R = +30; % reservoir radius
Nr = 25; % number of faces of the polygon approximating the circle
arc_r = 2.*pi/Nr; % angle related to each polygon face

% <<<DELETE
xcir_r = R * cos(0:arc_r:2*pi - arc_r)'; % x coordinate of each node
ycir_r = R * sin(0:arc_r:2*pi - arc_r)'; % y coordinate of each node
node_r = [xcir_r,ycir_r]; % matrix with (x,y) coordinates of nodes     

% edges of each polygon face
edge_r = zeros(Nr,2);
edge_r(:,1) = [1:1:Nr]'; 
edge_r(:,2) = [2:1:Nr,1]'; 
% DELETE>>>>>>

% Inner circle - Wellbore boundary 
% Here you have to define the inner boundary of the mesh 

rw = +1.; % wellbore radius
Nw = 25; % number of faces of the polygon approximating the circle
arc_w = 2.*pi/Nw; % angle related to each polygon face

% <<<DELETE
xcir_w = rw * cos(0:arc_w:2*pi - arc_w)'; % x coordinate of each node
ycir_w = rw * sin(0:arc_w:2*pi - arc_w)'; % y coordinate of each node
node_w = [xcir_w,ycir_w]; % matrix with (x,y) coordinates of nodes         

% edges of each polygon face
edge_w = zeros(Nw,2);
edge_w(:,1) = [1:1:Nw]'; 
edge_w(:,2) = [2:1:Nw,1]'; 
edge_w = edge_w+size(node_r,1);
% DELETE>>>>>>

% Finally, we put together all the boundaries of the domain
node = [node_r; node_w];
edge = [edge_r; edge_w];

% call mesh-generator MESH2D - download it at https://ch.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation
opts.kind = 'delfront';
% opts.rho2 = +1.0 ;
% opts.siz1 = 1.33;
% opts.siz2 = 1.3;
 
h_x = 2; %% Max elt area -> control the refinement of the mesh here
[mesh.nodes,mesh.edge,mesh.connectivity,mesh.id] = refine2(node,edge,[],opts,h_x) ;
% This time we use smooth2 function to "smooth" the mesh
[mesh.nodes,mesh.edge,mesh.connectivity,mesh.id] = smooth2(mesh.nodes,mesh.edge,mesh.connectivity,mesh.id);

% Graph to plot the mesh
figure (1);
    patch('faces',mesh.connectivity(:,1:3),'vertices',mesh.nodes, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;

%% Boundary conditions

% Build the force vector related to the boundary condition of constant 
% flux at the wellbore
% HINT: you do not need to perform integration since the specific discharge
% is contant along the perimeter of the wellbore.

% First, find the nodes along the perimeter of the wellbore
well_edge = find(abs(mesh.nodes(:,1).^2+mesh.nodes(:,2).^2 - rw^2) < .0001);

% Use this figure to check that your nodes are correctly founded
figure (2);
plot(mesh.nodes(well_edge,1),mesh.nodes(well_edge,2),'ob');

% Now, you have to compute the force vector
Q = 1; % wellbore injection rate [m3/s]
q = Q/(2*pi*rw); % specific discharge [m/s] (per unit thickness)

n_nodes = length(mesh.nodes); % number of nodes
f_q = zeros(n_nodes,1); % initialization of force vector

% <<<DELETE
f_q(well_edge) = q*(2*pi*rw)/Nw; % same amount of lumped flux at each node
% DELETE>>>
%% Finite element solution + Theta-method

% Conductivity matrix
perm = 1; % Permeability coefficient 
[C] = AssembleConductivityMatrix(mesh,perm,'2D');

% Mass matrix
S = 1; % Specific storage
[M] = AssembleMassMatrix(mesh,S,'2D');

% Solving the system of equations at each time step

theta=.5; % theta parameter - time integration scheme choice (theta in [0,1])
tMax=500; % maximum time up to which we seek the solution
iter_max=2000; % maximum number of iterations 
time_step=tMax/iter_max; % time step

po = zeros(n_nodes,1);  % initial pore pressure
pressure = [ ]; % Initialization of pore pressure solution matrix (with the 
% solution of each time step (1 step = 1 row)
pressure(1,:) = po; % Initial condition

tn=0.; % Initial time = 0
time=[tn]; % Vector where all the time steps will be stored

Id =speye(n_nodes,n_nodes); % Identity matrix necessary to solve by theta-method

j=0;
while tn<tMax && j<=iter_max
    j=j+1;
    % <<<DELETE
    tn=tn+time_step; % Increase the time by delta t
    dp=(M+theta*time_step*C)\(f_q*time_step-C*po*time_step); % Compute the increment of pressure dp
    po=po+dp; % Adding to the previous pressure
    pressure(j+1,:)=po'; % Store the results in the corresponding pressure matrix
    time(j+1)=tn; % Store the time in the corresponding vector
    % DELETE>>>
end

%% Graphs

% Use this graph to see spatio-temporal evolution of the pressure
figure(3)
times2plot = 1:round(length(time)/10/2):length(time)/2;
for i = 1 : length(times2plot)
  trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),pressure(times2plot(i),:));
  title( sprintf('t = %.5g', time(times2plot(i))));
  hold all
  pause( 0.5 );
  zlim([0 max(pressure(end,:))]);
end
hold off

%% Comparison with analytical solution

% Define points to plot
points2plot = [0 rw; 0 R/2; 0 R]; % points to plot --> wellbore, middle of
% domain, and reservoir boundary

% Find the nearest node to the points2plot
closest_node = zeros(length(points2plot),1);
for i = 1:length(points2plot)
    distances = sqrt(sum((points2plot(i,:)-mesh.nodes)'.^2));
    closest_node(i) = find(distances==min(distances));
end

% Plot to check that points were correctly found
figure(4)
patch('faces',mesh.connectivity(:,1:3),'vertices',mesh.nodes, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
hold on;
patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
hold on;
plot(mesh.nodes(closest_node,1),mesh.nodes(closest_node,2),'ob');
axis image off;

% Now you have to compare the numerical and analytical solutions

% Compute here the analytical solutions
% <<<DELETE
line_source = line_source_solution(mesh.nodes,time,Q,perm,S); % early-time line-source solution
cylindrical_source = cylindrical_source_solution(mesh.nodes,time,Q,perm,S,rw,R); % steady-state cylindrical-source solution
% DELETE>>>

% Create here the 3 graphs requested and plot at each graph the numerical 
% and both analytical solutions

% <<<<DELETE
% First graph, at wellbore
nn = 1;
figure(5)
title('Pressure versus time at the wellbore'); hold on;
plot(time,pressure(:,closest_node(nn)),'.-k') ; hold on;
plot(time,line_source(:,closest_node(nn)),'-r');hold on;
plot(time,cylindrical_source(:,closest_node(nn)),'-b');
xlabel(' time');
ylabel('Pressure at wellbore');
legend('Numerical','Early-time','Pseudo Steady-state','Location','southeast');
xlim([0 tMax])

% Second graph, at the reservoir boundary
nn = 3;
figure(6)
title('Pressure versus time at the reservoir boundary'); hold on;
plot(time,pressure(:,closest_node(nn)),'.-k') ; hold on;
plot(time,line_source(:,closest_node(nn)),'-r');hold on;
plot(time,cylindrical_source(:,closest_node(nn)),'-b');
xlabel(' time');
ylabel('Pressure at the reservoir boundary');
legend('Numerical','Early-time','Pseudo Steady-state','Location','southeast');
xlim([0 tMax])

% Third graph, at the middle of wellbore and reservoir boundary 
nn = 2;
figure(7)
title('Pressure versus time at the middle'); hold on;
plot(time,pressure(:,closest_node(nn)),'.-k') ; hold on;
plot(time,line_source(:,closest_node(nn)),'-r');hold on;
plot(time,cylindrical_source(:,closest_node(nn)),'-b');
xlabel(' time');
ylabel('Pressure at the middle');
legend('Numerical','Early-time','Pseudo Steady-state','Location','southeast');
xlim([0 tMax])
% DELETE>>>>>