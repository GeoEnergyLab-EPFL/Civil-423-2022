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

xcir_r = ----------; % x coordinate of each node
ycir_r = ----------; % y coordinate of each node
node_r = [xcir_r,ycir_r]; % matrix with (x,y) coordinates of nodes     

% edges of each polygon face
edge_r = zeros(Nr,2);
edge_r(:,1) = ------; 
edge_r(:,2) = ------; 

% Inner circle - Wellbore boundary 
% Here you have to define the inner boundary of the mesh 

rw = +1.; % wellbore radius
Nw = 25; % number of faces of the polygon approximating the circle
arc_w = 2.*pi/Nw; % angle related to each polygon face

xcir_w = ---------; % x coordinate of each node
ycir_w = ---------; % y coordinate of each node
node_w = [xcir_w,ycir_w]; % matrix with (x,y) coordinates of nodes         

% edges of each polygon face
edge_w = zeros(Nw,2);
edge_w(:,1) = -----; 
edge_w(:,2) = -----; 

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

% First, we find the nodes along the perimeter of the wellbore
well_edge = find(abs(mesh.nodes(:,1).^2+mesh.nodes(:,2).^2 - rw^2) < .0001);

% Use this figure to check that your nodes are correctly founded
figure (2);
plot(mesh.nodes(well_edge,1),mesh.nodes(well_edge,2),'ob');

% Now, you have to compute the force vector
Q = 1; % wellbore injection rate [m3/s]
q = Q/(2*pi*rw); % specific discharge [m/s] (per unit thickness)

n_nodes = length(mesh.nodes); % number of nodes

% Now you have to compute the force vector
f_q = zeros(n_nodes,1); % initialization of force vector

f_q(well_edge) = --------;
%% Finite element solution + Theta-method

% Conductivity matrix
perm = 1; % Permeability coefficient 
[C] = AssembleConductivityMatrix(mesh,perm,'2D');

% Mass matrix
S = 1; % Specific storage
[M] = AssembleMassMatrix(mesh,S,'2D');

% We solve the system of equations at each time step

theta=.8; % theta parameter - time integration scheme choice (theta in [0,1])
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
    % Here you have to solve the system of equations at each step
    tn = -----; % Increase the time by delta t
    dp = -----; % Compute the increment of pressure dp
    po = -----; % Adding to the previous pressure
    pressure = -----; % Store the results in the corresponding pressure matrix
    time = -----; % Store the time in the corresponding vector
end

%% Graphs

% Use this graph to see spatio-temporal evolution of the pressure
figure(3)
times2plot = 1:round(length(time)/10):length(time);
for i = 1 : length(times2plot)
  trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),pressure(times2plot(i),:));
  title( sprintf('t = %.5g', time(times2plot(i))));
  hold all
  pause( 0.5 );
  zlim([0 max(pressure(end,:))]);
end
hold off

%% Comparison with analytical solution

% We define the points to plot
points2plot = [0 rw; 0 R/2; 0 R]; % points to plot --> wellbore, middle of
% domain, and reservoir boundary

% We find the nearest node to the points2plot
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

line_source = --------; % early-time line-source solution
cylindrical_source = --------; % steady-state cylindrical-source solution

% Create here the 3 graphs asked in the exercise and plot at each graph 
% the numerical and both analytical solutions

% First graph, at wellbore
figure(5)

% Second graph, at the reservoir boundary
figure(6)

% Third graph, at the middle of wellbore and reservoir boundary
figure(7)
