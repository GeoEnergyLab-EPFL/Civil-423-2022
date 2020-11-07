 %% Exercise 2.1: Sheet Pile Wall problem - Confined fluid flow
 
 % Solution for uniform permeability case. 
 
% The goal of this exercise is to understand how to solve a problem of 
% confined fluid flow by the finite element method
 
%% Mesh 

% parameters that define the domain
H=.2;T=.4;S=.5;B=.5;D=1; % see figure 1, exercise #2. We are solving 
% the problem for any of these values, you will have to analyze two 
% particular cases by the end of the exercise.

% coordinates of the different vertices that define the domain boundary.
% The origin is in the bottom-left corner.
x_exc=S+B; % middle point of excavation, x coordinate (symmetrical problem)
y_exc=D-H; % excavation height, y coordinate
x_wall=S; % horizontal position of sheet pile wall, x coordinate
y_wall=D-(H+T); % height where sheet pile wall bottom is located, y coordinate
y_ground=D; % ground level, y coordinate
t_wall=0.001; % thickness of the wall (take a small value) 

% defining outer boundary of the domain
node_coor = [ 0 0 ; x_exc 0 ; x_exc y_exc ; x_wall+t_wall/2 y_exc ;...
    x_wall+t_wall/2 y_wall ; x_wall-t_wall/2 y_wall ; x_wall-t_wall/2 ...
    y_ground; 0 y_ground ];

edge = [];

for e=1:7
edge = [ edge ; e e+1 ];
end
edge = [ edge ; 8  1 ];

% call mesh-generator MESH2D - 
 opts.kind = 'delfront';
% opts.rho2 = +1.0 ;
% opts.siz1 = 1.33;
% opts.siz2 = 1.3;
 
h_x = 0.05; %% Max elt area -> control the refinement of the mesh here
[mesh.nodes,mesh.edge, mesh.connectivity,mesh.id] = refine2(node_coor,edge,...
    [],opts,h_x) ; 

% plotting mesh

figure(1);
    patch('faces',mesh.connectivity(:,1:3),'vertices',mesh.nodes, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge,'vertices',node_coor, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',2.0) ;
 
%% Boundary conditions

% Here you have to impose the boundary condition of the problem

% First, you have to find the nodes that are located along each boundary
% For example, the nodes related to the sheet pile wall are find in the
% following way...

sheet_pile=find((mesh.nodes(:,1)==x_wall-t_wall/2 & mesh.nodes(:,2)>=y_wall)...
    | (mesh.nodes(:,1)>=x_wall-t_wall/2 & mesh.nodes(:,1)<=x_wall+t_wall/2 ...
& mesh.nodes(:,2)==y_wall) | (mesh.nodes(:,1)==x_wall+t_wall/2 & mesh.nodes(:,2)>=y_wall));

% Now you have to complete the code by finding the nodes related to the
% remaining boundaries

left_edge=------
right_edge=------
bottom_edge=------
top_ground=------
top_excavation=------

% plotting boundary nodes

figure(2)
plot(mesh.nodes(sheet_pile,1),mesh.nodes(sheet_pile,2),'or');
hold on
plot(mesh.nodes(left_edge,1),mesh.nodes(left_edge,2),'ob');
hold on
plot(mesh.nodes(right_edge,1),mesh.nodes(right_edge,2),'ob');
hold on
plot(mesh.nodes(bottom_edge,1),mesh.nodes(bottom_edge,2),'ob');
hold on
plot(mesh.nodes(top_ground,1),mesh.nodes(top_ground,2),'ob');
hold on
plot(mesh.nodes(top_excavation,1),mesh.nodes(top_excavation,2),'ob');

% Second, you have to impose the boundary conditions

% In the case of left, bottom and right edges, and sheet pile wall as well,
% the no flow boundary condition is automatically satisfy (natural 
% boundary condition) (note that finding the nodes of these boundaries was 
% useless in the end). 

% <<<Reminder: h = p/gamma_w + y>>>

% You have to define now the boundary condition at the ground level 

h_top_ground=-------

% ... and now the boundary condition at the excavation bottom

h_top_excavation=--------

% Deleting duplication of fixed nodes. 
[nodes_fixed, ia, ic]=unique([top_ground;top_excavation]) ; 
h_aux=[h_top_ground;h_top_excavation];
h_fixed=h_aux(ia); 
% Note: Sometimes some boundaries where we fie values share some nodes 
% (usually at the corners of the domain), so we have to delete the 
% duplication. This is not our case, since top_ground and top_excavation 
% do not share actually any node.

%% Solution for piezometric head, h [m]

% computing conductivity matrix

K=[1]; % Hydraulic conductivity of each material, [L/T]
% In our case is an scalar, since we are using only one material and we
% assume that it is isotropic and homogeneous
[C] = AssembleConductivityMatrix(mesh,K,'2D'); % Note: the assembly of 
% the conductivity matrix follows the same procedure that you solved the 
% past week to assemble the mass matrix (projection of flux). Have a look
% at the routine nonetheless.

% Now you have to complete the code and solve the linear system of
% equations to obtain the piezometric head at all the nodes
% (See your course notes!!)

% First, we have to separate the nodes where the piezometric head is 
% unknown from nodes where the piezometric head is fixed (because of the 
% boundary conditions).

nodes_unknows=setdiff(1:length(mesh.nodes),nodes_fixed)';

% Now, you have to complete the code and solve the final system 

% Compute the force vector
f = ---------

% Compute the unknown piezometric heads
h_unknows = -----------

% And finally we build the solution for all the nodes by gluing the solved
% h and the dirichlet boundary conditions
h = zeros(length(mesh.nodes(:,1)),1);
h(nodes_unknows)=h_unknows;
h(nodes_fixed)=h_fixed;

% plotting solution
figure(4)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),h)

%% Projection for flux -> here we use what we did last week

% estimating flux at the nodes (You finished coding up this function past
% week)

Q =ProjectFlux(mesh,'2D',K,h); %[L/s]

% plotting flux_x

figure(5)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1))

% plotting flux_y
% You can see the value for the exit gradient in this graph

figure(6)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,2))

% plotting flux as vector

figure(7)
quiver(mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1),Q(:,2),'AutoScaleFactor',3)