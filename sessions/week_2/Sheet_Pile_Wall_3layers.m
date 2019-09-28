 %% Exercise 2.2: Sheet Pile Wall problem - Confined fluid flow
 
 % Solution for three soil layers
 
% The goal of this exercise is to understand how to extend the code in 
% order to consider more than one material in the domain, and to understand 
% the impact of a intermediate layer with higher permeability

%% Mesh 

% parameters that define the domain

H=.3;T=.3;S=.5;B=.5;D=1; % see figure 1, exercise #2
x_exc=S+B; % middle point of excavation, x coordinate (symmetrical problem)
y_exc=D-H; % excavation height, y coordinate
x_wall=S; % horizontal position of sheetpile wall, x coordinate
y_wall=D-(H+T); % height above sheetwall is located, y coordinate
y_ground=D; % ground level (side without excavation), y coordinate
t_wall=0.001; % thickness of the wall (take a small value) 

y_lay12=.6;
y_lay23=.2; % must be less than y_12

% defining outer boundary of the domain
node_outer = [ 0 0 ; x_exc 0 ; x_exc y_exc ; x_wall+t_wall/2 y_exc ;...
    x_wall+t_wall/2 y_wall ; x_wall-t_wall/2 y_wall ; x_wall-t_wall/2 y_ground; 0 y_ground ];
edge_outer = [];

for e=1:7
edge_outer = [ edge_outer ; e e+1 ];
end
edge_outer = [ edge_outer ; 8  1 ];

% defining internal constraints of the mesh (layers)

node_constraints = [ 0 y_lay12 ; x_wall-t_wall/2 y_lay12 ;... % boundary between layers 1 and 2, left side of wall
    x_wall+t_wall/2 y_lay12 ; x_exc y_lay12 ; ... % boundary between layers 1 and 2, right side of wall
    0 y_lay23 ; x_exc y_lay23 ]; % boundary between layers 2 and 3

edge_constraints = [ 9 10 ; 11 12 ; 13 14 ];

% outer + internal constraints (layers)

node = [node_outer;node_constraints];
edge = [edge_outer;edge_constraints];
part{1} = 1:8; % edges being part of outer boundary to call refine2 function

% call mesh-generator MESH2D - download it at https://ch.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation
 opts.kind = 'delfront';
% opts.rho2 = +1.0 ;
% opts.siz1 = 1.33;
% opts.siz2 = 1.3;
 
h_x = 0.05; %% Max elt area -> control the refinement of the mesh here
[mesh.nodes,mesh.edge, mesh.connectivity,mesh.id] = refine2(node,edge,...
    part,opts,h_x) ; 

% plotting mesh

figure(1);
    patch('faces',mesh.connectivity(:,1:3),'vertices',mesh.nodes, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge_outer,'vertices',node_outer, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',2.0) ;
    hold on;
    plot(node_constraints(1:2,1),node_constraints(1:2,2),'-',...
        'color','black','linewidth',2.0)
    hold on;
    plot(node_constraints(3:4,1),node_constraints(3:4,2),'-',...
        'color','black','linewidth',2.0)
    hold on;
    plot(node_constraints(5:6,1),node_constraints(5:6,2),'-',...
        'color','black','linewidth',2.0)
 
%% Boundary conditions

% Finding boundary nodes

left_edge=find(mesh.nodes(:,1)==0);
right_edge=find(mesh.nodes(:,1)==x_exc);
bottom_edge=find(mesh.nodes(:,2)==0);
top_ground=find(mesh.nodes(:,2)==y_ground);
top_exc=find(mesh.nodes(:,2)==y_exc & mesh.nodes(:,1)>=x_wall);
sheet_pile=find((mesh.nodes(:,1)==x_wall-t_wall/2 & mesh.nodes(:,2)>=y_wall)...
    | (mesh.nodes(:,1)>=x_wall-t_wall/2 & mesh.nodes(:,1)<=x_wall+t_wall/2 ...
& mesh.nodes(:,2)==y_wall) | (mesh.nodes(:,1)==x_wall+t_wall/2 & mesh.nodes(:,2)>=y_wall));

% plotting boundary nodes

figure(2)
plot(mesh.nodes(left_edge,1),mesh.nodes(left_edge,2),'ob');
hold on
plot(mesh.nodes(right_edge,1),mesh.nodes(right_edge,2),'ob');
hold on
plot(mesh.nodes(bottom_edge,1),mesh.nodes(bottom_edge,2),'ob');
hold on
plot(mesh.nodes(top_ground,1),mesh.nodes(top_ground,2),'ob');
hold on
plot(mesh.nodes(top_exc,1),mesh.nodes(top_exc,2),'ob');
hold on
plot(mesh.nodes(sheet_pile,1),mesh.nodes(sheet_pile,2),'or');

% Assigning boundary conditions

% left, bottom, right and sheetpile --> No flow (natural boundary condition)
% top_ground p = gamma_w h1 --> h = h1 + y_top
h1=0; % we can define any piezometric head h1
h_top_ground=h1+mesh.nodes(top_ground,2);
% top_exc p = gamma_w h2 --> h = h2 + y_exc
h2=0; % we can define any piezometric head h2
h_top_exc=h2+mesh.nodes(top_exc,2);

[nodes_fixed, ia, ic]=unique([top_ground;top_exc]) ; % deleting duplication of fixed ones (just in case)
h_aux=[h_top_ground;h_top_exc];
h_fixed=h_aux(ia);

%% Assigning different materials/permeabilities per layer

% Reminder: mesh structure = [mesh.nodes,mesh.edge, mesh.connectivity,mesh.id]

for e=1:length(mesh.connectivity(:,1)) % there should be a way to do it without loop, using find function for example
    n_e = mesh.connectivity(e,:); % label of nodes per each element
    coor = mesh.nodes(n_e,:); % coordinates of the respective nodes
    if coor(:,2)>=y_lay12;
        mesh.id(e)=1;
    elseif coor(:,2)>=y_lay23 & coor(:,2)<=y_lay12;
        mesh.id(e)=2;
    else coor(:,2)<=y_lay23;
        mesh.id(e)=3;
    end 
end
e1=find(mesh.id==1); % elements belonging to material 1
e2=find(mesh.id==2); % elements belonging to material 2
e3=find(mesh.id==3); % elements belonging to material 3

% plotting different materials (layers)

figure(3);
    patch('faces',mesh.connectivity(e1,:),'vertices',mesh.nodes, ...
        'facecolor','w','edgecolor','r') ;
    hold on; axis image off;
    patch('faces',mesh.connectivity(e2,:),'vertices',mesh.nodes, ...
        'facecolor','w','edgecolor','b') ;
    hold on;
    patch('faces',mesh.connectivity(e3,:),'vertices',mesh.nodes, ...
        'facecolor','w','edgecolor','black') ;
    hold on;
    patch('faces',edge_outer,'vertices',node_outer, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',2.0) ;

%% Solution for piezometric head, h [m]

% computing conductivity matrix

K=[1;10;1]; % Hydraulic conductivity of each layer, [L/s]
[C] = AssembleConductivityMatrix(mesh,K,'2D');

% separating unknowns from fixed nodes

unknowns=setdiff(1:length(mesh.nodes),nodes_fixed)';

% matrix solution

f = -C(unknowns,nodes_fixed)*h_fixed;
h_unknows = C(unknowns,unknowns)\f;
h = zeros(length(mesh.nodes(:,1)),1);
h(unknowns)=h_unknows;
h(nodes_fixed)=h_fixed;

% plotting solution
figure(4)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),h)

%% Projection for flux

% estimating flux at the nodes

Q =ProjectFlux(mesh,'2D',K,h); %[L/s]

% plotting flux_x

figure(5)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1))

% plotting flux_y

figure(6)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,2))

% plotting flux as vector

figure(7)
quiver(mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1),Q(:,2),'Color','black','AutoScaleFactor',4)

% Exit gradient

exitnode = find(mesh.nodes(:,1)==x_wall+t_wall/2 & mesh.nodes(:,2)==y_exc);
I_E=Q(exitnode,2);