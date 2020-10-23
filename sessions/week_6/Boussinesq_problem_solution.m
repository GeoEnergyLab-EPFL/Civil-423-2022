%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Boussinesq Problem - Exo week 6               %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Mesh

%------ Boundaries of the basic polygone (evaluated rectangle)
r_f = 1.; r_total = 30.; d = r_f*30;

% Define limits of the observed area
% The origin of the coordinate system is on the top-left corner and set as
% (0,0)
node = [0. 0.; r_f 0.; r_total 0.; r_total -d; 0. -d;]

edge = [];

for e=1:length(node(:,1))-1
edge = [ edge ; e e+1 ];
end
edge = [ edge ; length(node(:,1))  1 ];
part = {1:length(node(:,1))};

%------ generate a mesh fan around the critical point (border of
%foundation)

% represents the basic polygon (outer boundaries)
base.node = node;
base.edge = edge;

% set properties of the mesh fan
det.center = [r_f,0.]; % center point
det.n_fan = 30;        % number of spikes
det.l_max = 0.03;      % guides the number of circles of the mesh fan
det.angle =[pi,2*pi];  % rotation angle (around the center for spikes)
det.r_fan = 1.;        % radius of the fan (total length of spike)

base = meshFan(base,det); % generates a fan returning a set of nodes and
                          % edges including internal constraints for the
                          % fan.

% call mesh-generator MESH2D
opts.kind = 'delfront';
h_x = 3; %% Max elt area -> control the refinement of the mesh here
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = refine2(base.node,...
    base.edge,part,opts,h_x);
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = smooth2(meshL.nodes,...
    meshL.edge, meshL.connectivity,meshL.id); % Refine the mesh

% change from constant strain triangle (CST) to linear strain triangle
% (LST) (for the details see thefunction)
[meshQ.nodes,meshQ.connectivity] = Tri3ToTri6(meshL);
meshQ.id = meshL.id;
meshQ.edge = meshL.edge;

%------ Decide which mesh to use
%mesh = meshL; % Comment line 64 to get a Tri3 mesh
 mesh = meshQ; % Comment line 63 to get a Tri6 mesh

%------ plotting mesh
% Blue dotes are the additional nodes of the quadratic element.
figure(1);
plotmesh(meshL.nodes,meshL.connectivity(:,1:3),[.2 .2 .2],'w')
hold on
plot(meshQ.nodes(:,1),meshQ.nodes(:,2),'b.')
plot(meshL.nodes(:,1),meshL.nodes(:,2),'r.')
 
%% Boundary conditions

%------ Traction:
%%A linear load is applied in vertical direction on the radius r_f

% find nodes on the traction boundary
tract_boundary = find(mesh.nodes(:,1) <= r_f & mesh.nodes(:,2) == 0);

%%Get the force (f_s is a force vector for all the nodes [zero except on
%%the nodes where the force is applied])

% the applied load is 100 [kN/m] in downward z-direction ... outward normal
t_s = [0;-100*10^3]; 
 
f_s=AssembleTractionsOverLine(mesh,tract_boundary,t_s,'Axis');

% check here sum(f_s) == t_s(2)*pi (no error in the projections)
sum(f_s) - t_s(2)*pi

%------ Displacement
%%The horizontal displacement on the axis of symmetry is blocked (note that
%%this is automatically verified by the configuration of the axissymmetric
%%problem [check it if you like!]). Additionally we restrict movement on
%%the bottom and the right side (Note that this is only true if we are out 
%%of the zone of influence of the load!).

% Find the corresponding boundaries
bottom = find(mesh.nodes(:,2) == -d);       % bottom of the area
right = find(mesh.nodes(:,1) == r_total);   % right side of the area

top=find(mesh.nodes(:,2) == 0);

%%Only the horizontal movement (ui_x) is restricted on the left boundary.
%%All movements are restricted on the other two boundaries. Remeber that we
%%have defined the DOF as u1_x = 1, u1_y = 2, u2_x = 3, u2_y = 4 etc. 
fixed_nodes = unique([bottom;right]);
fixed_dof = unique([2*(fixed_nodes-1)+1;2*fixed_nodes]);

%% Assembling the system
%%We neglect the weight of the soil for the moment (ass not to calculate
%%body forces on the domain). Recall that we seek to solve the solution of 
%%the system F=Ku. In a first step we will assemble the stiffness matrix 
%%(recall that we have an axissymmetric problem).

% First we define the soil parameters
soil = [20*10^6, 0.3]; % E [Pa], nu [-]
% We consider here typical values for a loose sand

% Then we assemble the stiffness matrix
[K] = AssembleStiffnessMatrix(mesh,soil,'Axis');

%% Solving the system
%%We need to define the fixed displacements (solution to the problem which
%%we can thus eliminiate from the equation to solve). Recall that fixed
%%(non-zero) displacements generate a force term on other nodes.

u_set=fixed_dof*0;
% prepare the system to solve
eq_to_solve=setdiff([1:length(mesh.nodes)*2],fixed_dof)';
% prepare RHS
f = -K(eq_to_solve,fixed_dof)*u_set+f_s(eq_to_solve);
% Note: we need to ad here any force we know (nodal forces, surfaces forces
% and body forces) to the force vectore when we solve the system.

% We solve now the system for the unknown displacements
u_aux = K(eq_to_solve,eq_to_solve)\f;

% And finally re-assemble the full solution
u_res=0*[1:length(mesh.nodes)*2]';
u_res(fixed_dof)=u_set; 
u_res(eq_to_solve)=u_aux;

%% Plotting displacements

figure(2)
trisurf(mesh.connectivity(:,1:3),mesh.nodes(:,1),mesh.nodes(:,2),u_res(1:2:end))
title(' Radial displacements: u_r ');
xlabel(' r [m] ');
ylabel(' z [m] ');
zlabel(' u_r [m] ');

%% Plot of a deformed shape

% Reshape the displacement solution (Usol_u)
udisp = [u_res(1:2:2*length(mesh.nodes(:,1))) u_res(2:2:2*length(mesh.nodes(:,1)))];

amp_factor=1e2; % To highlight the displacements
figure(3) 
plotmesh(mesh.nodes,mesh.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(mesh.nodes+udisp*amp_factor,mesh.connectivity(:,1:3),[.8 .2 .2],'none')
title(' Deformed mesh in red, original one in black. ');
xlabel(' r [m] ');
ylabel(' z [m] ');

%% Comparison with the point load solution along r for z=0
% First we calculate the analytical solution for a point load at the origin
% (top left corner) with coordinates (0,0).
[ur_top,uz_top]=BoussinesqSolution_PointForce(mesh.nodes(top,1),...
    mesh.nodes(top,2),soil(1),soil(2));

% Extract the corresponding displacements (recall the system of DOF!).
ur_res_top=u_res((2*(top-1)+1))/(-t_s(2));
uz_res_top=u_res((2*(top)))/(-t_s(2));

%------ Plotting the solution
figure(4)
plot(mesh.nodes(top,1),ur_top,'.'); hold on
plot(mesh.nodes(top,1),ur_res_top,'.r');
title(' Compare radial displacement: red = calculation, blue = approximation ');
xlabel(' z [m] ');
ylabel(' u_r [m] ');

figure(5)
plot(mesh.nodes(top,1),uz_top,'.'); hold on
plot(mesh.nodes(top,1),-uz_res_top,'.r');
title(' Compare vertical displacement: red = calculation, blue = approximation ');
xlabel(' r [m] ');
ylabel(' u_z [m] ');

%% Comparison with the point load solution along z for r=0
% This will be strictly valid only in the far field - here the mesh is too 
% small. 

% We extract the axis of symmetry
left_boundary = find(mesh.nodes(:,1) == 0 & mesh.nodes(:,2) < 0);% axis of
                                                                 % symmetry

[ur_l,uz_l]=BoussinesqSolution_PointForce(mesh.nodes(left_boundary,1),...
    mesh.nodes(left_boundary,2),soil(1),soil(2)); % Calculate the solution

% Extract the corresponding displacements (recall the system of DOF!).
ur_res_l=u_res((2*(left_boundary-1)+1))/(-t_s(2));
uz_res_l=u_res((2*(left_boundary)))/(-t_s(2));

%------ Plotting the solution
figure(6)
plot(mesh.nodes(left_boundary,2),ur_l,'.'); hold on
plot(mesh.nodes(left_boundary,2),ur_res_l,'.r');
title(' Compare horizontal displacement: red = calculation, blue = approximation ');
xlabel(' z [m] ');
ylabel(' u_r [m] ');

figure(7)
plot(mesh.nodes(left_boundary,2),uz_l,'.'); hold on
plot(mesh.nodes(left_boundary,2),-uz_res_l,'.r');
title(' Compare vertical displacement: red = calculation, blue = approximation ');
xlabel(' z [m] ');
ylabel(' u_z [m] ');

%% Comparison with circular footing - displacement at r=0,z=0
% Analytical solution under the center of the circular footing. Note that
% there exist solution for differnt points under the footing by adopting
% numerical factors.

% Solution at z and r = 0
uz_0_0 =2*(1-soil(2)^2)*t_s(2)*r_f/soil(1)

% Compare the displacements
[ uz_0_0/t_s(2)  max(-uz_res_l) ]

% Relative error directly beneath the footing
rel_uz_0=[ uz_0_0/t_s(2) - max(-uz_res_l) ]/ (uz_0_0/t_s(2))

% Vertical displacemen at the edge of the footing.
% uz (r=x_f,z=0) = 2/pi  uz_0_0
n_a=find((mesh.nodes(:,1)== r_f) & (mesh.nodes(:,2)==0));

% Compare the displacements
[ 2/pi*(uz_0_0/t_s(2))  u_res(2*n_a)/t_s(2) ]

% Relative error at the edge of the footing
rel_uz_a=[ 2/pi*(uz_0_0/t_s(2)) - u_res(2*n_a)/t_s(2) ]/ (2/pi*(uz_0_0/t_s(2)))

%% Getting the stresses
% Similarly to the calculation of the flux we will now get the stresses at
% the centroid by averaging their value over the gauss points.
S = GetStress(mesh,soil,u_res,'Axis');
%% Project the stresses
% As the further is once more only a approximation and it is more usefull
% to get the stresses directly at the nodes. We need to project the
% stresses from the gauss points to the element nodes.
Sp = ProjectStress(mesh,'Axis',soil,u_res);

% Analytical solution along the line r=0 
% sigma_z = q*(1 - 1/(1 + (R/z)^2)^3/2 
Sp_anal = t_s(2).*(1-1./(1+(r_f./mesh.nodes(left_boundary,2)).^2).^(3/2));

%% Plot the results

figure(8)
trisurf(mesh.connectivity(:,1:3),mesh.nodes(:,1),mesh.nodes(:,2),Sp(:,1))
title(' Radial stress over the domain ( Projected S_{rr} ) ');
xlabel(' r [m] ');
ylabel(' z [m] ');
zlabel(' S_{rr} [Pa] ');

figure(9)
trisurf(mesh.connectivity(:,1:3),mesh.nodes(:,1),mesh.nodes(:,2),Sp(:,2))
hold on
plot3(mesh.nodes(left_boundary,1)...
    ,mesh.nodes(left_boundary,2),Sp_anal,'*r')
plot3(mesh.nodes(left_boundary,1)...
    ,mesh.nodes(left_boundary,2),Sp(left_boundary,2),'*b')
title(' Vertical stress over the domain ( Projected S_{zz} ) ');
xlabel(' r [m] ');
ylabel(' z [m] ');
zlabel(' S_{zz} [Pa] ');

%% Comparison between vertical displacements along z with r = 0
figure(10)
plot(mesh.nodes(left_boundary,2),Sp_anal,'*r'); hold on
plot(mesh.nodes(left_boundary,2),Sp(left_boundary,2),'.b')
title(' Vertical stress along r = 0 ( Projected S_{zz} ) ');
xlabel(' r [m] ');
ylabel(' S_{zz} [Pa] ');
