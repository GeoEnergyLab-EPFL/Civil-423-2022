%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Poroelastic sphere - Exo week 7               %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all
close all
%% Mesh generation
Radius = 1.;    % Radius of the sphere
nrad = 80;      % number of segments of the circular face used
thetas=linspace(0.,pi/2.,nrad);
node = [ 0., 0. ;
    Radius*cos(thetas)' , Radius*sin(thetas)'
    ];

edge = [];
for e=1:nrad,
    edge = [ edge ; e e+1 ];
end
edge = [ edge ; nrad+1 1 ];

%-- call mesh-gen.
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
    refine2(node,edge,[],[],0.08) ;
%-- Transform into quadratic mesh
[meshQ.nodes,meshQ.connectivity] = Tri3ToTri6(meshL);
meshQ.edge = meshL.edge;
meshQ.id = meshL.id;

% Decide which mesh to use for elasticity (meshE) and pressure (meshP)
meshE = meshL;
meshP = meshL;

% Plot the mesh and higlight the boundary.
figure(1);
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
%% Material properties
%%All stiffnesses are given in [MPa] and correspond to a ohio Sandstone.

k=8.4e3;            % elastic drained bulk modulus [MPa]
g=6.8e3;            % shear modulus [MPa]
b=0.707692;         % Biot coefficient 
M=9.18478e3;        % Biot Modulus
k_u=k+b^2*M         % undrained bulk modulus [MPa]
perm =0.137549e-3;  % permeability value adjusted such that c= k /(mu S)=1 
                    % , where S is storage coefficient

mu_f=1;             % fluid viscosity [Pa s]
rho=1;              % density
kappa = perm / mu_f;% conductivity 

nu_u=(1-2*(g)/(3*k_u))/(2+2*(g)/(3*k_u));   % undrained poisson's ratio
E_u = (9*(k_u))/(1+3*(k_u)/(g));            % undrained youngs modulus           
params_u = [E_u, nu_u];

nu=(1-2*(g)/(3*k))/(2+2*(g)/(3*k));         % poisson's ratio
E=(9*(k))/(1+3*(k)/(g));                    % youngs modulus
params = [E, nu];

%% Boundary conditions
%%We will impose that there are no radial displacements on the left
%%boundary and no vertical displacements on the bottom boundary.
%%Additionally, the pressure on the outer boundary (circle) will be fixed).

%----- Displacements
% find the corresponding nodes
bottom = find(meshE.nodes(:,2)==0.);
left = find(round(meshE.nodes(:,1),4) ==0.);

% adjust the fixed dofs and nodes
fixed_dof_left = [2*(left-1)+1];
fixed_dof_bottom = [2*bottom];
fixed_nodes = unique([bottom,left]);
fixed_dof = unique([fixed_dof_left;fixed_dof_bottom]);

%----- Pressure
% nodes on radius. Note that it will be necessary to introduce a rounding
% as to realy capture all the points. This is because we need to account
% for the numerical imprecision of the meshing.
circle=find(round(sqrt(meshP.nodes(:,1).^2.+meshP.nodes(:,2).^2),4)==Radius);

%----- Plotting the boundarys onto the initial plot of the mesh
hold on;
plot(meshE.nodes(bottom,1),meshE.nodes(bottom,2),'or')
hold on;
plot(meshE.nodes(left,1),meshE.nodes(left,2),'or')
hold on;
plot(meshP.nodes(circle,1),meshP.nodes(circle,2),'*b')

%% Initial condition
%%As an initial condition we will set a unifrom pressure field withou
%%deviatoric stresses all equal to -1. Note that it would also be possible
%%to define the initial state by a traction applied on the boundary which
%%results in the same final result.

%----- Initial stress field
% Set the initial stressfield
Sig_unif = [-1;-1;0;-1];
% Assemble the stress vector at the nodes
Sig_o = SetStressField(meshE,'Axis',...
    {[1:length(meshE.nodes(:,1))],Sig_unif});
%% Assembling the different matrices

% elasticity
[K]=AssembleStiffnessMatrix(meshE,params,'Axis');

% mass matrix term
[Mass]=AssembleMassMatrix(meshP,1/M,'Axis');

% laplacian term
[C]=AssembleConductivityMatrix(meshP,kappa,'Axis');

% Coupling term 
[A]=AssembleCouplingMatrix(meshE,meshP,b,'Axis');

%% Solving directly for the undrained solution at t=0 ! 
%%We dont fix the pore pressure to zero on the outer radius for the
%%undrained response.
dt= 0.;         % time step
AA=Mass+dt*C;   % lower right part of the total matrix

%----- Getting the different number of DOF and assembling the total matrix
% getting DOF
ntot_P = length(C(:,1));
ntot_E = length(K(:,1));
ntot=ntot_E+ntot_P;

% Assemble the matrix
TotMat=sparse(ntot,ntot);
TotMat(1:ntot_E,1:ntot_E)=K;
TotMat(ntot_E+1:ntot,ntot_E+1:ntot)=-AA;
TotMat(ntot_E+1:ntot,1:ntot_E)=-A';
TotMat(1:ntot_E,ntot_E+1:ntot)=-A;
%----- Setting the force vector to the initial conditions
Ftot=sparse(ntot,1);
Ftot(1:ntot_E) = Sig_o;

%----- Preparing the system
% Get free equations
eq_free_u = setdiff([1:ntot_E],fixed_dof)';
eq_free_p=[1:ntot_P]';
eq_free= [eq_free_u;eq_free_p+ntot_E];

%----- Solving the system
Undrained_sol=sparse(ntot,1);
Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

%----- Seperate displacements from pressure
pp_o=Undrained_sol(ntot_E+1:end);
U_o=sparse(ntot_E,1);
U_o(eq_free_u) = Undrained_sol(eq_free_u);

%----- Defining the analytical solution along r   for the undrained displacement u_r
AnalyticSolUr= @(r) -(r*(1-2* nu_u)/(2.* g*(1.+ nu_u)));

%----- Comparing results
% radial displacement at z=0 and at r=0. They should be the same and
% corresponding to the Analytical solution.
figure(2)
plot(meshE.nodes(left,2), U_o([2*left]),'*k' )
hold on
plot(meshE.nodes(bottom,1), U_o([2*(bottom-1)+1]),'sb','MarkerSize',12.5)
plot(meshE.nodes(left,2),AnalyticSolUr(meshE.nodes(left,2)),'-r');

%----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [U_o(1:2:end) U_o(2:2:end)];

figure(3)
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*1e3,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title('deformed mesh')

%----- t=0+ solution is the undrained response
%pp_o(circle)=0.; % NOW we fix the pore pressure to zero. 

% plotting the pressure profile
figure(4) 
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2),full(pp_o))
xlabel(' r ');
ylabel(' z ');
zlabel(' p_o '); 
%% Solving now the drained response for t > 0

%
pp_o(circle)=0.; % NOW we fix the pore pressure to zero on the outer surface... 


dt=0.002; % setting a the constant time step

%----- matrix for fluid flow for ct time step
% this would have to be adapted whenever we have a non- constant time step.
AA=Mass+dt*C;

%----- Assembly of the matrix only requires a change to the changed part
TotMat(ntot_E+1:ntot,ntot_E+1:ntot)=-AA;

%----- Redefining the free equations (as it changes in pressure)
eq_free_p=setdiff([1:ntot_P],circle)';
eq_free= [eq_free_u;eq_free_p+ntot_E]';

%----- initialize solution vectors
t_k(1)=0.;              % time
hist_pp_1(1)=pp_o(1);   % history of pp at the sphere center
n_step=501;             % numer of time steps

%----- time loop (from two to the number of iterations chosen 
for i=2:n_step
    
    t_k(i)=t_k(i-1)+dt; % store the time
    
    Dp=0.*pp_o;         % at nodes for fluid flow
    D_U=0.*U_o;         % at nodes for elasticity
     
    D_x=[D_U; Dp];      % Assemble solution vector
    
    % calculate the change in flux
    po_lhs=dt*C(eq_free_p,:)*pp_o(:);
    
    % Set the right hand side
    P_s=sparse(ntot,1);
    P_s(eq_free_p+ntot_E)=po_lhs;
    
    % Solve the system
    D_x(eq_free)=TotMat(eq_free,eq_free)\P_s(eq_free);
    
    % Seperate elastic and flow components
    Dp(eq_free_p)=D_x(eq_free_p+(ntot_E));
    D_U(eq_free_u)=D_x(eq_free_u);
    
    % set new pressure and displacement
    pp_o=pp_o+Dp  ;
    U_o=U_o+D_U;
    
    % save pressure history at node ID 1
    hist_pp_1(i)=pp_o(1);
end

%% Visualizing the solution

%----- Load the analytical solution for the pore pressure
ResA = csvread("PP-Sphere-ohio.csv");

%----- plot the pressure evolution at the center of the sphere
figure(5)
plot(t_k,hist_pp_1,'k') ; hold on;
plot(ResA(:,1),ResA(:,2),'r') 
title(' Pore pressure evolution at the sphere center');
xlabel('time');
ylabel(' pore pressure / applied load');



%----- plot the error on the pressure calculation
% calculate the error
rel_error = abs(ResA(1:n_step,2)-hist_pp_1(1:end)')./ResA(1:n_step,2);
figure(6)
plot(t_k,rel_error);
xlabel('time');
ylabel(' relative error on pressure '); 

%----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [U_o(1:2:end) U_o(2:2:end)];

figure(7) 
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*1e3,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title([' Deformed mesh at time t=', num2str(t_k(end))]);

%----- Plotting the pore pressure
figure(8) 
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2)...
    ,full(pp_o))
title([' Pore pressure distribution at time t=', num2str(t_k(end))]);
xlabel('r');
ylabel('z');
ylabel('p');