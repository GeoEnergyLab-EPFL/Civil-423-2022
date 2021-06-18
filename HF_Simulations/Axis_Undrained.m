%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Axissymmetry poroelastic fracture             %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all
close all

%% Material properties
%%All stiffnesses are given in [MPa] and correspond to a ohio Sandstone.

k=8.4e3;            % elastic drained bulk modulus [MPa]
g=6.8e3;            % shear modulus [MPa]
b=0.707692;         % Biot coefficient 
M=9.18478e3;        % Biot Modulus
k_u=k+b^2*M;        % undrained bulk modulus [MPa]
perm = 25*0.137549e-3;  % permeability value adjusted such that c= k /(mu S)=1 
                    % , where S is storage coefficient
B =(k_u-k)/(b*k_u);
mu_f=1;             % fluid viscosity [Pa s]
rho=1;              % density
kappa = perm / mu_f;% conductivity 

nu_u=(1-2*(g)/(3*k_u))/(2+2*(g)/(3*k_u));   % undrained poisson's ratio
E_u = (9*(k_u))/(1+3*(k_u)/(g));            % undrained youngs modulus           
params_u = [E_u, nu_u];

nu=(1-2*(g)/(3*k))/(2+2*(g)/(3*k));         % poisson's ratio
E=(9*(k))/(1+3*(k)/(g));                    % youngs modulus
params = [E, nu];

c = 2 * kappa * g * (1-nu)*(nu_u-nu) / ...
    (b^2 * (1-2*nu)^2 * (1-nu_u)); % diffusivity

eta = b*(1-2*nu_u)/(2*(1-nu_u)); % poroelastic eta parameter

beta = 1e-5; % fluid compressibility
%% Calculate initial fracture variables
Vo = 0.01; % Fracture volume
KIc = 1.5; % Fracture toughness
Kp = 4*(2/pi())^(1/2)* KIc; 
Ep = params(1) / (1-params(2)^2); % Drained Plane strain modulus
Ep_u = params_u(1) / (1-params_u(2)^2); % Undrained Plane strain modulus

Radius = 2/pi()^(2/3)*(Ep * Vo / Kp)^(2/3); % Fracture radius from volume

p_f = pi()^(1/3)/8 * (Kp^4 / (Vo * Ep))^(1/3); % Fracture pressure
%% Mesh generation
l_area = floor(5*Radius);    % length of observed area
d_area = floor(5*Radius);    % Depth of observed area

nrad = 500;     % number of segments in the tip
nfr = 2500;     % number of segments in the fracture

tip_zone = flip([Radius-flip(logspace(-5,-1,nrad));
            -d_area*ones(1,nrad)],2)'; % nodes of the tip zone
external_approach = [Radius+flip(logspace(-5,-1,nrad));
            -d_area*ones(1,nrad)]'; % nodes ahead of the tip zone
fracture_zone = flip([linspace(0,Radius-0.1,nfr);
            -d_area*ones(1,nfr)],2)';   % nodes from injection point to the
                                        % tip zone

node = unique([ 0., 0. ;
    l_area, 0. ;
    l_area, -d_area;
    external_approach;
    Radius, -d_area;
    tip_zone;
    fracture_zone],'rows','stable');
node(end,1) = 0.; % combine all nows

edge = [];
for e=1:length(node)
    edge = [ edge ; e e+1 ];
end
edge(end,2) = 1; % get the edges

%-- call mesh-gen.
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
    refine2(node,edge,[],[],4) ;
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
    smooth2(meshL.nodes,meshL.edge, meshL.connectivity,meshL.id); 
    % Refine the mesh

%-- Transform into quadratic mesh
[meshQ.nodes,meshQ.connectivity] = Tri3ToTri6(meshL);
meshQ.edge = meshL.edge;
meshQ.id = meshL.id;

% Decide which mesh to use for elasticity (meshE) and pressure (meshP)
meshE = meshQ; % Quadratic elasticity
meshP = meshL; % Linear pressure

% Plot the mesh and higlight the boundary.
fig_ind = 1;
figure(fig_ind);
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
fig_ind = fig_ind + 1;

% Controlling the CFL condition of the numerical solution
hmin=(pi/2)/nrad ;
dt_cfl=hmin^2

%% Boundary conditions

%----- Boundaries
% find the corresponding nodes
bottom = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)>=Radius);
bottom_P = find(meshP.nodes(:,2)==-d_area & meshP.nodes(:,1)>=Radius);
left = find(meshE.nodes(:,1)==0. & meshE.nodes(:,2)>-d_area );
top = find(meshE.nodes(:,2)==0);
right = find(meshE.nodes(:,1)==l_area);

%----- Displacements
% adjust the fixed dofs and nodes
fixed_dof_left = [2*(left-1)+1];
fixed_dof_bottom = [2*bottom];
fixed_nodes = unique([bottom;left]);
fixed_dof = unique([fixed_dof_left;fixed_dof_bottom]);

%----- Pressure
fracture_E = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)<Radius);
fracture_P = find(meshP.nodes(:,2)==-d_area & meshP.nodes(:,1)<=Radius);

%----- Plotting the boundarys onto the initial plot of the mesh
% Boundaries where displacements are fixed will be red, such that are fixed
% in pressure blue and the ones fixed in tractions green
hold on;
plot(meshE.nodes(bottom,1),meshE.nodes(bottom,2),'or')
hold on;
plot(meshE.nodes(left,1),meshE.nodes(left,2),'or')
hold on;
plot(meshP.nodes(fracture_P,1),meshP.nodes(fracture_P,2),'*b')
hold on;
plot(meshE.nodes(right,1),meshE.nodes(right,2),'sg')
hold on;
plot(meshE.nodes(top,1),meshE.nodes(top,2),'sg')
%% Initial condition
%%As an initial condition we will set a unifrom pressure field with
%%deviatoric stresses where s_x  and s_y =  but without shear
%%stresses. Additionally we will define an initial pore pressure field
%%which is uniform as well.

%----- Initial stress field
% Set the initial stressfield in MPa
%
Po = 30; % mean compressive 
So = 0; % deviatoric stress
Sig_unif = -[(Po-So);(Po+So);0;0];  %  [MPa]
% Assemble the stress vector at the nodes
Sig_o = SetStressField(meshE,'Axis',...
    {[1:length(meshE.nodes(:,1))],Sig_unif});

%----- Initial pore pressure field
% Define the fixed and free nodes
eq_free_p = setdiff([1:length(meshP.nodes(:,1))],fracture_P)';
eq_fixed_p = setdiff([1:length(meshP.nodes(:,1))],eq_free_p)';

% Define the  initial pore pressure field
p_o = 5; %[MPa]
p_f = Po + p_f;
p_fixed = sparse(length(meshP.nodes(:,1)),1);
p_fixed(eq_free_p) = p_o;
p_fixed(eq_fixed_p) = p_f;
%% Set boundary tractions
%%We will need to set the corresponding surface tractions on the "free"
%%(i.e. no displacment fixed) edges on the top and right. We will perform
%%this by calculating the equivalent surface tractions from the initial
%%stress field and apply them. Note that you cansimply check if they
%%correspond to the stresses given in the section before.

%----- Get the corresponding tractions
f_y_top = sum(Sig_o(top*2))/l_area;
f_x_right = sum(Sig_o(right*2-1))/d_area;
f_y_fracture = p_f;

%----- Assemble the force vector
f_s=AssembleTractionsOverLine(meshE,right,[f_x_right,0],'Axis') + ...
    AssembleTractionsOverLine(meshE,top,[0,f_y_top],'Axis') + ...
    AssembleTractionsOverLine(meshE,fracture_E,[0,f_y_fracture],'Axis') - ...
    Sig_o;

%% Assembling the different matrices and vectors

% Stiffness matrix
[K]=AssembleStiffnessMatrix(meshE,params,'Axis');

% mass matrix term
[Storage]=AssembleMassMatrix(meshP,1/M,'Axis');

% laplacian term
[C]=AssembleConductivityMatrix(meshP,kappa,'Axis');

% Coupling term 
[Ce]=AssembleCouplingMatrix(meshE,meshP,b,'Axis');

%% Fracture terms

% Volume compressibility
f_V = 1/2*beta*Vo;

% Dirichlet BC terms
I_pf = ones(length(fracture_P),1);

C_diri = C(:,fracture_P)*I_pf;

Ce_diri = Ce(:,fracture_P)*I_pf;

% Opening change fracture
f_udot = sparse(2*length(meshE.nodes(:,1)),1);
f_udot(fracture_E.*2)=1/2;

% Surface traction condition
f_p = AssembleTractionsOverLine(meshE,fracture_E,[0,1],'Axis');

%% Solving directly for the undrained solution at t=0 ! 
%%First we solve the undrained response. Where we already set the boundary
%%condition within the borehole to zero pressure and zero stress.
dt= 0.;         % time step

AA=Storage+dt*C;   % lower right part of the system matrix

%----- Getting the different number of DOF and assembling the total matrix
% getting DOF
ntot_P = length(C(:,1));
ntot_E = length(K(:,1));
ntot=ntot_E+ntot_P;

% Assemble the matrix
TotMat=sparse(ntot,ntot);
TotMat(1:ntot_E,1:ntot_E)=K;
TotMat(ntot_E+1:ntot,ntot_E+1:ntot)=-AA;
TotMat(ntot_E+1:ntot,1:ntot_E)=-Ce';
TotMat(1:ntot_E,ntot_E+1:ntot)=-Ce;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%----- Preparing the system
% Get free equations.
eq_free_u = setdiff([1:ntot_E],fixed_dof)';
eq_fixed_u = setdiff([1:ntot_E],eq_free_u)';
eq_free= [eq_free_u;eq_free_p+ntot_E];
eq_fixed = setdiff([1:ntot],eq_free)';

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%----- Setting the force vector to the initial conditions
% Note that you need to account for the initial pore pressure field which
% will influence the initial force vector.
Ftot=sparse(ntot,1);
Ftot(1:ntot_E) = f_s-Ce*p_fixed;
Ftot(ntot_E+1:ntot) = -AA*p_fixed;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%----- Solving the system
Undrained_sol=sparse(ntot,1);
Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

%----- Seperate displacements from pressure
U_u=sparse(ntot_E,1);
U_u(eq_free_u) = Undrained_sol(eq_free_u);
pp_u = sparse(ntot_P,1);
pp_u(eq_free_p) = Undrained_sol(ntot_E+eq_free_p);
pp_u(eq_fixed_p) = p_fixed(eq_fixed_p);

%% ----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [U_u(1:2:end) U_u(2:2:end)];
 
figure(fig_ind)
amp_fact=1e2;
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*amp_fact,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title('deformed mesh - undrained response')
fig_ind = fig_ind + 1;

% plotting the pressure profile
figure(fig_ind) 
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2),full(pp_u))
title('Pore pressure distribution - undrained response')
xlabel(' x [m] ');
ylabel(' y [m] ');
zlabel(' p_u [MPa]');
fig_ind = fig_ind + 1;

%% Calculate the projected stress
S = ProjectStress(meshE,'2D',[E_u nu_u],U_u);
%%
figure(fig_ind)
trisurf(meshE.connectivity(:,1:3),meshE.nodes(:,1),meshE.nodes(:,2),S(:,2))
title(' Horizontal stress over the domain ( Projected S_{xx} ) ');
xlabel(' x [m] ');
ylabel(' y [m] ');
zlabel(' S_{xx} [MPa] ');
fig_ind = fig_ind + 1;

%% Plotting the opening profile along the bottom line
ind_tip = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)==Radius);
L_f = meshE.nodes(ind_tip,1) + U_u(ind_tip * 2 - 1);

op_prof = (meshE.nodes(fracture_E,:)+...
            udisp(fracture_E,:))';
save_val = [L_f*(1-logspace(-9,0,500));
    (2/(E_u))*L_f*(p_f-Po)*(1-(1-logspace(-9,0,500)).^2).^(1/2)];
figure(fig_ind)
scatter(op_prof(1,:),op_prof(2,:)+d_area,'Displayname',...
    'Numerical solution (course code)')
% hold on;
% plot(L_f*(1-logspace(-9,0,500)),...
%     (2/(E_u))*L_f*(p_f-Po)*(1-(1-logspace(-9,0,500)).^2).^(1/2),...
%     'Displayname','Analytical solution','LineWidth',3)
title('Fracutre opening at the undrained state')
legend
xlabel(' x [m]');
ylabel(' w_u [m]');
fig_ind = fig_ind + 1;

% %% Compare undrained pore pressure singularity
% 
% %----- pressure undrained response \theta=0
% % define the analytical solution
% Analyticp_u = @(x) -B*KIc./((x-L_f)*pi()*2).^(1/2)+p_o;
% 
% % extract and re-arrange the numerical solution
% bottom_p = find(meshP.nodes(:,2)==-d_area);
% [p_bottom,ia]=sort(meshP.nodes(bottom_p,1));
% 
% % plotting the pore pressure
% figure(fig_ind)
% plot(p_bottom, pp_u(bottom_p(ia)),'-*b')
% hold on
% plot(p_bottom(p_bottom > L_f),Analyticp_u(p_bottom(p_bottom > L_f)),'-r')
% title('Pore pressure along x=-d_{area} - undrained response')
% xlabel(' x [m]');
% ylabel(' p_u [MPa]');
% fig_ind = fig_ind + 1;
% 
% % plotting the error on the pore pressure along the axis.
% figure(fig_ind)
% semilogy([L_f L_f],[10^(-3), 10^(3)],'--k','LineWidth',2.5,...
%     'Displayname','Crack front location')
% hold on;
% semilogy(p_bottom,abs(pp_u(bottom_p(ia))-Analyticp_u(p_bottom))./...
%     abs(Analyticp_u(p_bottom)),'b','LineWidth',2.5,...
%     'Displayname','Error')
% title('Relative error on p_u along x=-d_{area} - undrained response')
% legend
% xlabel(' x [m] ');
% ylabel(' E(p_u)_{rel} ');
% fig_ind = fig_ind + 1;
% 
% %% Compare undrained stress singularity ahead of crack tip
% 
% %----- pressure undrained response \theta=0
% % define the analytical solution
% AnalyticSxx_u = @(x) KIc./((x-L_f)*pi()*2).^(1/2);
% 
% % extract and re-arrange the numerical solution
% bottom_S = find(meshE.nodes(:,2)==-d_area);
% [S_bottom,is]=sort(meshE.nodes(bottom_S,1));
% 
% % plotting the pore pressure
% figure(fig_ind)
% plot(S_bottom, S(bottom_S(is),2),'-*b')
% hold on
% plot(S_bottom(S_bottom > L_f),AnalyticSxx_u(S_bottom(S_bottom > L_f)),'-r')
% title('Horizontal stress along x=-d_{area} - undrained response')
% xlabel(' x [m]');
% ylabel(' S_{xx} [MPa]');
% fig_ind = fig_ind + 1;
% 
% % plotting the error on the pore pressure along the axis.
% figure(fig_ind)
% semilogy([L_f L_f],[10^(-3), 10^3],'--k','LineWidth',2.5,...
%     'Displayname','Crack front location')
% hold on;
% semilogy(S_bottom,abs(S(bottom_S(is),1)-AnalyticSxx_u(S_bottom))./...
%     abs(Analyticp_u(S_bottom)),'b','LineWidth',2.5,...
%     'Displayname','Error')
% title('Relative error on S_{xx} along x=-d_{area} - undrained response')
% legend
% xlabel(' x [m] ');
% ylabel(' E(S_{xx})_{rel} ');
% fig_ind = fig_ind + 1;