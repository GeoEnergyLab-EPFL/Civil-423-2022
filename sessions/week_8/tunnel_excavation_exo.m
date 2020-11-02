%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Tunnel Excavation - Exo week 8                %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all
close all
%% Mesh generation
Radius = 1.;    % Radius of the tunnel
nrad = 90;       % number of segments of the circular face used
l_area = 50;    % length of observed area
d_area = 50;    % Depth of observed area
c_t = [0,-d_area];  % Center of the tunnel 

thetas=linspace(0.,1/2*pi,nrad);
node = [ 0., 0. ;
    l_area, 0. ;
    l_area, -d_area ;
    c_t(1)+Radius*cos(thetas)' , c_t(2)+Radius*sin(thetas)'];
node(end,1) = 0.;

edge = [];
for e=1:nrad+2,
    edge = [ edge ; e e+1 ];
end
edge = [ edge ; nrad+3 1 ];

%-- call mesh-gen.
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
    refine2(node,edge,[],[],4) ;
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = smooth2(meshL.nodes,...
    meshL.edge, meshL.connectivity,meshL.id); % Refine the mesh
%-- Transform into quadratic mesh
[meshQ.nodes,meshQ.connectivity] = Tri3ToTri6(meshL);
meshQ.edge = meshL.edge;
meshQ.id = meshL.id;

% Decide which mesh to use for elasticity (meshE) and pressure (meshP)
meshE = meshQ;
meshP = meshQ;

% Plot the mesh and higlight the boundary.
figure(1);
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;

% Controlling the CFL condition of the numerical solution
hmin=(pi/2)/nrad ;
dt_cfl=hmin^2

%% Material properties
%%All stiffnesses are given in [MPa] and correspond to a ohio Sandstone.

k=8.4e3;            % elastic drained bulk modulus [MPa]
g=6.8e3;            % shear modulus [MPa]
b=0.707692;         % Biot coefficient 
M=9.18478e3;        % Biot Modulus
k_u=k+b^2*M;        % undrained bulk modulus [MPa]
perm =0.137549e-3;  % permeability value adjusted such that c= k /(mu S)=1 
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

eta = b*(1-2*nu_u)/(2*(1-nu_u)); % poroelastic eta parameter

%% Boundary conditions
%%We will impose that there are no horizontal displacements on the left
%%boundary and no vertical displacements on the bottom boundary.
%%Additionally, the pressure within the tunnel will be set to zero. And
%%surface tractions on the right and the top will be set accordingly to the
%%initial stress field.

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%----- Boundaries
% find the corresponding nodes
bottom
left
top
right

%----- Displacements
% adjust the fixed dofs and nodes
fixed_dof_left
fixed_dof_bottom
fixed_nodes
fixed_dof

%----- Pressure
% nodes on tunnel. Note that it will be necessary to introduce a rounding
% as to really capture all the points. This is because we need to account
% for the numerical imprecision of the meshing.
tunnel

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%----- Plotting the boundarys onto the initial plot of the mesh
% Boundaries where displacements are fixed will be red, such that are fixed
% in pressure blue and the ones fixed in tractions green
hold on;
plot(meshE.nodes(bottom,1),meshE.nodes(bottom,2),'or')
hold on;
plot(meshE.nodes(left,1),meshE.nodes(left,2),'or')
hold on;
plot(meshP.nodes(tunnel,1),meshP.nodes(tunnel,2),'*b')
hold on;
plot(meshE.nodes(right,1),meshE.nodes(right,2),'sg')
hold on;
plot(meshE.nodes(top,1),meshE.nodes(top,2),'sg')
%% Initial condition
%%As an initial condition we will set a unifrom pressure field with
%%deviatoric stresses where s_x = -20 and s_y = -40 but without shear
%%stresses. Additionally we will define an initial pore pressure field
%%which is uniform as well.

%----- Initial stress field
% Set the initial stressfield in MPa
%
Po = 30; %mean compressive 
So = 10; % deviatoric compressive
Sig_unif = -[(Po-So);(Po+So);0];  %  [MPa]
% Assemble the stress vector at the nodes
Sig_o = SetStressField(meshE,'2D',...
    {[1:length(meshE.nodes(:,1))],Sig_unif});

%----- Initial pore pressure field
% Define the fixed and free nodes
eq_free_p = setdiff([1:length(meshP.nodes(:,1))],tunnel)';
eq_fixed_p = setdiff([1:length(meshP.nodes(:,1))],eq_free_p)';


% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Define the initial pore pressure field
p_o = 3; %[MPa]
p_fixed = sparse(length(meshP.nodes(:,1)),1);
p_fixed(eq_free_p)
p_fixed(eq_fixed_p)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%% Set boundary tractions
%%We will need to set the corresponding surface tractions on the "free"
%%(i.e. no displacment fixed) edges on the top and right. We will perform
%%this by calculating the equivalent surface tractions from the initial
%%stress field and apply them. Note that you cansimply check if they
%%correspond to the stresses given in the section before.

%----- Get the corresponding tractions
f_y_top = sum(Sig_o(top*2))/l_area;
f_x_right = sum(Sig_o(right*2-1))/d_area;


% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%----- Assemble the force vector
f_s

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%%Note: As we switch the stress field as a force we need to switch the
%%side of the equation. This means that we got a minus signe in front. 
%% Test on deformation
%%We will balance now the entire stress field by the corresponding boundary
%%tractions. As such one should observe that no displacements take place
%%and thus nothing happens.

% We need to redefine the boundary at the tunnel in function of the elastic
% nodes.
tunnel_E=find(meshE.nodes(:,2)<0 & ...
    round(sqrt((meshE.nodes(:,1)-c_t(1)).^2.+...
    (meshE.nodes(:,2)-c_t(2)).^2),4)==Radius);

% Define the corresopnding stresses (stress field + all the boundary
% tractions at free surfaces).
f_test=AssembleTractionsOverLine(meshE,right,[f_x_right,0],'2D') + ...
    AssembleTractionsOverLine(meshE,top,[0,f_y_top],'2D') - ...
    Sig_o + ...
    AssembleTractionsOverLine(meshE,tunnel_E,-Sig_unif(1:2),'2D');
%%Note 2: we need the minus in the the forces over the tunnel because the
%%outward directing normal there is negative.

% Stiffness matrix
[K]=AssembleStiffnessMatrix(meshE,params,'2D');

eq_to_solve=setdiff([1:length(meshE.nodes)*2],fixed_dof)';
% We get all the DOF's to solve for. Note that we only solve for elasticity
% so we solve for all DOF's without prescribed displacements.

% We solve now the system for the unknown displacements
u_aux = K(eq_to_solve,eq_to_solve)\f_test(eq_to_solve);

% And finally re-assemble the full solution
u_res=0*[1:length(meshE.nodes)*2]';
u_res(fixed_dof)=0; 
u_res(eq_to_solve)=u_aux;

%----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [u_res(1:2:end) u_res(2:2:end)];

%------ Plotting the mesh to check for displacements
figure(2)
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title('deformed mesh (equilibrium test)')

% check for small displacements
max_disp = max(u_res)
%% Assembling the different matrices

% elasticity
%%Already defined previously.

% mass matrix term
[Storage]=AssembleMassMatrix(meshP,1/M,'2D');

% laplacian term
[C]=AssembleConductivityMatrix(meshP,kappa,'2D');

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Coupling term  -> you have to modify ElementCouplingMatrix_exo
[Ce]=AssembleCouplingMatrix_exo(meshE,meshP,b,'2D');

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


%% Solving directly for the undrained solution at t=0 ! 
%%First we solve the undrained response. Where we already set the boundary
%%condition within the borehole to zero pressure and zero stress.
dt= 0.;         % time step



% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

AA  % lower right part of the system matrix

%----- Getting the different number of DOF and assembling the total matrix
% getting DOF
ntot_P
ntot_E
ntot

% Assemble the matrix
TotMat

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
Ftot(1:ntot_E) 
Ftot(ntot_E+1:ntot)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%----- Solving the system
Undrained_sol=sparse(ntot,1);
Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

%----- Seperate displacements from pressure
U_u=sparse(ntot_E,1);
U_u(eq_free_u) = Undrained_sol(eq_free_u);
pp_u = Undrained_sol(ntot_E+1:end);

%----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [U_u(1:2:end) U_u(2:2:end)];
 
figure(3)
amp_fact=1e2;
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*amp_fact,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title('deformed mesh - undrained response')

% plotting the pressure profile
figure(4) 
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2),full(pp_u))
title('Pore pressure distribution - undrained response')
xlabel(' r ');
ylabel(' z ');
zlabel(' p_u ');
%% Compare undrained displacements and pore pressure field
%%In this section we will compare the displacements at the top node of the
%%tunnel (i.e. at (0,-d_area+radius)) and the pore pressure profile along
%%the axis x = -d_area.

%---- Displacements
% Get the corresopnding displacements
top_node = find(meshE.nodes(:,1)==0 & meshE.nodes(:,2)==-d_area + Radius);


% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Get the undrained solution
AnalyticSolUx

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Calculate the error
Error_undrained_ur=abs((AnalyticSolUx - U_u(top_node*2))/AnalyticSolUx)

%----- pressure undrained response \theta=0
% define the analytical solution
Analyticp_u = @(x) 4*B*(1+nu_u)/3*Radius.^2./x.^2*So+p_o;

% extract and re-arrange the numerical solution
bottom_p = find(meshP.nodes(:,2)==-d_area);
[p_bottom,ia]=sort(meshP.nodes(bottom_p,1));

% plotting the pore pressure
figure(5)
plot(p_bottom, pp_u(bottom_p(ia)),'-*b')
hold on
plot(p_bottom,Analyticp_u(p_bottom),'-r')
title('Pore pressure along x=-d_{area} - undrained response')
xlabel(' r ');
ylabel(' p_u ');

% plotting the error on the pore pressure along the axis.
figure(6)
plot(p_bottom,abs(pp_u(bottom_p(ia))-Analyticp_u(p_bottom))./Analyticp_u(p_bottom))
title('Relative error on p_u along x=-d_{area} - undrained response')
xlabel(' r ');
ylabel(' E(p_u)_{rel} ');
%% Solving now the drained response for t > 0

% As the problem is sensitive to the time scale in the beginning we cannot
% use a constant time step here. Instead, we will predefine times when we
% want to get the solution and get the time step from the difference in the
% times. Investigate how the time step is built up to get a feeling how we
% split up the resolution.
t_k=[0]; 
for i = 1:7
    t_k = [t_k;
        logspace((i-6),(i-5),40)'];
end
t_k =unique(t_k);
final_time = max(t_k) % Final time

% prepare rhs
F_n=sparse(ntot,1);
%----- initialize solution vectors
n_step= length(t_k);        % numer of time steps
u_num_90=zeros(n_step,1);   % displacements at top node
u_num_90(1) = U_u(top_node*2);

U_n=U_u;
pp_n=pp_u;

%----- time loop (from two to the number of iterations chosen 
for i=2:n_step


% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    % calculate the time step
    dt
    
    % changes of the the total matrix
    TotMat

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    
    Dp=0.*pp_n;        % at nodes for fluid flow
    DU=0.*U_n;         % at nodes for elasticity
    D_x=[DU; Dp];      % solution vector
    
    % calculate the change in flux
    pn_lhs=dt*C(eq_free_p,:)*pp_n(:); 
    
    % Set the right hand side
    F_n(eq_free_p+ntot_E)=pn_lhs;
    
    % Solve the system
    D_x(eq_free)=TotMat(eq_free,eq_free)\F_n(eq_free);
    
    % Seperate elastic and flow components
    Dp(eq_free_p)=D_x(eq_free_p+(ntot_E));
    DU(eq_free_u)=D_x(eq_free_u);
    
    % Update pressure and displacement
    pp_n=pp_n+Dp;
    U_n=U_n+DU;
    
    % save displacement along uy of the top_node
    u_num_90(i) = U_n(top_node*2);

end
%% Visualizing the solution

%----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [U_n(1:2:end) U_n(2:2:end)];

figure(7) 
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*amp_fact,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title(' Deformed mesh at the final simulation time (drained) ');

  
%----- Plotting the pore pressure
figure(8) 
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2)...
    ,full(pp_n))
title(' Pore pressure distribution at the final simulation time');
xlabel('r');
ylabel('z');
ylabel('p');
%% Compare the stresses at large time
%%We will compare the stresses along the bottom boundary to the known
%%analytical solution to see if our approximation is valid. Additionally,
%%we will compare the evolutino of the displacement of the top node to the
%%large and early time solutions.

%----- Getting the analytical solutions (displacements)
% For the undrained displ we already obtained it previously.
AnalyticSolUx_d=-(Po.*Radius./(2.*g)+(3-4*nu).*Radius.*So./(2.*g));

%----- Plotting disp evolution
figure(9) 
semilogx(t_k(1:n_step),u_num_90(1:end),'.-b')
hold on;
semilogx([10^(-6),t_k(n_step)],[AnalyticSolUx_d,AnalyticSolUx_d],'--k')
hold on
semilogx([10^(-6),t_k(n_step)],[AnalyticSolUx,AnalyticSolUx],'-.r')
title('Time evolution of the vertical displacements at the borehole wall at \theta = 90^\circ node');
xlabel('t');
ylabel('u_r');

%----- Calculating the errors
% undrained
error_undrained=(u_num_90(1)-AnalyticSolUx)/AnalyticSolUx
% drained
error_drained=(u_num_90(end)-AnalyticSolUx_d)/AnalyticSolUx_d

%---- Get the corresponding stresses

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%---- Compute the corresponding total stresses 
Sigma

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%---- Define the analytical solution (stresses)
% Define the solution
AnalyticSolSx= @(x) -Po*(1-Radius^2./x.^2)+So*(1-4.*Radius^2./x.^2 ...
+3.*Radius^4./x.^4) + eta.*p_o*(1-Radius^2./x.^2);
% Getting numerical values
A_sol = AnalyticSolSx(sort(meshE.nodes(bottom,1)));

 
% Extract numerical solution
[bootom_stress,ia]=sort(meshE.nodes(bottom,1));

% Plotting the stresses
figure(10)
plot(bootom_stress,Sigma(bottom(ia),1),'-b')  % here you need to add Sig_unif
hold on
plot(sort(meshE.nodes(bottom,1)),...
     A_sol,'-r');
title('Horizontal stresses along x=-d_{area} at the final simulation time');
xlabel('x');
ylabel('S_{xx}');