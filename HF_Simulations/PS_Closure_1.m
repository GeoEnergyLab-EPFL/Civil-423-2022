%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Plain strain poroelastic closure              %%%%%%%
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
Sig_unif = -[(Po-So);(Po+So);0];  %  [MPa]
% Assemble the stress vector at the nodes
Sig_o = SetStressField(meshE,'2D',...
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
f_s=AssembleTractionsOverLine(meshE,right,[f_x_right,0],'2D') + ...
    AssembleTractionsOverLine(meshE,top,[0,f_y_top],'2D') + ...
    AssembleTractionsOverLine(meshE,fracture_E,[0,f_y_fracture],'2D') - ...
    Sig_o;

%% Assembling the different matrices and vectors

% Stiffness matrix
[K]=AssembleStiffnessMatrix(meshE,params,'2D');

% mass matrix term
[Storage]=AssembleMassMatrix(meshP,1/M,'2D');

% laplacian term
[C]=AssembleConductivityMatrix(meshP,kappa,'2D');

% Coupling term 
[Ce]=AssembleCouplingMatrix(meshE,meshP,b,'2D');

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
f_p = AssembleTractionsOverLine(meshE,fracture_E,[0,1],'2D');

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
hold on;
plot(L_f*(1-logspace(-9,0,500)),...
    (2/(E_u))*L_f*(p_f-Po)*(1-(1-logspace(-9,0,500)).^2).^(1/2),...
    'Displayname','Analytical solution','LineWidth',3)
% hold on;
% scatter(sim_1(1,:),...
%     sim_1(2,:),...
%     'Displayname','Numerical solution (course code) 2','LineWidth',3)
% hold on;
% plot(sim_1_ana(1,:),...
%     sim_1_ana(2,:),...
%     'Displayname','Analytical solution 2','LineWidth',3)
% hold on;
% scatter(sim_2(1,:),...
%     sim_2(2,:),...
%     'Displayname','Numerical solution (course code) 3','LineWidth',3)
% hold on;
% plot(sim_2_ana(1,:),...
%     sim_2_ana(2,:),...
%     'Displayname','Analytical solution 3','LineWidth',3)
title('Fracutre opening at the undrained state')
legend
xlabel(' x [m]');
ylabel(' w_u [m]');
fig_ind = fig_ind + 1;

%% Compare undrained pore pressure singularity

%----- pressure undrained response \theta=0
% define the analytical solution
Analyticp_u = @(x) -B*KIc./((x-L_f)*pi()*2).^(1/2)+p_o;

% extract and re-arrange the numerical solution
bottom_p = find(meshP.nodes(:,2)==-d_area);
[p_bottom,ia]=sort(meshP.nodes(bottom_p,1));

% plotting the pore pressure
figure(fig_ind)
plot(p_bottom, pp_u(bottom_p(ia)),'-*b')
hold on
plot(p_bottom(p_bottom > L_f),Analyticp_u(p_bottom(p_bottom > L_f)),'-r')
title('Pore pressure along x=-d_{area} - undrained response')
xlabel(' x [m]');
ylabel(' p_u [MPa]');
fig_ind = fig_ind + 1;

% plotting the error on the pore pressure along the axis.
figure(fig_ind)
semilogy([L_f L_f],[10^(-3), 10^(3)],'--k','LineWidth',2.5,...
    'Displayname','Crack front location')
hold on;
semilogy(p_bottom,abs(pp_u(bottom_p(ia))-Analyticp_u(p_bottom))./...
    abs(Analyticp_u(p_bottom)),'b','LineWidth',2.5,...
    'Displayname','Error')
title('Relative error on p_u along x=-d_{area} - undrained response')
legend
xlabel(' x [m] ');
ylabel(' E(p_u)_{rel} ');
fig_ind = fig_ind + 1;

%% Compare undrained stress singularity ahead of crack tip

%----- pressure undrained response \theta=0
% define the analytical solution
AnalyticSxx_u = @(x) KIc./((x-L_f)*pi()*2).^(1/2);

% extract and re-arrange the numerical solution
bottom_S = find(meshE.nodes(:,2)==-d_area);
[S_bottom,is]=sort(meshE.nodes(bottom_S,1));

% plotting the pore pressure
figure(fig_ind)
plot(S_bottom, S(bottom_S(is),2),'-*b')
hold on
plot(S_bottom(S_bottom > L_f),AnalyticSxx_u(S_bottom(S_bottom > L_f)),'-r')
title('Horizontal stress along x=-d_{area} - undrained response')
xlabel(' x [m]');
ylabel(' S_{xx} [MPa]');
fig_ind = fig_ind + 1;

% plotting the error on the pore pressure along the axis.
figure(fig_ind)
semilogy([L_f L_f],[10^(-3), 10^3],'--k','LineWidth',2.5,...
    'Displayname','Crack front location')
hold on;
semilogy(S_bottom,abs(S(bottom_S(is),1)-AnalyticSxx_u(S_bottom))./...
    abs(Analyticp_u(S_bottom)),'b','LineWidth',2.5,...
    'Displayname','Error')
title('Relative error on S_{xx} along x=-d_{area} - undrained response')
legend
xlabel(' x [m] ');
ylabel(' E(S_{xx})_{rel} ');
fig_ind = fig_ind + 1;

%% Solving now the drained response for t > 0

t_k=[0]; 
for i = 1:10
    t_k = [t_k;
        logspace((i-8),(i-7),20)'];
end
t_k =unique(t_k);
final_time = max(t_k) % Final time

%----- initialize solution vectors
n_step= length(t_k);        % numer of time steps
r_change = zeros(n_step,1);
r_change(1) = meshE.nodes(ind_tip,1) + U_u(ind_tip * 2 - 1);

step_index = 1;
opening_profile = zeros((floor(n_step/step_index)+2)*2,length(fracture_E));
udisp = [U_u(1:2:end) U_u(2:2:end)];
opening_profile(1:2,:) = (meshE.nodes(fracture_E,:)+udisp(fracture_E,:))';

pressure_bottom = zeros(n_step,length(bottom_p));
pressure_bottom(1,:) = pp_u(bottom_p(ia));

inj_point_op = zeros(n_step,1);
ind_inj = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)==0);
inj_point_op(1) = U_u(ind_inj * 2);

fracture_pressure = zeros(n_step,1);
fracture_pressure(1) = p_f;
p_along_frac = zeros(n_step,length(fracture_P));
p_along_frac(1,:) = pp_u(fracture_P);

U_n=U_u;
pp_n=pp_u;
pp_n(ntot_P+1) = p_f;

V_frac = zeros(n_step,1);
V_frac(1) = Vo;

disp("t_(lim) = " + num2str((l_area-Radius)^2/c))
disp("t^*_(lim) = " + num2str((l_area-Radius)^2/Radius^2))

%% Adapt the matrix system

% Assemble the matrix
TotMat=sparse(ntot+1,ntot+1);
TotMat(1:ntot_E,1:ntot_E)=K;
TotMat(ntot_E+1:ntot,ntot_E+1:ntot)=-AA;
TotMat(ntot_E+1:ntot,1:ntot_E)=-Ce';
TotMat(1:ntot_E,ntot_E+1:ntot)=-Ce;
% Add the closure specific terms
TotMat(ntot+1,1:ntot_E) = -f_udot';
TotMat(ntot+1,ntot+1) = -f_V;
TotMat(1:ntot_E,ntot+1) = -Ce_diri - f_p;
TotMat(ntot_E+1:ntot,ntot+1) = -C_diri;

% Generate force vector
F_n = sparse(ntot+1,1);

%% Adapt free equations
% Impose pressure on boundary nodes
eq_free_p = setdiff([1:ntot_P+1]',fracture_P);
%eq_free_p = [1:ntot_P+1]';
helper = length(eq_free_p)-1;
eq_free = [eq_free_u;ntot_E+eq_free_p];
%% ----- time loop (from two to the number of iterations chosen 
i = 1;
while i<n_step && V_frac(i)>0
    
    i = i+1;
    
    % calculate the time step
    dt=t_k(i)-t_k(i-1);
    
    % changes of the the total matrix
    TotMat(ntot_E+1:ntot,ntot_E+1:ntot)=-(Storage+dt*C);
    TotMat(ntot_E+1:ntot,ntot+1) = -dt*C_diri;
    TotMat(ntot+1,ntot+1) = -f_V;
   
    Dp=0.*pp_n;        % at nodes for fluid flow
    DU=0.*U_n;         % at nodes for elasticity
    D_x=[DU; Dp];      % solution vector
    
    % calculate the change in flux
    pn_lhs=dt.*[C(eq_free_p(1:helper),:) ...
        C_diri(eq_free_p(1:helper))]*pp_n(:); 
    
    % Set the right hand side
    %F_tot(eq_free_p+ntot_E)=pn_lhs;
    F_n([eq_free_p(1:helper)+ntot_E])=pn_lhs;
    
    % Solve the system
    D_x(eq_free)=TotMat(eq_free,eq_free)\F_n(eq_free);
    
    % Seperate elastic and flow components
    Dp(eq_free_p)=D_x(eq_free_p+ntot_E);
    DU(eq_free_u)=D_x(eq_free_u);
    
    % Update pressure and displacement
    pp_n=pp_n+Dp;
    p_along_frac(i,:) = pp_n(fracture_P);
    pp_n(fracture_P) = pp_n(end);
    U_n=U_n+DU;
    
    % Fluid pressure in the fracture
    fracture_pressure(i) = pp_n(end);
    
    % Fracture volume recalculation
    udisp = [U_n(1:2:end) U_n(2:2:end)];
    opening = (meshE.nodes(fracture_E,:)+...
            udisp(fracture_E,:));
    opening(:,2) = opening(:,2) + d_area;
    opening = sortrows(opening,1);
    %To adapt here with trapz
    id = find(opening(:,2)>0);
    V_frac(i) = 4*trapz(opening(id,1)-opening(id(1),1),...
            opening(id,2));
    if V_frac(i) > Vo
        disp("Error")
    end
    
    figure(20)
    plot(opening(:,1),opening(:,2))
    pause(2)
    %----- Plotting the pore pressure
    figure(9) 
    trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2)...
        ,full(pp_n(1:ntot_P)))
    title(' Pore pressure distribution at the final simulation time');
    xlabel('x');
    ylabel('y');
    zlabel('p');
    view(120,45)
    pause(1)

    
    % Save the change in radius
    r_change(i) = meshE.nodes(ind_tip,1) + U_n(ind_tip * 2 - 1);
    
    % Save the displacement at the tip
    inj_point_op(i) = meshE.nodes(ind_inj,2)+...
            udisp(ind_inj,2) + d_area;
    
    % Save disp profile
    if mod(i,step_index) == 0
        opening_profile((floor(i/step_index)+1)*2-1:(floor(i/step_index)+1)*2,:) = ...
            (meshE.nodes(fracture_E,:)+...
            udisp(fracture_E,:))';
    elseif i == n_step
        opening_profile(end-1:end,:) = ...
            (meshE.nodes(fracture_E,:)+...
            udisp(fracture_E,:))';
    end

    f_V = 1/2*beta*V_frac(i);
    disp("The remaining volume is: " + num2str(V_frac(i)))
    disp("The fracture pressure is is: " + num2str(fracture_pressure(i)))
    disp("The displacement of the injection point is: " + ...
        inj_point_op(i))
    
    pressure_bottom(i,:) = pp_n(bottom_p(ia));

    disp("Solved for time t = " + num2str(t_k(i)))
    
end

%% Trying different stuff

Q = ProjectFlux(meshP,'2D',perm,pp_n(1:end-1));

figure(fig_ind)
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2),Q(:,1))
fig_ind = fig_ind + 1;

figure(fig_ind)
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2),Q(:,2))
fig_ind = fig_ind + 1;


%% Visualizing the solution

%----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [U_n(1:2:end) U_n(2:2:end)];

figure(8) 
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*amp_fact,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title(' Deformed mesh at the final simulation time (drained) ');

  
%----- Plotting the pore pressure
figure(9) 
trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2)...
    ,full(pp_n(1:ntot_P)))
title(' Pore pressure distribution at the final simulation time');
xlabel('r');
ylabel('z');
ylabel('p');

%%
[~,ind] = min(abs(inj_point_op))
% figure(10)
% semilogx(t_k(1:ind).*c./r_change(1:ind).^2, r_change(1:ind)./ L_f,...
%     'LineWidth',2.5)
% title('Evaluation of radius in time (disp of tip point)');
% xlabel('t^*');
% ylabel(' R(t) / R(t=0) ');
ind = length(inj_point_op);
ind_pic = [1 (1:floor(n_step/step_index))*step_index n_step];

figure(11)
for i = 1:length(opening_profile(:,1))/2
    scatter(opening_profile(i*2-1,:)/r_change(ind_pic(i)),...
        (opening_profile(i*2,:)+d_area)*g/(p_f*r_change(ind_pic(i))),...
        'DisplayName',"t^* = " + num2str(c * t_k(ind_pic(i)) / ...
        r_change(ind_pic(i))^2))
    hold on;
end
hold off;
legend
title('Normalized opening for different dimensionless times t^* = ct/R(t)^2');
ylabel('u_y * G / (p_f R(t))');
xlabel(' r / R(t) ');

figure(12)
for i = 1:length(opening_profile(:,1))/2
    [~, argmax] = max(abs(opening_profile(i*2,:)+d_area));
    scatter(opening_profile(i*2-1,:)./r_change(ind_pic(i)),...
        (opening_profile(i*2,:)+d_area)./(opening_profile(i*2,argmax)+d_area),...
        'DisplayName',"t^* = " + num2str(c * t_k(ind_pic(i)) / ...
        r_change(ind_pic(i))^2))
    hold on;
end
hold off;
legend
title('normalized opening for different dimensionless times t^* = ct/R(t)^2');
ylabel('u_y / max(u_y)');
xlabel(' r / R(t) ');

figure(13)
semilogx(t_k(1:ind-1).*c./r_change(1:ind-1).^2,...
    inj_point_op(1:ind-1)*g./(p_f*r_change(1:ind-1)),'LineWidth',2.5)
title('normalized displacement at the injection point');
ylabel('u_y(0,0,t)*g/(p_f R(t))');
xlabel('t^*');

figure(14)
semilogx(t_k(1:ind-1).*c./r_change(1:ind-1).^2,...
    fracture_pressure(1:ind-1),'LineWidth',2.5)
title('normalized displacement at the injection point');
ylabel('p_f');
xlabel('t^*');

figure(15)
semilogx(t_k(1:ind-1).*c./r_change(1:ind-1).^2,...
    V_frac(1:ind-1),'LineWidth',2.5)
title('normalized displacement at the injection point');
ylabel('V_frac');
xlabel('t^*');
%%
tstarlast = t_k(ind-2:ind-1).*c./r_change(ind-2:ind-1).^2
save('opening_profile.mat','opening_profile');
save('r_change.mat','r_change');
save('last_time.mat','tstarlast');
% 
% %% Compare the stresses at large time
% %%We will compare the stresses along the bottom boundary to the known
% %%analytical solution to see if our approximation is valid. Additionally,
% %%we will compare the evolutino of the displacement of the top node to the
% %%large and early time solutions.
% 
% %----- Getting the analytical solutions (displacements)
% % For the undrained displ we already obtained it previously.
% AnalyticSolUx_d=-(Po.*Radius./(2.*g)+(3-4*nu).*Radius.*So./(2.*g));
% 
% %----- Plotting disp evolution
% figure(9) 
% semilogx(t_k(1:n_step),u_num_90(1:end),'.-b')
% hold on;
% semilogx([10^(-6),t_k(n_step)],[AnalyticSolUx_d,AnalyticSolUx_d],'--k')
% hold on
% semilogx([10^(-6),t_k(n_step)],[AnalyticSolUx,AnalyticSolUx],'-.r')
% title('Time evolution of the vertical displacements at the borehole wall at \theta = 90^\circ node');
% xlabel('t');
% ylabel('u_r');
% 
% %----- Calculating the errors
% % undrained
% error_undrained=(u_num_90(1)-AnalyticSolUx)/AnalyticSolUx
% % drained
% error_drained=(u_num_90(end)-AnalyticSolUx_d)/AnalyticSolUx_d
% 
% %---- Get the corresponding stresses
% % 
% Sigma = ProjectStress(meshE,'2D',params,U_n);
% 
% %---- Define the analytical solution (stresses)
% % Define the solution
% AnalyticSolSx= @(x) -Po*(1-Radius^2./x.^2)+So*(1-4.*Radius^2./x.^2 ...
% +3.*Radius^4./x.^4) + eta.*p_o*(1-Radius^2./x.^2);
% % Getting numerical values
% A_sol = AnalyticSolSx(sort(meshE.nodes(bottom,1)));
% 
% % undrained solution to code !
% % AnalyticSolSx_u= @(x) -Po*(1-Radius^2./x.^2)+So*(1-4.*Radius^2./x.^2 ...
% % +3.*Radius^4./x.^4) + eta.*p_o*(1-Radius^2./x.^2);
% 
% % Extract numerical solution
% [bootom_stress,ia]=sort(meshE.nodes(bottom,1));
% 
% % Plotting the stresses
% figure(10)
% plot(bootom_stress,Sigma(bottom(ia),1)+Sig_unif(1)-b*(pp_n(bottom(ia))-p_o),'-b')  % here you need to add Sig_unif
% hold on
% plot(sort(meshE.nodes(bottom,1)),...
%      A_sol,'-r');
% title('Horizontal stresses along x=-d_{area} at the final simulation time');
% xlabel('x');
% ylabel('S_{xx}');
% 


%%
%     %----- Plotting the deformed mesh
%     % we first need to reshape the mesh
%     udisp = [U_n(1:2:end) U_n(2:2:end)];
%     
%     if mod(i,25) == 0
%         figure(11)
%         amp_fact=1e2;
%         plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
%         hold on;
%         plotmesh(meshE.nodes+udisp*amp_fact,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
%         title('deformed mesh - undrained response')
% 
%         % plotting the pressure profile
%         figure(12) 
%         trisurf(meshP.connectivity(:,1:3),meshP.nodes(:,1),meshP.nodes(:,2),full(pp_n))
%         title('Pore pressure distribution - undrained response')
%         xlabel(' r ');
%         ylabel(' z ');
%         zlabel(' p_u ');      
%         
%         %pause;
%     end