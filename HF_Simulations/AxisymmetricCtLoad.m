%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Axissymmetry Elasticity fracture              %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all
close all

%% Material properties
%%All stiffnesses are given in [MPa] and correspond to a ohio Sandstone.

k=22e3;            % elastic drained bulk modulus [MPa]
g=6.e3;            % shear modulus [MPa]

nu=(1-2*(g)/(3*k))/(2+2*(g)/(3*k));         % poisson's ratio
E=(9*(k))/(1+3*(k)/(g));                    % youngs modulus
params = [E, nu];

Ep=E/(1-nu*nu);

%% Mesh generation
Radius = 2.0
l_area = floor(10);    % length of observed area
d_area = floor(10);    % Depth of observed area

nrad = 100;     % number of segments in the tip
nfr = 20;     % number of segments in the fracture

tip_zone = flip([Radius-flip(logspace(-8,-1,nrad));
            -d_area*ones(1,nrad)],2)'; % nodes of the tip zone

        external_approach = [Radius+flip(logspace(-8,0.,nrad));
            -d_area*ones(1,nrad)]'; % nodes ahead of the tip zone
        
fracture_zone = flip([linspace(0,Radius-0.1,nfr);
            -d_area*ones(1,nfr)],2)';   % nodes from injection point to the
                                        % tip zone

node = unique([
     0., 0. ;
    l_area, 0. ;
    l_area, -d_area;
    external_approach;
    Radius, -d_area;
    tip_zone;
    fracture_zone],'rows','stable');
node(end,1) = 0.; % combine all nows
% SWITCH ALL COORDINATES by d_area up
node(:,2)=node(:,2)+d_area


edge = [];
for e=1:length(node)
    edge = [ edge ; e e+1 ];
end
edge(end,2) = 1; % get the edges


hmax=2.; hmin=0.01;
hfun=@(x) hmax-(hmax-hmin)*exp(-0.02*((x(:,1)-Radius).^2+x(:,2).^2))

%-- call mesh-gen.
[meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
    refine2(node,edge,[],[],hfun) ;
 [meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
     smooth2(meshL.nodes,meshL.edge, meshL.connectivity,meshL.id); 
    % Refine the mesh

%-- Transform into quadratic mesh
[meshQ.nodes,meshQ.connectivity] = Tri3ToTri6(meshL);
meshQ.edge = meshL.edge;
meshQ.id = meshL.id;

% Decide which mesh to use for elasticity (meshE) and pressure (meshP)
meshE = meshQ; % Quadratic elasticity

% Plot the mesh and higlight the boundary.
fig_ind = 1;
figure(fig_ind);
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
fig_ind = fig_ind + 1;


%% Boundary conditions

%----- Boundary condition on displacement 
% find the corresponding nodes
bottom = find(meshE.nodes(:,2)==0. & meshE.nodes(:,1)>=Radius);
bottom_crack= find(meshE.nodes(:,2)==0. & meshE.nodes(:,1)<Radius);
%bottom_P = find(meshP.nodes(:,2)==-d_area & meshP.nodes(:,1)>=Radius);
left = find(meshE.nodes(:,1)==0.);
top = find(meshE.nodes(:,2)==d_area);
right = find(meshE.nodes(:,1)==l_area);

%----- Displacements
% adjust the fixed dofs and nodes
fixed_dof_left = [2*(left-1)+1]; % u_x
fixed_dof_bottom = [2*bottom];  %u_y
fixed_nodes = unique([bottom;left]);
fixed_dof = unique([fixed_dof_left;fixed_dof_bottom]);


%----- Plotting the boundarys onto the initial plot of the mesh
% Boundaries where displacements are fixed will be red, such that are fixed
% in pressure blue and the ones fixed in tractions green
hold on;
plot(meshE.nodes(bottom,1),meshE.nodes(bottom,2),'or')
hold on;
plot(meshE.nodes(left,1),meshE.nodes(left,2),'or')
hold on;
plot(meshE.nodes(right,1),meshE.nodes(right,2),'sg')
hold on;
plot(meshE.nodes(top,1),meshE.nodes(top,2),'sg')

hold on;
plot(meshE.nodes(bottom_crack,1),meshE.nodes(bottom_crack,2),'sb')


%% 

%----- Far-field stress field
% Set the initial stressfield in MPa
%
% here we will just superimpose later on.
S_rr = 0;
S_zz = 0; 
Sig_unif = -[S_rr;S_zz;0;S_zz];  %  [MPa]
% Assemble the stress vector at the nodes
Sig_o = SetStressField(meshE,'Axis',...
    {[1:length(meshE.nodes(:,1))],Sig_unif});



%% Set boundary tractions

p_f = 1

p_net=p_f-S_zz
%----- Get the corresponding tractions

f_y_fracture = p_f;

%----- Assemble the force vector
f_s= AssembleTractionsOverLine(meshE,bottom_crack,[0,f_y_fracture],'Axis');
 

%% Assembling the different matrices and vectors

% Stiffness matrix
[K]=AssembleStiffnessMatrix(meshE,params,'Axis');



%% Solve elastic problem
ntot_E=2*length(meshE.nodes)

eq_free_u = setdiff([1:ntot_E],fixed_dof)';
eq_fixed_u = setdiff([1:ntot_E],eq_free_u)';

Ftot = f_s+Sig_o;

solU=sparse(ntot_E,1);
solU(eq_free_u)=K(eq_free_u,eq_free_u)\Ftot(eq_free_u);

%% ----- Plotting the deformed mesh
% we first need to reshape the mesh
udisp = [solU(1:2:end) solU(2:2:end)];
 
figure(fig_ind)
amp_fact=1e2;
plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
hold on;
plotmesh(meshE.nodes+udisp*amp_fact,meshE.connectivity(:,1:3),[.8 .2 .2],'none')
title('deformed mesh ')
fig_ind = fig_ind + 1;

%% plot crack width

rc=meshE.nodes(bottom_crack,1);

w_n=solU(2*bottom_crack);

w_true=4*p_net/Ep/pi*(sqrt(Radius^2.0-rc.^2.0))
figure(fig_ind)
plot(rc,w_true,'.r')
hold on
plot(rc,w_n,'.b')
fig_ind = fig_ind + 1;

%% plot shear
u_x_n=solU(2*(bottom_crack-1)+1);
figure(fig_ind)
semilogx( Radius-rc, (Ep*(u_x_n-u_x_n(end))./(Radius-rc)),'.r')
ylim([-20,5])
%loglog( Radius-rc, Ep*(u_x_n-u_x_n(end))./(Radius-rc),'.r')

%plot( Radius-rc,-u_x_n,'.')
fig_ind = fig_ind + 1;


rb=meshE.nodes([bottom; bottom_crack],1);
u_x_all=solU(2*([bottom; bottom_crack]-1)+1);
% figure(fig_ind)
 %plot(rb,u_x_all,'.r')
% 
% fig_ind = fig_ind + 1;

%% rel error

rel_error=abs(w_n./w_true-1);
semilogy(rc,rel_error,'.b')

figure(fig_ind)
abs_error=abs(w_n-w_true);
semilogy(rc,abs_error,'.b')
fig_ind=fig_ind+1
%% Calculate the projected stress
S = ProjectStress(meshE,'Axis',[E nu],solU);


%% plot stress
xx=meshE.nodes([bottom_crack; bottom],1);

s_zz_n=S([bottom_crack; bottom],2)
s_rr_n=S([bottom_crack; bottom],1)
s_rz_n=S([bottom_crack; bottom],3)
s_tt_n=S([bottom_crack; bottom],4)

figure(fig_ind)
plot(xx,s_rr_n,'.')
hold on
plot(xx,s_zz_n,'.')
fig_ind=fig_ind+1

%% just ahead of tip 
x_ahead=meshE.nodes(bottom,1)-Radius;

s_zz_n=S([bottom],2)
s_rr_n=S([bottom],1)
s_rz_n=S([bottom],3)
s_tt_n=S([ bottom],4)

figure(fig_ind)
loglog(x_ahead,s_zz_n,'.')
hold on
loglog(x_ahead,s_rr_n,'.')

k_i=(2*p_net*(Radius^0.5))/sqrt(pi);
s_exp=k_i./sqrt(2.0*pi.*sort(x_ahead));
loglog(sort(x_ahead/Radius),s_exp,'-r')
fig_ind=fig_ind+1

%% plotting s_rr-s_zz
figure(fig_ind)
 semilogx(x_ahead/Radius,s_rr_n-s_zz_n,'.')
 ylim([-10,10])
hold on
fig_ind=fig_ind+1

%% plottinf s_xx inside crack
figure(fig_ind)
semilogx(-(rc-Radius)/Radius,S([bottom_crack],1),'.')
ylim([-20,20])
hold on
fig_ind=fig_ind+1

%%

%%
figure(fig_ind)
trisurf(meshE.connectivity(:,1:3),meshE.nodes(:,1),meshE.nodes(:,2),S(:,1)-S(:,2))
title(' S_{rr}-S_{zz} over the domain ( Projected) ');
xlabel(' x [m] ');
ylabel(' y [m] ');
zlabel(' S_{rr}-S_{zz} [MPa] ');
fig_ind = fig_ind + 1;
%%
figure(fig_ind)
trisurf(meshE.connectivity(:,1:3),meshE.nodes(:,1),meshE.nodes(:,2),S(:,2))
title(' S_zz over the domain ( Projected S_{zz} ) ');
xlabel(' x [m] ');
ylabel(' y [m] ');
zlabel(' S_{zz} [MPa] ');
fig_ind = fig_ind + 1;
