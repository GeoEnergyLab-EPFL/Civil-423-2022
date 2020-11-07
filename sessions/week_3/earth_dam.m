% Earth dams unconfined flow
% with fixed point iterations with under-relaxation
%% Mesh

% geometry according to figures 1 and 2
b=8;hd=40;hw=35;m=2;m1=2; 

% define the dam vertices coordinates
% The origin of the coordinate system is on the bottom-left corner
xtl=m*hd;  % crest left x
xtr=m*hd+b; % crest right x
x_ext=m*hd+b+m1*hd; % toe dam x

node = [0,0 ;xtl,hd ;xtr,hd ;x_ext,0 ];

edge = [];
for e=1:3,
edge = [ edge ; e e+1 ];
end
edge = [ edge ; 4  1 ];
 
% call mesh-generator MESH2D - download it at https://ch.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation
 opts.kind = 'delfront';
% opts.rho2 = +1.0 ;
% opts.siz1 = 1.33;
% opts.siz2 = 1.3;
 
 h_x = 2; %% Max elt area -> control the refinement of the mesh here
 [mesh.nodes,mesh.edge, mesh.connectivity,mesh.id] = refine2(node,edge,[],opts,h_x) ; 

 % plotting mesh
 figure(1);
 plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
 ne_t=length(mesh.connectivity(:,1)); 
 
%%  Boundary Conditions 

% find the nodes belonging to each boundary
left_edge=find(abs(mesh.nodes(:,1)/m-mesh.nodes(:,2))<0.001);
right_edge=find(abs((x_ext-mesh.nodes(:,1))/m1-mesh.nodes(:,2))<0.001);
bottom_edge=find(mesh.nodes(:,2)==0);
top_edge=find(mesh.nodes(:,2)==hd);

% plotting boundary nodes
figure(2)
plot(mesh.nodes(left_edge,1),mesh.nodes(left_edge,2),'ob');
hold on
plot(mesh.nodes(right_edge,1),mesh.nodes(right_edge,2),'ob');
hold on
plot(mesh.nodes(bottom_edge,1),mesh.nodes(bottom_edge,2),'ob');
hold on
plot(mesh.nodes(top_edge,1),mesh.nodes(top_edge,2),'ob');

% Here you have to impose the boundary conditions of the problem

% <<<Reminder: h = p/gamma_w + y>>>

% You have to define now the boundary condition at the upstream face (left side)

h_left=--------

% Now the boundary condition at the downstream face (right side)

h_right=---------

% And now the boundary condition at the crest of the dam (top)

h_top=----------

% Deleting duplication of fixed nodes (just in case). 
[nodes_fixed, ia, ic]=unique([left_edge;right_edge;top_edge]) ; % be careful there is duplicate-> use unique

% Defining fixed (set) piezometric head
h_all=[h_left ; h_right; h_top];
h_set=h_all(ia);

%% Fixed Point Interation Scheme
 
k_o=[1.]; % Hydraulic conductivity of material, [L/s]
k_current=k_o; % Initialization of hydraulic conductivity vector 
h_res=0*[1:length(mesh.nodes)]'; % Initialization of piezometric head solution 
beta=0.75; % Under relaxation coefficient (1 if none)
eps_perm=0.01; % Parameter that controls sharpness of relative permeability function
tolerance = 1.e-6; % Tolerance for iteration
err=1; % Initialization of error
iter_max=30; % Max. number of iterations
k=0; % Itertion step

while (k<iter_max) && (err>tolerance)
 k=k+1;
 
[C] = AssembleConductivityMatrix(mesh,k_current,'2D');

% First, we have to separate the nodes where the piezometric head is 
% unknown from nodes where the piezometric head is fixed (because of the 
% boundary conditions).

eq_to_solve=setdiff([1:length(mesh.nodes)],nodes_fixed)';

% Now, you have to complete the code and solve the nonlinear system 

% Compute the force vector
f = -------

% Compute the unknown piezometric heads at each step
h_aux = -------

% Defining solution of the previous step
h_res_1=h_res;

% Compute the solution of the current step
h_res(nodes_fixed)=------
h_res(eq_to_solve)=-------

% Relax the solution of the current step by using under-relaxation
h_res=--------

% Compute pressure from solution of head
gammaw=10;
p_res = ---------

% Compute the relative permeability of each element
[kr]=ComputeRelPerm(mesh,p_res,eps_perm);

% Find the non-fully saturated elements
[r,c,v]=find(------);

% Assign a new tag to each non-fully saturated element
new_mat=1+[1:length(r)]';
mesh.id=ones(length(mesh.connectivity(:,1)),1);
mesh.id(r)=------

% Update hydraulic conductivity vector 
k_current=---------

% Compute the error of the current step
err=max(abs(h_res-h_res_1)./h_res);

disp(['max relative change at its ', num2str(k), ' is : ',num2str(err)]);

end
 

%% plotting solution on the grid.

figure(3)
title(' Num results - head');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),h_res)
  daspect([1 1 1])
 xlabel('x');
 ylabel('y');

figure(4)
title(' Num results - p');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),p_res)
 daspect([1 1 5])
 xlabel('x');
 ylabel('y');
 
figure(5)
title(' Num results - perm');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),(1+tanh(p_res/eps_perm))/2.)
 daspect([1 1 1])
  
%mat2cell(kr)

%%
nn_nodes=find(p_res<0.04);
 figure(1);
plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
plot(mesh.nodes(nn_nodes,1),mesh.nodes(nn_nodes,2),'or');
 
%% FLUXES at centroid 
phi_res=h_res;
Qelt=zeros(length(mesh.connectivity),2);
centroids =zeros(length(mesh.connectivity),2);
% loop over all the element
for e=1:length(mesh.connectivity)
    
     % get nodes of element e
     n_e = mesh.connectivity(e,:);
     %corresponds coordinates
     coor = mesh.nodes(n_e,:);
     centroids(e,:)=mean(coor);
     % creat Element object  
     
     % here we would have to do a switch for Tri6 mesh
     local_elt=ElementTri3(coor,'2D');

     % corresponding local solution - scalar dof mapping
     usol = phi_res(n_e);
     [qel]=GetElementFlux(local_elt,1,usol);
     Qelt(e,:)=qel';

end 

% scatter3(centroids(:,1),centroids(:,2),Qelt(:,1),...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75])
% 
    scatter3(centroids(:,1),centroids(:,2),sqrt(Qelt(:,2).^2+Qelt(:,1).^2),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])

%% PROJECTION FOR FLUX

Q =ProjectFlux(mesh,'2D',k_current,h_res);

figure(6)
title(' Proj Q_x - linear elt');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1))

figure(7)
title(' Proj Q_y - linear elt');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,2))
 