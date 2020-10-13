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

h_x = 2; %% Max elt area -> control the refinement of the mesh here
[mesh.nodes,mesh.edge, mesh.connectivity,mesh.id] = refine2(node,edge,[],opts,h_x) ;

% plotting mesh
figure(1);
plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
ne_t=length(mesh.connectivity(:,1));

%  Boundary Conditions
% left  :: y= x/m
% right edge ::  y= (x_ext-x)/m1
% bottom edge :: y=0 x
% top edge  :: y=hd x

% find the nodes belonging to each boundaries
left_edge=find(abs(mesh.nodes(:,1)/m-mesh.nodes(:,2))<0.001);
right_edge=find(abs((x_ext-mesh.nodes(:,1))/m1-mesh.nodes(:,2))<0.001);
bottom_edge=find(mesh.nodes(:,2)==0);
top_edge=find(mesh.nodes(:,2)==hd);

% list of mesh on the possinle unsaturated boundary

[il]=find(mesh.nodes(left_edge,2)>=hw);
left_edge_unsat=left_edge(il);

% list nodes on the upstram side of the dam
[il]=find(mesh.nodes(left_edge,2)<hw);
left_edge_sat=left_edge(il);

% unsat nodes
[nodes_unsat, ia_unsat, ic]=unique([left_edge_unsat;top_edge;right_edge;]);

% plotting boundary nodes
plot(mesh.nodes(left_edge_unsat,1),mesh.nodes(left_edge_unsat,2),'ob');
hold on
plot(mesh.nodes(nodes_unsat,1),mesh.nodes(nodes_unsat,2),'or');
hold on

% Here you have to impose the boundary conditions of the problem
 
% left unsat hydrostatic -> p=gamma_w (hw-y) -> h= (hw-y) + y == hw

h_left_unsat=------

% saturated upstrtam  face of the dam
h_left_sat=---------
%
% unsaturated face - prepare constraints if the fluid is seeping through
% this face
%
% Now the Dirichlet boundary condition at the downstream face (right side)
% right p = 0-> h=y
h_right=----------
% And now the boundary condition at the crest of the dam (top)
% top -> p=0 -> h=y
h_top=-------------
h_all_unsat=[h_left_unsat;h_top;h_right];
h_unsat=h_all_unsat(ia_unsat);  % we use here ia_unsat to re-order properly


% Deleting duplication of fixed nodes (just in case).
[nodes_fixed, ia, ic]=unique([left_edge_sat; ]) ; % be careful there is duplicate-> use unique
h_set=h_left_sat;

% Fixed Point Interation Scheme preparation

k_o=[1.]; % Hydraulic conductivity of material, [L/s]
mesh.id=[1:length(mesh.connectivity(:,1))]';
% current perm.
k_current=k_o+0.*[1:length(mesh.connectivity(:,1))]';
% Initialization of hydraulic conductivity vector
h_res=0.+0*[1:length(mesh.nodes)]'; % Initialization of piezometric head solution

h_res(nodes_fixed)=h_set;  %% enforce BC

beta=1; % Under relaxation coefficient (1 if none)
eps_perm=0.01; % Parameter that controls sharpness of relative permeability function
tolerance = 1.e-6; % Tolerance for iteration
err=1; % Initialization of error
iter_max=20; % Max. number of iterations
k=0; % Iteration count

while (k<iter_max) && (err>tolerance)
    k=k+1;
     % start the under-relaxation after first iteration only as the firt
     % iterate is to get an full initial guess
    if k==1
     beta=1;
    else
        beta = 0.9;
    end
    
    [C] = AssembleConductivityMatrix(mesh,k_current,'2D');
    
    % First, we have to separate the nodes where the piezometric head is
    % unknown from nodes where the piezometric head is fixed (because of the
    % boundary conditions).
    
    eq_to_solve=setdiff([1:length(mesh.nodes)],nodes_fixed)';
        
    % Compute the force vector
    f = ----------
    
    % Compute the unknown piezometric heads at each step
    h_aux = ------------
    
    % Solution of the previous step
    h_res_1=h_res;
    
    % glue back fixed and other nodes the solution of the current step
    h_res(nodes_fixed)=h_set;  % not really needed as outsie the loop
    h_res(eq_to_solve)=h_aux;
    
    % Under-Relax the solution of the current step by using under-relaxation
    h_res=-----------------
    % Compute pressure from head
    gammaw=10;
    p_res = gammaw*(h_res - mesh.nodes(:,2));
    
    % check constraints on pressure on -- unsat boundary
    % check via head 
    % set of nodes on the unsat boundary where pressure is >0 
    % will need to be fixed to p=0
    ic = -----------
    
    % check  flux are < 0 on seeping part of the unsat boundary 
    ff=C*h_res; 
    % get set of nodes where p> 0 and qn <0
    if_o = -----------------
    i_fix_u=ic(if_o);
    %i_fix_u=find( ((h_res(nodes_unsat)>mesh.nodes(nodes_unsat,2)) && (ff(nodes_unsat)<0))  );

    nodes_fixed=[left_edge_sat; nodes_unsat(i_fix_u)];
    h_set=[h_left_sat;h_unsat(i_fix_u)];

    % Compute new the relative permeability of each element
    [k_current,pm]=ComputeRelPerm(mesh,p_res,eps_perm);
    
    [jk]=find(isnan(k_current));
    if ~isempty(jk)
        disp('error - Nan Perm');
    end
    
    
    % Compute the error of the current step
    err=median(abs(h_res-h_res_1)./h_res);
    
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
ff=C*h_res;
[fi]=find(ff(nodes_unsat)>0);
figure(1);
plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
%plot(mesh.nodes(nodes_unsat,1),mesh.nodes(nodes_unsat,2),'og');

%plot(mesh.nodes(nodes_unsat(ic),1),mesh.nodes(nodes_unsat(ic),2),'or');

plot(mesh.nodes(nodes_unsat(i_fix_u),1),mesh.nodes(nodes_unsat(i_fix_u),2),'ob');

%plot(mesh.nodes(nodes_unsat(fi),1),mesh.nodes(nodes_unsat(fi),2),'oy');

plot(mesh.nodes(left_edge_sat,1),mesh.nodes(left_edge_sat,2),'*k');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),(1+tanh(p_res/eps_perm))/2.)


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
    % create Element object
    
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
