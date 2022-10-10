% Exercice 3.1 -- unconfined flow in an earth dam
% 
% The goal of this exercise is to understand how to calculate an unconfined
% flow in an earth dam using a relative permeability.
% 
% For the solution we use a fixed point iterations scheme with 
% under-relaxation
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

% Call mesh-generator MESH2D
opts.kind = 'delfront';

h_x = 1.; %% Max elt area -> control the refinement of the mesh here
[mesh.nodes,mesh.edge, mesh.connectivity,mesh.id] = ...
    refine2(node,edge,[],opts,h_x) ;

% plotting mesh
figure(1);
plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
ne_t=length(mesh.connectivity(:,1));

%  Boundary Conditions
% left  :: y= x/m
% right edge ::  y= (x_ext-x)/m1
% bottom edge :: y=0 x
% top edge  :: y=hd x

% find the nodes belonging to each boundary
left_edge=find(abs(mesh.nodes(:,1)/m-mesh.nodes(:,2))<0.001);
right_edge=find(abs((x_ext-mesh.nodes(:,1))/m1-mesh.nodes(:,2))<0.001);
bottom_edge=find(mesh.nodes(:,2)==0);
top_edge=find(mesh.nodes(:,2)==hd);

% Nodes on the left side above the water level (unsaturated;b.c. zero flow)
[il]=find(mesh.nodes(left_edge,2)>=hw);
left_edge_unsat=left_edge(il);

% Nodes on the left side below the water level (saturated; hydrostatic)
[il]=find(mesh.nodes(left_edge,2)<hw);
left_edge_sat=left_edge(il);

% We now combine all nodes where we are either sure or do not know yet if
% they are unsaturated. This is the combination of the nodes above water
% level on the left side (unsaturated), the ones on top (unsaturated), and
% the nodes on the right/downstream side. The nodes on the right side can
% be saturated or not (depending on where the water level passes in the
% earth dam).
[nodes_unsat, ia_unsat, ic]=unique([left_edge_unsat;top_edge;right_edge;]);

% plotting boundary nodes
plot(mesh.nodes(left_edge_sat,1),mesh.nodes(left_edge_sat,2),'ob');
hold on
plot(mesh.nodes(nodes_unsat,1),mesh.nodes(nodes_unsat,2),'or');
hold on

%% Here you have to impose the boundary conditions of the problem
 
% Nodes on the left side above the water level (unsaturated;b.c. zero flow)
% left unsaturated -> p=0 -> h=y (we are sure this is unsaturated)
h_left_unsat=mesh.nodes(left_edge_unsat,2);

% Nodes on the left/upstream side below the water level (saturated; 
% hydrostatic)
% hydrostatic -> p=gamma_w (hw-y) -> h= (hw-y) + y == hw
h_left_sat=hw*ones(length(left_edge_sat),1);

% Now the Dirichlet boundary condition at the downstream face (right side)
% right p = 0-> h=y (we estimate zero pressure everywhere there)
h_right=mesh.nodes(right_edge,2);

% And now the boundary condition at the crest of the dam (top)
% top -> p=0 -> h=y (we are sure this is unsaturated)
h_top=mesh.nodes(top_edge,2);

% Combining all b.c. on unsaturated parts
h_all_unsat=[h_left_unsat;h_top;h_right];
h_unsat=h_all_unsat(ia_unsat);  % we use here ia_unsat to re-order properly
                                % such that we have correspondance between
                                % nodes_unsat and h_unsat


% Deleting duplication of fixed nodes (just in case).
[nodes_fixed, ia, ic]=unique([left_edge_sat; ]) ;
h_set=h_left_sat; % The known hydraulic heads are said to be the ones below
                  % water level on the upstream side.

% Fixed Point Interation Scheme preparation

k_o=[1.]; % Hydraulic conductivity of material, [L/s]
mesh.id=[1:length(mesh.connectivity(:,1))]';
% current permeability
k_current=k_o+0.*[1:length(mesh.connectivity(:,1))]';
% Initialization of hydraulic conductivity vector
h_res=0.+0*[1:length(mesh.nodes)]'; % Initialization of piezometric head 
                                    % solution

h_res(nodes_fixed)=h_set;  % enforce BC on the upstraem side below the 
                           % the water level.

gammaw=10; % Weight of the water
beta_fix=0.9; % Under relaxation coefficient (1 if none)
eps_perm=0.01; % Parameter that controls sharpness of relative permeability
               % function
tolerance = 2.e-5; % Tolerance for iteration
err=1; % Initialization of error
iter_max=20; % Max. number of iterations
k=0; % Iteration count

while (k<iter_max) && (err>tolerance)
    k=k+1;
     % start the under-relaxation after first iteration only as the first
     % iterate is to get an full initial guess
    if k==1
        beta=1;
    else
        beta = beta_fix;
    end
    
    % We assemble the conductivity matrix with the current permeability
    [C] = AssembleConductivityMatrix(mesh,k_current,'2D');
    
    % First, we have to separate the nodes where the piezometric head is
    % unknown from nodes where the piezometric head is fixed (because of 
    % the boundary conditions).   
    eq_to_solve=setdiff([1:length(mesh.nodes)],nodes_fixed)';
        
    % Compute the force vector
    f = -C(eq_to_solve,nodes_fixed)*h_set ;
    
    % Compute the unknown piezometric heads at each step
    h_aux = C(eq_to_solve,eq_to_solve)\f;
    
    % Solution of the previous step
    h_res_1=h_res;
    
    % glue back fixed and solved nodes to get the head after the current
    % iteration
    h_res(nodes_fixed)=h_set;  % not really needed as outsie the loop (just
                               % to be sure)
    h_res(eq_to_solve)=h_aux;  % Solution of the yet unknown piezometric 
                               % heads
    
    % Under-Relax the solution of the current step
    h_res=(1-beta)*h_res_1+beta*h_res;
    
    % Compute pressure from head
    p_res = gammaw*(h_res - mesh.nodes(:,2));
    
    % check constraints on pressure on the unsaturated boundary
    % (if really unsaturated the pressure should be 0)
    % We perform the check via the value of the hydraulic head
    % (which should thus be below the y coordinate)
    % set of nodes on the unsat boundary where pressure is >0 
    % will need to be fixed to p=0
    ic = find(h_res(nodes_unsat)>mesh.nodes(nodes_unsat,2)); 
    
    % check if the flux are < 0 on seeping part of the unsat boundary 
    ff=C*h_res; 
    % get set of nodes where p > 0 and qn <0
    if_o = find(ff(nodes_unsat(ic))<0);
    i_fix_u=ic(if_o); % These are the indices of the nodes with a pressure
                      % above zero and a flux.
                      % We will fix there the pressure to be zero.
    
    % We adapt the fixed nodes by adding the ones just found.                  
    nodes_fixed=[left_edge_sat; nodes_unsat(i_fix_u)];
    
    % We fix the hydraulic heads at the locations of the previously found 
    % points.
    h_set=[h_left_sat;h_unsat(i_fix_u)];

    % Compute the new relative permeability of each element
    [k_current,pm]=ComputeRelPerm(mesh,p_res,eps_perm);
    
    % We check if we got a valid new permeability.
    [jk]=find(isnan(k_current));
    if ~isempty(jk)
        disp('error - Nan Perm');
    end
    
    % Compute the error of the current step by comparing the solution of 
    % the hydraulic head to the one of the previous step.
    err=median(abs(h_res-h_res_1)./h_res);   
    disp(['median of relative change at iteration ', num2str(k),...
        ' is : ', num2str(err)]);
    
end


%% plotting solution on the grid.

figure(3)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),h_res)
daspect([1 1 1])
title(' Num results - head');
xlabel('x');
ylabel('y');

figure(4)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),p_res)
daspect([1 1 5])
title(' Num results - p');
xlabel('x');
ylabel('y');

figure(5)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),(1+tanh(p_res/eps_perm))/2.)
daspect([1 1 .1])
title(' Num results - perm');

%% Visualize buondaries
% You can visualize here by commenting the different parts of the 
% unsaturated boundary
ff=C*h_res;
[fi]=find(ff(nodes_unsat)>0);
figure(1);
plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
hold on
plot(mesh.nodes(nodes_unsat,1),mesh.nodes(nodes_unsat,2),'og');
plot(mesh.nodes(nodes_unsat(ic),1),mesh.nodes(nodes_unsat(ic),2),'or');
plot(mesh.nodes(nodes_unsat(i_fix_u),1),...
    mesh.nodes(nodes_unsat(i_fix_u),2),'ob');
plot(mesh.nodes(nodes_unsat(fi),1),mesh.nodes(nodes_unsat(fi),2),'oy');
plot(mesh.nodes(left_edge_sat,1),mesh.nodes(left_edge_sat,2),'*k');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2)...
    ,(1+tanh(p_res/eps_perm))/2.)


%% Calculate the fluxes at the centroid
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

figure(2);
scatter3(centroids(:,1),centroids(:,2),Qelt(:,1),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
title('   q_x (centroids) - ');

figure(3);
scatter3(centroids(:,1),centroids(:,2),sqrt(Qelt(:,2).^2+Qelt(:,1).^2),...
   'MarkerEdgeColor','k',...
   'MarkerFaceColor',[0 .75 .75])
title('   q_y (centroids) - ');

%% PROJECTION FOR FLUX

Q =ProjectFlux(mesh,'2D',k_current,h_res);

figure(4)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1))
daspect([1 1 .01])
title('   q_x (projection) - ');

figure(5)
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,2))
daspect([1 1 .01])
title(' q_y (projection) -');
