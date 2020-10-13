% Exercice 1.1  - Computational Geomechanics  
% 
% The goal of this exercise is to understand the so-called "projection"
% procedure which consist in obtaining/projecting the flux at the nodes
% from the knowledge of the solution at the nodes
% We assume Darcy's law : q_i = - k h_,i
%
% 
% On the unit square, let's assume that the head h has the following form :
%     h(x,y) = sin ( pi x) sinh(pi y) / sinh(pi) 
%  (actually it is solution of the Laplacian)
  
%% mesh the unit square with triangle 

% define the unit square 
% the coordinates of the 4 nodes
node_coor = [   0  0 ; 0 1 ; 1 1; 1 0];

% the connectivity of the segment defining the closed boundary
edge = [];
for e=1:3
    edge = [ edge ; e e+1 ];
end
edge = [ edge ; 4  1 ];
 
% call mesh-generator MESH2D (included)
% opts.kind = 'delfront';

 h_x = 0.05; %% Max elt area -> control the refinement of the mesh here
 [mesh.nodes,mesh.edge, mesh.connectivity,mesh.id] = refine2(node_coor,edge,[],[],h_x) ; 

 figure(1);
 plotmesh(mesh.nodes,mesh.connectivity,[.2 .2 .2],'w')
 ne_t=length(mesh.connectivity(:,1)) 
 
 %% compute the head at the location of all nodes
 % in matlab it's good practice to avoid loops and use vectorized code
 xx=mesh.nodes(:,:); % matrix of all nodes
 h_true=sin(pi*xx(:,1)).*sinh(pi*xx(:,2))/sinh(pi);

% let's plot it on this unstructured mesh-> function trisurf
figure(2)
title(' Exact head at the nodes');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),h_true)
 

%% PROJECTION FOR FLUX

% call to ProjectFlux 
Q =ProjectFlux(mesh,'2D',1.,h_true);
 
 % for comparison, estimate the exact flux at the nodes (in a similar way
 % than for the head above
 % analytical 


q_x_true= -pi*cos(pi*xx(:,1)).*sinh(pi*xx(:,2))/sinh(pi);
q_y_true= -pi*sin(pi*xx(:,1)).*cosh(pi*xx(:,2))/sinh(pi);

% absolute error
abs_q_x = abs((Q(:,1)-q_x_true) ) ; 
abs_q_y = abs((Q(:,2)-q_y_true) ) ; 

% relative error
% in order to avoid blow if the true solution is zero, at this points, the
% rel error = abs error
[kx]=find(abs(q_x_true)<1.e-8); %cutting off at 1.e-8
rel_q_x= abs_q_x./q_x_true;  %  numerical blow up due to 1/0
rel_q_x(kx)=abs_q_x(kx);

[ky]=find(abs(q_y_true)<1.e-8);
rel_q_y= abs_q_y./q_y_true;
rel_q_y(ky)=abs_q_y(ky);


% Plots 

figure(3)
title(' Proj Q_x - linear elt - numerics vs true');
subplot(2,1,1), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,1))
subplot(2,1,2), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),q_x_true(:,1)) 


figure(4)
title(' Proj Q_y - linear elt - numerics vs true');
subplot(2,1,1),trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),Q(:,2))
subplot(2,1,2), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),q_y_true(:)) 
 

figure(5)
title(' Absolute error on flux (q_x & q_y )');

subplot(2,1,1), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),abs_q_x)
subplot(2,1,2), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),abs_q_y) 


figure(6)
title(' relative error on flux  (q_x & q_y )');

subplot(2,1,1), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),rel_q_x)
subplot(2,1,2), trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),rel_q_y) 
