% Exercice 1.1  - Computational Geomechanics 2019-20
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

% this is equivalent to the following few lines involving a for loop ;(
 
%  h_true=0*[1:length(mesh.nodes)]';
%for i=1:length(h_true)
%      xx=mesh.nodes(i,:); % coor of nodes i
%     h_true(i)= sin(pi*xx(1))* sinh(pi*xx(2))/sinh(pi);
% end

% let's plot it on this unstructured mesh-> function trisurf
figure(2)
title(' Exact head at the nodes');
trisurf(mesh.connectivity,mesh.nodes(:,1),mesh.nodes(:,2),h_true)
 

%% PROJECTION FOR FLUX

% Here you need to finish coding up the function ProjectFlux !!!
 Q =ProjectFlux(mesh,'2D',1.,h_true);
 
 % plot numerical and exact solution for flux components at the nodes
 
 % compute abs and relative error
 
 
 
 