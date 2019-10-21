function [K] = AssembleStiffnessMatrix(mesh,params,simultype)
%AssembleStiffnessMatrix :: assemble the stiffness 
%matrix for a mesh structure
%   mesh structure with the following fields (node, connectivity, edge, id)
%   params = matrix with youngs modulus and poissons ratio
%   simultype ='2D' or 'Axis'
%   Detailed explanation goes here

 n_u = length(mesh.nodes); % total number of nodes
 
 K=sparse(2*n_u,2*n_u); % prepae the stiffness matrix

[~, nc]=size(mesh.connectivity); % switch on element type
switch nc
    case 3
        eltype='Tri3';
    case 6
        eltype='Tri6';
end

% transforming youngs modulus and poisson's ratio into bulk and shear
% modulus
k = params(1)/(3*(1-2*params(2)));
g = params(1)/(2*(1+params(2)));

% create the constitutive matrix D
if strcmp(simultype,'2D')
    D = Elastic_Isotropic_Stiffness(k,g,'PlaneStrain');
elseif strcmp(simultype,'Axis')
    D = Elastic_Isotropic_Stiffness(k,g,'Axisymmetry');
else
    error('Approximation not yet implemented')
end

%------ Loop over the elements to fil the matrix
for e=1:length(mesh.connectivity(:,1))
       
     % get nodes of element e
     n_e = mesh.connectivity(e,:);
     % get the DOF numbers
     n_dof = reshape([2.*n_e-1; 2*n_e],[1,nc*2]);
     %corresponds coordinates
     coor = mesh.nodes(n_e,:);
     % creat Element object  
     % --- Switch in function of element type
     switch eltype
         case 'Tri3'
             local_elt=ElementTri3(coor,simultype);
         case 'Tri6'
             local_elt=ElementTri6(coor,simultype);
     end
     
     % get element stiffness matrix
     Kel=ElementStiffnessMatrix(local_elt,D);

    % Assemble the matrix
    K(n_dof,n_dof)=K(n_dof,n_dof)+Kel;
    
end
end

