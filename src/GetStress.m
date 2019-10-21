function [Selt] = GetStress(mesh,params,u,simultype)
% Function to get the stress (derivatives of the solution) at the centroid
% from the knowledge of the solution at the nodes
%   sigma == D B u
% inputs:
% mesh :: mesh structure
% simultype :: '2D' or 'Axis'
% u :: displacement field at the nodes
% params :: young's modulus and poissons ratio
%
% output:
%   Selt :: vector of length n_d(=3 ('2D') or 4 ('Axis')) \times n_elements

%-- Asses element type
[~, nc]=size(mesh.connectivity);
switch nc
    case 3
        eltype='Tri3';
    case 6
        eltype='Tri6';
end  

%-- Transform young's modulus and poisson's ratio into bulk and shear
%modulus
k = params(1)/(3*(1-2*params(2)));
g = params(1)/(2*(1+params(2)));

%-- calculate constitutive matrix and prepare output
if strcmp(simultype,'2D')
    D = Elastic_Isotropic_Stiffness(k,g,'PlaneStrain');
    Selt=zeros(length(mesh.connectivity(:,1)),3);
elseif strcmp(simultype,'Axis')
    D = Elastic_Isotropic_Stiffness(k,g,'Axisymmetry');
    Selt=zeros(length(mesh.connectivity(:,1)),4);
else
    error('Approximation not yet implemented')
end

% loop over all the element
for e=1:length(mesh.connectivity(:,1))
    
     % get nodes of element e
     n_e = mesh.connectivity(e,:);
     % get corresponding DOF
     n_dof = reshape([2.*n_e-1; 2*n_e],[1,nc*2]);
     % get corresponds coordinates
     coor = mesh.nodes(n_e,:);
     
     % creat Element object  
     % --- Switch on element type
     switch eltype
         case 'Tri3'
             local_elt=ElementTri3(coor,simultype);
         case 'Tri6'
             local_elt=ElementTri6(coor,simultype);
     end
     
     % corresponding local solution - scalar dof mapping
     sol = u(n_dof);
     [Sel]=GetElementStress(local_elt,D,sol);
     Selt(e,:)=Sel;

end

end

