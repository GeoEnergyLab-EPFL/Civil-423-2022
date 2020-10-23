function [S]=ProjectStress(mesh,simultype,params,solution)
% Function to project the stresses (derivatives of the solution) at the 
% nodes of the mesh from the knowledge of the solution at the nodes 
%   sigma ==  D B solution
% inputs:
% mesh :: mesh structure
% simultype :: '2D' or 'Axis'
% solution :: displacement field at the nodes
% params :: young's modulus and poissons ratio
%
% output:
%   S :: vector of length n_d(=3 ('2D') or 4 ('Axis')) \times n_nodes

% we solve the global projection problem via 3 ('2D') respectively 4
% ('Axis') sub-problmes for each component of the stress vector. 

%  Union_elt (\int N^T N   ) \times S_x =
%  Union_elt \int (N^T (D b \times solution)_x ) 
%  and same for the other stresses

%----- Step 1 : creation of the mesh and constitutive matrix 
[M]= AssembleMassMatrix(mesh,1,simultype); % rho=1

%Transform young's modulus and poisson's ratio into bulk and shear modulus
k = params(1)/(3*(1-2*params(2)));
g = params(1)/(2*(1+params(2)));

%-- calculate constitutive matrix and prepare output
switch simultype
    case '2D'
        D = Elastic_Isotropic_Stiffness(k,g,'PlaneStrain');
        %----- Step 2 : Assembly of the element force vectors
        f_x = zeros(length(mesh.nodes),1);
        f_y=f_x;
        tau = f_x;
    case 'Axis'
        D = Elastic_Isotropic_Stiffness(k,g,'Axisymmetry');
        %----- Step 2 : Assembly of the element force vectors
        f_r = zeros(length(mesh.nodes),1);
        f_z=f_r;
        f_theta = f_r;
        tau = f_r;
    otherwise
        error('Approximation not yet implemented')
end
 
n_mat=length(unique(mesh.id)); % check number of mat in the mesh
 
%-- Asses element type
[~, nc]=size(mesh.connectivity);
switch nc
    case 3
        eltype='Tri3';
    case 6
        eltype='Tri6';
end

if length(params(:,1)) ~= n_mat
        error('number of mat and size of cond list not equal ');       
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
     
     % --- Switch in function of element type
     switch eltype
         case 'Tri3'
             local_elt=ElementTri3(coor,simultype);
         case 'Tri6'
             local_elt=ElementTri6(coor,simultype);
     end

     % corresponding local solution - scalar dof mapping
     usol = solution(n_dof);
    
     switch simultype
        case '2D' 
            [fel_x,fel_y,tauel,~]=ProjectElementStress(local_elt,D,usol);

            f_x(n_e) = f_x(n_e)+fel_x ;
            f_y(n_e) = f_y(n_e)+fel_y ;
            tau(n_e) = tau(n_e)+tauel ;
        case 'Axis'
            [fel_r,fel_z,tauel,fel_theta]=...
                ProjectElementStress(local_elt,D,usol);

            f_r(n_e) = f_r(n_e)+fel_r ;
            f_z(n_e) = f_z(n_e)+fel_z ;
            f_theta(n_e) = f_theta(n_e)+fel_theta ;
            tau(n_e) = tau(n_e)+tauel ;
    end
    
end

% Step 3 Solution of the projection problem

switch simultype
    case '2D' 
        S_x = M\f_x;
        S_y = M\f_y;
        S_tau = M\tau;
        
        S = [ S_x S_y S_tau];
    case 'Axis'
        S_r = M\f_r;
        S_z = M\f_z;
        S_tau = M\tau;
        S_theta = M\f_theta;
        
        S = [ S_r S_z S_tau S_theta];
end

end