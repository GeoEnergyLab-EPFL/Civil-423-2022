function [q]=ProjectFlux(mesh,simultype,cond,solution)
% Function to project the flux (derivatives of the solution) at the nodes of the mesh from the knowledge
% of the solution at the nodes 
%   q== - Cond Grad solution
% inputs:
% mesh :: mesh structure
% simultype: '2D' or 'Axis'
% Cond : scalr or matrix for conductivity
% solution : vector containing the solution at the nodes of the whole mesh
%
% outputs:
%   q :: vector of length n_d(=2 here) \times n_nodes

% we solve the global projection problem via 2 sub-problmes for each
% component of the flux, 

%  Union_elt (\int N^T N   ) \times Q_x =
%  Union_elt \int (-N^T (Cond.Grad N \times solution)_x ) 
% and same for y

% Step 1 : creation of the mesh matrix 
[M]= AssembleMassMatrix(mesh,1,simultype); % rho=1

% Step 2 : Assembly of the element force vectors for the 2 compo
f_x = zeros(length(mesh.nodes),1);
f_y=f_x;

% check 
n_mat=length(unique(mesh.id)); % number of mat in the mesh
 
[~, nc]=size(mesh.connectivity);
switch nc
    case 3
        eltype='Tri3';
    case 6
        eltype='Tri6';
end

if length(cond(:,1)) ~= n_mat
        error('number of mat and size of cond list not equal ');       
end

% loop over all the element
for e=1:length(mesh.connectivity)
    
     % get nodes of element e
     n_e = mesh.connectivity(e,:);
     %corresponds coordinates
     coor = mesh.nodes(n_e,:);
     
     m_id = mesh.id(e);
     cond_e=cond(m_id);
      
     % creat Element object  
     
     % --- Switch in function of element type
     switch eltype
         case 'Tri3'
             local_elt=ElementTri3(coor,simultype);
         case 'Tri6'
             local_elt=ElementTri6(coor,simultype);
     end

     % corresponding local solution - scalar dof mapping
     usol = solution(n_e);
     
    [fel_x,fel_y]=ProjectElementFlux(local_elt,cond_e,usol);
     
    f_x(n_e) = f_x(n_e)+fel_x ;
    f_y(n_e) = f_y(n_e)+fel_y ;
    
end
% chech nan for unconfined flow
[k]=find(isnan(f_x)==1);
f_x(k)=0.;
[k]=find(isnan(f_y)==1);
f_y(k)=0.;


% Step 3 Solution of the projection problem

% 3.a for x

q_x = M\f_x;

% 3.b for y;

q_y = M\f_y;


q = [ q_x q_y];


end