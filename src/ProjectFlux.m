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
%            Union_elt \int (-N^T (Cond.Grad N \times solution)_x ) 
% and same for y

% Step 1 : creation of the mesh matrix 
[M]= AssembleMassMatrix(mesh,1,simultype); % rho=1

% Step 2 : Assembly of the element force vectors for the 2 compo
f_x = zeros(length(mesh.nodes),1);
f_y=f_x;


% loop over all the element
for e=1:length(mesh.connectivity)
    
     % get nodes of element e
     n_e=-------
     %corresponds coordinates
    corr==-------
    %  mat_id  of the element....
    m_id=-------
    % get element conductivity
    cond_e=-------
    % create local Element object   (use the functino ElementTri3
     
     local_elt=ElementTri3( -------)

     % get the corresponding local solution - scalar dof mapping
     usol = solution(n_e);
     
     %  get the projected element flux
    [fel_x,fel_y]=ProjectElementFlux(local_elt,cond_e,usol);
     
    % add the contribution to the global vector
    f_x(-------) = f_x(-------)+fel_x ;
    f_y(-------) = f_y(-------)+fel_y ;
    
end
% chech NaN  (this is for unconfined flow / seepage boundary)
[k]=find(isnan(f_x)==1);
f_x(k)=0.;
[k]=find(isnan(f_y)==1);
f_y(k)=0.;


% Step 3 Solution of the projection problem

% 3.a for x

q_x = -------

% 3.b for y;

q_y = -------

% return both component
q = [ q_x q_y];


end