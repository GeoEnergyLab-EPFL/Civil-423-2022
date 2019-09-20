function [C] = AssembleConductivityMatrix(mesh,cond,simultype)
%AssembleConductivityMatrix :: assemble the conductivity (aka Laplacian
%type) matrix for a mesh structure
%   mesh structure with the following fields (node, connectivity, edge, id)
%   cond = (value of 'conductivity' / perm) - scalar (can be extended to
%   list)
%   simultype ='2D' or 'Axis'
%   eltype = 'Tri3' or 'Tri6' CAN BE access from size of the connectivity
%   array of the mesh ;)
%   Detailed explanation goes here

 n_u = length(mesh.nodes);
 
 C=sparse(n_u,n_u);

 % check 
 n_mat=length(unique(mesh.id)); % number of mat in the mesh
 

if length(cond(:,1)) ~= n_mat
        error('number of mat and size of cond list not equal ');       
end
 
[nr nc]=size(mesh.connectivity);
switch nc
    case 3
        eltype='Tri3';
    case 6
        eltype='Tri6';
end
    
 
for e=1:length(mesh.connectivity)
       
     % get nodes of element e
     
     n_e = mesh.connectivity(e,:);
     %corresponds coordinates
     coor = mesh.nodes(n_e,:);
     
     m_id = mesh.id(e);
     
     % create Element object  
     % here would have to do a switch for Tri6 mesh
     local_elt=ElementTri3(coor,simultype);
     % creat Element object
     % --- Switch in function of element type
     switch eltype
         case 'Tri3'
             local_elt=ElementTri3(coor,simultype);
         case 'Tri6'
             local_elt=ElementTri6(coor,simultype);
     end
     
     % get element conductivity matrix
     cond_e=cond(m_id);
     Cel=ElementConductivityMatrix(local_elt,cond_e);

    % 
    C(n_e,n_e)=C(n_e,n_e)+Cel;
    
 end
 

end

