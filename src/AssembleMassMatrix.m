function [M] = AssembleMassMatrix(mesh,rho,simultype)
%AssembleConductivityMatrix :: assemble the conductivity (aka Laplacian
%type) matrix for a mesh structure
%   mesh structure with the following fields (node, connectivity, edge, id)
%   rgo = (value of 'density' / storage parameter) - scalar (can be extended to
%   list)
%   simultype ='2D' or 'Axis'
%   Detailed explanation goes here

 n_u = length(mesh.nodes);
 
 M=sparse(n_u,n_u); 

[~, nc]=size(mesh.connectivity);
switch nc
    case 3
        eltype='Tri3';
    case 6
        eltype='Tri6';
end
   

for e=1:length(mesh.connectivity(:,1))
       
     % get nodes of element e
     
     n_e = mesh.connectivity(e,:);
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
     
     % get element conductivity matrix
     Mel=ElementMassMatrix(local_elt,rho);

    % 
    M(n_e,n_e)=M(n_e,n_e)+Mel;
    
 end
 

end

