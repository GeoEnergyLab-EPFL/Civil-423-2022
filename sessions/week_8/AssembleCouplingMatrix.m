function [Ce] = AssembleCouplingMatrix_exo(meshE,meshP,alpha,simultype)
%AssembleCouplingMatrix :: assemble the coupling matrix for poroelastic
% problems for a mesh structure.
%   mesh structure with the following fields (node, connectivity, edge, id)
%   alpha = value of 
%   simultype ='2D' or 'Axis'
%   Detailed explanation goes here

 n_E = length(meshE.nodes(:,1));
 n_p = length(meshP.nodes(:,1));
 
%----- This is a rectangular matrix (relating displacements [2*n_nodes DOF]
%      to the pore pressure [n_nodes DOF])
 Ce=sparse(n_E*2,n_p); 

[~, ncP]=size(meshP.connectivity);
[~, ncE]=size(meshE.connectivity);
switch ncP
    case 3
        eltypeP='Tri3';
    case 6
        eltypeP='Tri6';
end

switch ncE
    case 3
        eltypeE='Tri3';
    case 6
        eltypeE='Tri6';
end
%----- Loop over the elements
for e=1:length(meshP.connectivity(:,1))

     % get nodes of element e     
     n_e_E = meshE.connectivity(e,:);
     n_e_P = meshP.connectivity(e,:);
     % get number of DOF's
     n_dof_e = reshape([2.*n_e_E-1; 2*n_e_E],[1,ncE*2]);
     n_dof_p = n_e_P;
     %corresponds coordinates
     coorE = meshE.nodes(n_e_E,:);
     coorP = meshP.nodes(n_dof_p,:);
     % creat Element object  

     % --- Switch in function of element type
     switch eltypeE
         case 'Tri3'
             local_elt_E=ElementTri3(coorE,simultype);
         case 'Tri6'
             local_elt_E=ElementTri6(coorE,simultype);
     end
     
     switch eltypeP
         case 'Tri3'
             local_elt_P=ElementTri3(coorP,simultype);
         case 'Tri6'
             local_elt_P=ElementTri6(coorP,simultype);
     end

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     % get element conductivity matrix
     Ceel=ElementCouplingMatrix_exo(local_elt_E,local_elt_P,alpha);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
    % 
    Ce(n_dof_e,n_dof_p)=Ce(n_dof_e,n_dof_p)+Ceel;
    
 end
 

end


