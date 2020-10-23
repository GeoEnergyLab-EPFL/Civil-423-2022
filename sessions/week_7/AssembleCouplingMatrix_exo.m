function [Ce] = AssembleCouplingMatrix_exo(meshE,meshP,alpha,simultype)
%AssembleCouplingMatrix :: assemble the coupling matrix for poroelastic
% problems for a mesh structure.
%   mesh structure with the following fields (node, connectivity, edge, id)
%   alpha = value of 
%   simultype ='2D' or 'Axis'
%   Detailed explanation goes here
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 n_E
 n_p

%----- This is a rectangular matrix (relating displacements [2*n_nodes DOF]
%      to the pore pressure [n_nodes DOF])
Ce

%decide for pressure element
switch 
    case 3
        eltypeP='Tri3';
    case 6
        eltypeP='Tri6';
end

% decide for elasticity element
switch 
    case 3
        eltypeE='Tri3';
    case 6
        eltypeE='Tri6';
end


%----- Loop over the elements
for e=1:length(meshP.connectivity(:,1))

     % get nodes of element e     
     n_e_E 
     n_e_P 
     % get number of DOF's
     n_dof_e  
     n_dof_p 
     %corresponds coordinates
     coorE 
     coorP
     % creat Element object  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
     
     % get element conductivity matrix
     Ceel=ElementCouplingMatrix_exo(local_elt_E,local_elt_P,alpha);
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % 
    Ce(xx,xx)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
 end
 

end

