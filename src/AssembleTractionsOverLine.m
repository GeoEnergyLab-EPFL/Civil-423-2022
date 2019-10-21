function [f]=AssembleForceOverLine(mesh,nodelist,t,Type)
% t is assumed constant array of length 2
% t_n, t_s (normal and shear force -- here in 2D)
%
% i   particularize  for segment along x or z only here....


il = ismember(mesh.connectivity,nodelist);

% switch on mesh here
if length(mesh.connectivity(1,:))==3
    elt_line=find(sum(il,2)==2);
elseif length(mesh.connectivity(1,:))==6
    elt_line=find(sum(il,2)==3);
else
    error(" Error " ); 
end

f=zeros(2*length(mesh.nodes(:,1)),1);

for e=1:length(elt_line)
    
    % create 1D elt
    nn_l=il(elt_line(e),:)==1;
    
    glob_nodes=mesh.connectivity(elt_line(e),nn_l)';
    
    glob_dof = [ glob_nodes*2-1 glob_nodes*2];
    
    XY=mesh.nodes(glob_nodes,:);
    
    % this is where one need to generalize to change to coordinate s-n
    [seg_nodes,I_nodes]=sort(XY(:,1));
    switch length(XY(:,1))
        case 2
        elt_obj=ElementSeg2(seg_nodes,Type); % here I use sort... to get outward normal
        case 3
        elt_obj=ElementSeg3(seg_nodes,Type);
    end
    
        [fel_s]=ElementNeumann(elt_obj,t(1));
        [fel_n]=ElementNeumann(elt_obj,t(2));
        
        % here again here it is not at all general
        f(glob_dof(I_nodes,1))= f(glob_dof(I_nodes,1))+fel_s;
        f(glob_dof(I_nodes,2))= f(glob_dof(I_nodes,2))+fel_n;
        
    end
    
    
    
    
end
