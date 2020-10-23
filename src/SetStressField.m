function S = SetStressField(varargin)

%---------------------------------------------- extract args
    body_force = cell(1,2);
    
    if (nargin<3), error('SetStressField: to few input arguments'); end
    if (nargin>=+1), mesh = varargin{1}; end
    if (nargin>=+2), simultype = varargin{2}; end
    if (nargin>=+3), body_force = varargin{3}; end
    
%---------------------------------------------- check parameter input         
    
    % prepare the output vector of stresses
    n_u = length(mesh.nodes);
 
    S=zeros(2*n_u,1); 
    
    % decide on element type
    [~, nc]=size(mesh.connectivity);
    switch nc
        case 3
            eltype='Tri3';
        case 6
            eltype='Tri6';
    end
    
    
    if ~isempty(body_force{1}) 
        il = ismember(mesh.connectivity,body_force{1});
        % switch on mesh here to get elements with applied stresses
        if length(mesh.connectivity(1,:))==3
            elt_line=find(sum(il,2)==3);
        elseif length(mesh.connectivity(1,:))==6
            elt_line=find(sum(il,2)==6);
        else
            error(" Error " ); 
        end
        S_field = body_force{2};
    end
    
    %---------------------------------------------- loop
        
    for i = 1:length(elt_line)
         e = elt_line(i);
         % get nodes of element e
         n_e = mesh.connectivity(e,:);
         n_dof = reshape([2.*n_e-1; 2*n_e],[1,nc*2]);
         %corresponds coordinates
         coor = mesh.nodes(n_e,:);
         
         % get initial fiel
         S_el_set = S_field;

         % create Element object  
         % --- Switch in function of element type
         switch eltype
             case 'Tri3'
                 local_elt=ElementTri3(coor,simultype);
             case 'Tri6'
                 local_elt=ElementTri6(coor,simultype);
         end

         % get element body force vector
         Sel = SetElementStressField(local_elt,S_el_set);

        % 
        S(n_dof)=S(n_dof)+Sel;

     end
    
end

