function [newNodes , newConnectivity] = Tri3ToTri6(meshTri3)
    % Tri3ToTri6 is a function changing a Tri3 element into a Tri6 element by
    % addition of the lateral mid point.
    % --- Initiating necessary elements
    
    n_el = length(meshTri3.connectivity(:,1)); % number of elements
    n_nodes = length(meshTri3.nodes(:,1)); % number of nodes
    
    nodesSorted = meshTri3.nodes(meshTri3.connectivity(:),:);
    
    nodesToAdd = [(nodesSorted(1:n_el,:)+nodesSorted(n_el+1:2*n_el,:))/2;
    (nodesSorted(n_el+1:2*n_el,:)+nodesSorted(2*n_el+1:3*n_el,:))./2;
    (nodesSorted(1:n_el,:)+nodesSorted(2*n_el+1:3*n_el,:))./2];

    [nodesToAdd,~,ic] = unique(nodesToAdd,'rows','stable');

    newNodes = [meshTri3.nodes; nodesToAdd];
    newConnectivity = [meshTri3.connectivity, n_nodes+ic(1:n_el),...
        n_nodes+ic(n_el+1:2*n_el), n_nodes+ic(2*n_el+1:end)];

end