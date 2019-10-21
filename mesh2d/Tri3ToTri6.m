function [newNodes , newConnectivity] = Tri3ToTri6(meshTri3)
% Tri3ToTri6 is a function changing a Tri3 element into a Tri6 element by
% addition of the lateral mid point.
    % todo: this needs to be vectorized -> loop are super slow in matlab
    % --- Initiating necessary elements
    n_el = length(meshTri3.connectivity(:,1));
    newConnectivity = zeros(n_el,6);
    intermediateNodes = zeros(6,2);
    newNodes = [];
    % --- Subsequent iterations
    for i = 1:n_el
        intermediateNodes(1:3,:)= meshTri3.nodes(...
            meshTri3.connectivity(i,:),:); % nodes of Tri3 Element
        intermediateNodes(6,:) = (meshTri3.nodes(...
            meshTri3.connectivity(i,1),:)+meshTri3.nodes(...
            meshTri3.connectivity(i,3),:))./2; % node between node 1 and 3 
        % of the element
        for j = 1:2 % Loop to get remaining mid line nodes
            intermediateNodes(j+3,:) = (meshTri3.nodes(...
            meshTri3.connectivity(i,j),:)+meshTri3.nodes(...
            meshTri3.connectivity(i,j+1),:))./2;
        end
        newNodes = unique([newNodes;intermediateNodes],'rows','stable');
        % Get new nodes by unique to avoid double storage
        for j = 1:6 % Loop to set up new connectivity
            newConnectivity(i,j) = find(newNodes(:,1) == ...
                intermediateNodes(j,1) & newNodes(:,2) == ...
                intermediateNodes(j,2));
        end
    end
end

