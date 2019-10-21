function [newNodes , newConnectivity] = Tri6ToTri3(meshTri6)
% Tri3ToTri6 is a function changing a Tri3 element into a Tri6 element by
% addition of the lateral mid point.

%----------------------------------------------    

    % --- Initiating necessary elements
    n_el = length(meshTri6.connectivity(:,1));
    newConnectivity = zeros(n_el,3);
    newNodes = [];
    if isempty(solution)
        % --- Iterate on elements
        for e = 1:n_el
            newNodes = unique([newNodes;
    meshTri6.nodes(meshTri6.connectivity(e,[1 3 5]),:)],'rows','stable');

            oldConnectivity = meshTri6.connectivity(e,[1 3 5]);
            for n = 1:3 % Loop to set up new connectivity
                newConnectivity(e,n) = find(newNodes(:,1) == ...
                    meshTri6.nodes(oldConnectivity(n),1) & newNodes(:,2) == ...
                    meshTri6.nodes(oldConnectivity(n),2));
            end
        end
    else
                % --- Iterate on elements
        for e = 1:n_el
            [newNodes,~] = unique([newNodes;
    meshTri6.nodes(meshTri6.connectivity(e,[1 3 5]),:)],'rows','stable');
            
            oldConnectivity = meshTri6.connectivity(e,[1 3 5]);
            for n = 1:3 % Loop to set up new connectivity
                newConnectivity(e,n) = find(newNodes(:,1) == ...
                    meshTri6.nodes(oldConnectivity(n),1) & newNodes(:,2) == ...
                    meshTri6.nodes(oldConnectivity(n),2));
            end
        end
    end
end

