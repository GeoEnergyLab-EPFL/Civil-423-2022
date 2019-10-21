function [qel]=GetElementFlux(eltObj,Cond,usol)
 % Element level function
 % Compute the flux at integration points from the solution  
 %     - Cond. \grad usol
 % 
 % inouts:
 % eltObj is an element object
 % Cond: is the related conductivity (scalar or matrix if tensor) 
 % usol : is the solution (of the scalar laplacian) at the 3 nodes of the triangle
 % outputs:
 % qel : element flux vector
 
   if isa(eltObj,'ElementTri3')
        % Linear Triangle
        
        %  Single Gauss point (Constant gradient for linear triangle)
        xeta=[1./3 1./3];
        [DNaDx,j]=GradN(xeta,eltObj);
        
        qel=-Cond*DNaDx*usol;
        
   elseif isa(eltObj,'ElementTri6')
        % Quadratic Triangle
        
        %  three Gauss points (linear gradient for cubic triangle)
        xeta=[1./6 1./6;
            2./3 1./6;
            1./6 2./3];
        
        qel = zeros(length(xeta(:,1)),2);
        
        for i = 1:length(xeta(:,1))
            [DNaDx]=GradN(xeta(i,:),eltObj);
            qel(i,:)=-Cond*DNaDx*usol;
        end
        qel = mean(qel);
   else
        error('Element not yet implemented');
    end
    
end
