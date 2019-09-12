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
        Wl=1./2;
        [DNaDx,j]=GradN(xeta,eltObj);
        
        qel=-Cond*DNaDx*usol;
        
   elseif isa(eltObj,'ElementTri6')
        % Quadratic Triangle
        
        %  three Gauss points (linear gradient for cubic triangle)
        xeta=[1./6 1./6;
            2./3 1./6;
            1./6 2./3];
        wl=[1./3,1./3,1./3]; % Why one third?
        
        [DNaDx,j]=GradN(xeta(1,:),eltObj);
        qel=-wl(1)*Cond*DNaDx*usol;
        
        for i = 2:length(xeta(:,1))
            [DNaDx,j]=GradN(xeta(i,:),eltObj);
            qel=qel -wl(i)*Cond*DNaDx*usol;
        end
   else
        error('Element not yet implemented');
    end
    
end
