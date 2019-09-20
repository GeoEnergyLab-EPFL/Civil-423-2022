function [kel]=ElementConductivityMatrix(eltObj,Cond)
 % Element level function
 % Building the element conductivity matrix   
 %  \int_element (\grad N)^T Cond \grad N
 % 
 % inouts:
 % eltObj is an element object
 % Cond: is the related conductivity (scalar or matrix if tensor) 
 % outputs:
 % kel : element conductivity
 
   if isa(eltObj,'ElementTri3')
        % linear Triangle
        
        % Gauss Integration over the unit triangle
        % with a single Gauss point (Constant gradient for linear triangle)
        xeta=[1./3 1./3];
        Wl=1./2;
        [DNaDx,j]=GradN(xeta,eltObj);
          
        switch eltObj.type
            case '2D'
                kel=Wl*j*DNaDx'*Cond*DNaDx;
            case 'Axis'
                [X]=Mapx(xeta,eltObj);
                kel=(2*pi*X(1))*Wl*j*DNaDx'*Cond*DNaDx;
        end
        
    elseif isa(eltObj,'ElementTri6')
        % Linear Triangle
        
        % Gauss Integration over the unit triangle
        % with three Gauss point (linear gradient for quadratic triangle)
        xeta=[1./6 1./6;
            2./3 1./6;
            1./6 2./3];
        Wl=[1./6 1./6 1./6];
        [DNaDx,j]=GradN(xeta(1,:),eltObj);
        kel=Wl(1)*j*DNaDx'*Cond*DNaDx;
        
        switch eltObj.type
            case '2D'
                kel=Wl(1)*j*DNaDx'*Cond*DNaDx;
            case 'Axis'
                error('Axis-symmetrie not yet implemented');
        end
        
        for i = 2:length(Wl)
            [DNaDx,j]=GradN(xeta(i,:),eltObj);
          
            switch eltObj.type
                case '2D'
                    kel=kel + Wl(i)*j*DNaDx'*Cond*DNaDx;
                case 'Axis'
                    error('Axis-symmetrie not yet implemented');
            end
        end
    else
        error('Element not yet implemented');
    end
    
end
