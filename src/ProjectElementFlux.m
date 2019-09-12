function [fel_x,fel_y]=ProjectElementFlux(eltObj,Cond,usol)
% Element level function
% Building the element conductivity matrix
% - \int_element N^T . (Cond . \grad N . U )
%
% inouts:
% eltObj is an element object
% Cond: is the related conductivity (scalar or matrix if tensor)
% usol : is the solution (of the scalar laplacian) at the 3 nodes of the triangle
% outputs:
% fel : element flux force vector for projection system

if isa(eltObj,'ElementTri3')
    % Linear Triangle
    
    
    % Gauss Integration over the unit triangle
    %  Linear function -> 3 Gauss points here for the linear triangle
    
    xeta=[1./6 1./6;
        2./3 1./6;
        1./6 2./3.;
        ];
    
    Wl=[1./6 1./6 1./6];
    
    [j]=Jacobian(xeta(1,:),eltObj);
    
    
    switch eltObj.type
        case '2D'
            
            [N]=Na(xeta(1,:),eltObj);
            [DNaDx,j]=GradN(xeta,eltObj);
            qaux=-Cond*DNaDx*usol;
            fel_x=Wl(1)*j*(N'*qaux(1));
            fel_y=Wl(1)*j*(N'*qaux(2));
            
            for i=2:length(Wl)
                [N]=Na(xeta(i,:),eltObj);
                [DNaDx,j]=GradN(xeta,eltObj);
                qaux=-Cond*DNaDx*usol;
                %disp(size(qaux));
                
                fel_x=fel_x+Wl(i)*j*(N'*qaux(1));
                fel_y=fel_y+Wl(i)*j*(N'*qaux(2));
                
            end
            
        case 'Axis'
            
            [X]=Mapx(xeta(1,:),eltObj);
            [N]=Na(xeta(1,:),eltObj);
            [DNaDx,j]=GradN(xeta,eltObj);
            qaux=-Cond*DNaDx*usol;
            fel_x=(2*pi*X(1))*Wl(1)*j*(N'*qaux(1));
            fel_y=(2*pi*X(1))*Wl(1)*j*(N'*qaux(2));
            
            for i=2:3
                [N]=Na(xeta(i,:),eltObj);
                [DNaDx,j]=GradN(xeta,eltObj);
                [X]=Mapx(xeta(i,:),eltObj);
                qaux=-Cond*DNaDx*usol;
                fel_x=fel_x+(2*pi*X(1))*Wl(i)*j*(N'*qaux(1));
                fel_y=fel_y+(2*pi*X(1))*Wl(i)*j*(N'*qaux(2));
            end
    end
    
else
    error('Element not yet implemented');
end

end
