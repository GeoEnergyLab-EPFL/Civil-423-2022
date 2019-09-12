function [kel]=ElementMassMatrix(eltObj,rho)
% Element level function
% Building the element conductivity matrix
%  \int_element (  N)^T rho  N
%
% inouts:
% eltObj is an element object
% rho: is the related to density/storage (scalar)
% outputs:
% kel : element mass (not lumped)

if isa(eltObj,'ElementTri3')
    % Linear Triangle
    
    % Gauss Integration over the unit triangle
    %  Mass Matrix -> 3 Gauss points here for the linear triangle
    
    xeta=[1./6 1./6;
        2./3 1./6;
        1./6 2./3.;
        ];
    
    Wl=[1./6 1./6 1./6];
    
    [j]=Jacobian(xeta(1,:),eltObj);
    
    switch eltObj.type
        case '2D'
            
            [N]=Na(xeta(1,:),eltObj);
            kel=Wl(1)*j*(N'*rho*N);
            for i=2:length(Wl)
                [N]=Na(xeta(i,:),eltObj);
                kel=kel+Wl(i)*j*(N'*rho*N);
            end
            
        case 'Axis'
            
            [X]=Mapx(xeta(1,:),eltObj);
            [N]=Na(xeta(1,:),eltObj);
            
            kel=(2*pi*X(1))*Wl(1)*j*(N'*rho*N);
            for i=2:3
                [N]=Na(xeta(i,:),eltObj);
                [X]=Mapx(xeta(i,:),eltObj);
                kel=kel+(2*pi*X(1))*Wl(i)*j*(N'*rho*N);
            end
    end
    
else
    error('Element not yet implemented');
end

end



%  for TRI6
%         xeta = [ 1./3 1./3; 1./5 1./5;3./5 1/5; 1/5 3/5];
%         Wl= [-27/96; 25/96;25/96;25/96;];

