function [kel]=ElementStiffnessMatrix(eltObj,D)
% Element level function
% Building the element stifnness matrix
%  \int_element (  B)^T D  B
%
% inouts:
% eltObj is an element object
% D: elastic stiffness coefficient matrix
% outputs:
% kel : element stiffness  
    
    if isa(eltObj,'ElementTri3')
        % Constant Strain Triangle
        
        % Gauss Integration over the unit triangle
        % with a single Gauss point (Constant strain)
        xeta=[1./3 1./3];
        Wl=1./2;
        
        % get the B matrix
        [B,j]=Bamat(xeta,eltObj);
        
        switch eltObj.type
            case '2D'
                % compute element stiffness matrix
                kel=Wl*j*B'*D*B;
            case 'Axis'
                % compute element stiffness matrix
                [X]=Mapx(xeta,eltObj);
                kel=(2*pi*X(1))*Wl*j*B'*D*B;
        end
        
    elseif isa(eltObj,'ElementTri6')
        % Linear Strain Triangle
        
        % Gauss Integration over the unit triangle
        % with three Gauss point (Linear strain)
         xeta=[1/6 1/6;
             2/3 1/6;
             1/6 2/3];
         Wl=[1/6 1/6 1/6];
        
        switch eltObj.type
            case '2D'
                % get the b matrix
                [B,j]=Bamat(xeta(1,:),eltObj);
                % compute element stiffness matrix
                kel=Wl(1)*j*B'*D*B;
                % loop over the gauss points
                for i = 2:length(xeta(:,1))
                    % get the b matrix
                    [B,j]=Bamat(xeta(i,:),eltObj);
                    % compute element stiffness matrix
                    kel=kel + Wl(i)*j*B'*D*B;
                end
            case 'Axis'
                % get the b matrix
                [B,j]=Bamat(xeta(1,:),eltObj);
                % compute element stiffness matrix
                [X]=Mapx(xeta(1,:),eltObj);
                kel=(2*pi*X(1))*Wl(1)*j*B'*D*B;
                % loop over the gauss points
                for i = 2:length(xeta(:,1))
                    % get the b matrix
                    [B,j]=Bamat(xeta(i,:),eltObj);
                    [X]=Mapx(xeta(i,:),eltObj);
                    % compute element stiffness matrix
                    kel=kel + (2*pi*X(1))*Wl(i)*j*B'*D*B;
                end
        end
        
    else
        error('Element not yet implemented');
    end
    
end
