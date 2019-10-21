function [fel_x,fel_y,tauel,fel_theta]=ProjectElementStress(eltObj,D,usol)
% Element level function
% Building the element stress matrix
% \int (N^T (D b \times solution)_x )
%
% inouts:
% eltObj is an element object
% D: is the constitutive matris
% usol : is the solution of the displacement field at the nodes
% outputs:
% fel : element stress force vector for projection system

if isa(eltObj,'ElementTri3')
    % Constant strain Triangle

    % Gauss Integration over the unit triangle
    %  Linear function -> 3 Gauss points here for the CST
    xeta=[1./6 1./6;
        2./3 1./6;
        1./6 2./3.;
        ];
    
    Wl=[1./6 1./6 1./6];
    
    
    switch eltObj.type
        case '2D'
            % get the shape functions and the B matrix
            [N]=Na(xeta(1,:),eltObj);
            [B,j]=Bamat(xeta(1,:),eltObj);
            % solving the system at the gauss point
            saux= D*B*usol;
            % Project the gauss point solution onto the nodes
            fel_x=Wl(1)*j*(N'*saux(1));
            fel_y=Wl(1)*j*(N'*saux(2));
            tauel=Wl(1)*j*(N'*saux(3));
            
            %-- Loop over the gauss points
            for i=2:length(Wl)
                % get the shape functions and the B matrix
                [N]=Na(xeta(i,:),eltObj);
                [B,j]=Bamat(xeta(i,:),eltObj);
                % solving the system at the gauss point
                saux= D*B*usol;
                % Project the gauss point solution onto the nodes
                fel_x=fel_x+Wl(i)*j*(N'*saux(1));
                fel_y=fel_y+Wl(i)*j*(N'*saux(2));
                
            end
            fel_theta = 0.*tauel;
            
        case 'Axis'
            % Get the mean radius of the element
            [X]=Mapx(xeta(1,:),eltObj);
            % get the shape functions and the B matrix
            [N]=Na(xeta(1,:),eltObj);
            [B,j]=Bamat(xeta(1,:),eltObj);
            % solving the system at the gauss point
            saux= D*B*usol;
            % Project the gauss point solution onto the nodes
            fel_x=(2*pi*X(1))*Wl(1)*j*(N'*saux(1));
            fel_y=(2*pi*X(1))*Wl(1)*j*(N'*saux(2));
            tauel=(2*pi*X(1))*Wl(1)*j*(N'*saux(3));
            fel_theta=(2*pi*X(1))*Wl(1)*j*(N'*saux(4));
            % -- Loop over the gauss points
            for i=2:length(Wl)
                % get the shape functions and the B matrix
                [N]=Na(xeta(i,:),eltObj);
                [B,j]=Bamat(xeta(i,:),eltObj);
                [X]=Mapx(xeta(i,:),eltObj);
                % solving the system at the gauss point
                saux= D*B*usol;
                % Project the gauss point solution onto the nodes
                fel_x=fel_x+(2*pi*X(1))*Wl(i)*j*(N'*saux(1));
                fel_y=fel_y+(2*pi*X(1))*Wl(i)*j*(N'*saux(2));
                tauel=tauel+(2*pi*X(1))*Wl(i)*j*(N'*saux(3));
                fel_theta=fel_theta+(2*pi*X(1))*Wl(i)*j*(N'*saux(4));
            end
    end
elseif isa(eltObj,'ElementTri6')
    % Linear strain Triangle

    % Gauss Integration over the unit triangle
    %  quadratic function -> 4 Gauss points here for the LST
    xeta = [0.1666666667,0.7886751346;
         0.6220084679,0.2113248654;
         0.4465819874*10^(-1),0.7886751346;
         0.1666666667,0.2113248654];

    Wl= [0.5283121635*10^(-1);
             0.1971687836;
             0.5283121635*10^(-1);
             0.1971687836];

    switch eltObj.type
        case '2D'
            % get the shape functions and the B matrix
            [N]=Na(xeta(1,:),eltObj);
            [B,j]=Bamat(xeta(1,:),eltObj);
            % solving the system at the gauss point
            saux= D*B*usol;
            % Project the gauss point solution onto the nodes
            fel_x=Wl(1)*j*(N'*saux(1));
            fel_y=Wl(1)*j*(N'*saux(2));
            tauel=Wl(1)*j*(N'*saux(3));
            
            %-- Loop over the gauss points
            for i=2:length(Wl)
                % get the shape functions and the B matrix
                [N]=Na(xeta(i,:),eltObj);
                [B,j]=Bamat(xeta(i,:),eltObj);
                % solving the system at the gauss point
                saux= D*B*usol;
                % Project the gauss point solution onto the nodes
                fel_x=fel_x+Wl(i)*j*(N'*saux(1));
                fel_y=fel_y+Wl(i)*j*(N'*saux(2));
                
            end
            fel_theta = 0.*tauel;
            
        case 'Axis'
            % Get the mean radius of the element
            [X]=Mapx(xeta(1,:),eltObj);
            % get the shape functions and the B matrix
            [N]=Na(xeta(1,:),eltObj);
            [B,j]=Bamat(xeta(1,:),eltObj);
            % solving the system at the gauss point
            saux= D*B*usol;
            % Project the gauss point solution onto the nodes
            fel_x=(2*pi*X(1))*Wl(1)*j*(N'*saux(1));
            fel_y=(2*pi*X(1))*Wl(1)*j*(N'*saux(2));
            tauel=(2*pi*X(1))*Wl(1)*j*(N'*saux(3));
            fel_theta=(2*pi*X(1))*Wl(1)*j*(N'*saux(4));
            % -- Loop over the gauss points
            for i=2:length(Wl)
                % get the shape functions and the B matrix
                [N]=Na(xeta(i,:),eltObj);
                [B,j]=Bamat(xeta(i,:),eltObj);
                 [X]=Mapx(xeta(i,:),eltObj);
                % solving the system at the gauss point
                saux= D*B*usol;
                % Project the gauss point solution onto the nodes
                fel_x=fel_x+(2*pi*X(1))*Wl(i)*j*(N'*saux(1));
                fel_y=fel_y+(2*pi*X(1))*Wl(i)*j*(N'*saux(2));
                tauel=tauel+(2*pi*X(1))*Wl(i)*j*(N'*saux(3));
                fel_theta=fel_theta+(2*pi*X(1))*Wl(i)*j*(N'*saux(4));
            end
    end
    
else
    error('Element not yet implemented');
end

end
