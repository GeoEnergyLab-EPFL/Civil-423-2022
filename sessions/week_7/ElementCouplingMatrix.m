function [Ceel]=ElementCouplingMatrix_exo(eltObj_E,eltObj_P,alpha)
% Element level function
% Building the element coupling matrix
%  \int_element ( B)^T alpha  N
%
% inouts:
% eltObj is an element object
% rho: is the related to density/storage (scalar)
% outputs:
% kel : element mass (not lumped)

if isa(eltObj_E,'ElementTri3')
    % Linear Triangle
    
    % Gauss Integration over the unit triangle
    %  Mass Matrix -> 3 Gauss points here for the linear triangle

    xeta=[1/6 1/6;
        2/3 1/6;
        1/6 2/3;
        ];
    
    Wl=[1/6 1/6 1/6];
    
    switch eltObj_E.type
        case '2D'
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>           



%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>            
        case 'Axis'
            [B,j]=Bamat(xeta(1,:),eltObj_E);
            Baux = B(1,:) + B(2,:) + B(4,:);
            [X]=Mapx(xeta(1,:),eltObj_E);
            [N]=Na(xeta(1,:),eltObj_P);
            Ceel=(2*pi*X(1))*Wl(1)*j*(Baux'*alpha*N);
            for i=2:length(Wl)
                [B,j]=Bamat(xeta(i,:),eltObj_E);
                Baux = B(1,:) + B(2,:) + B(4,:);
                [X]=Mapx(xeta(i,:),eltObj_E);
                [N]=Na(xeta(i,:),eltObj_P);
                Ceel=Ceel+(2*pi*X(1))*Wl(i)*j*(Baux'*alpha*N);
            end
    end
    
elseif isa(eltObj_E,'ElementTri6')
    
    xeta = [1/6,0.7886751346;
         0.6220084679,0.2113248654;
         0.4465819874*10^(-1),0.7886751346;
         1/6,0.2113248654];

    Wl= [0.5283121635*10^(-1);
             0.1971687836;
             0.5283121635*10^(-1);
             0.1971687836];
    
    switch eltObj_E.type
        case '2D'
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>           




%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
        case 'Axis'
            [B,j]=Bamat(xeta(1,:),eltObj_E);
            Baux = B(1,:) + B(2,:) + B(4,:);
            [X]=Mapx(xeta(1,:),eltObj_E);
            [N]=Na(xeta(1,:),eltObj_P);
            Ceel=(2*pi*X(1))*Wl(1)*j*(Baux'*alpha*N);
            for i=2:length(Wl)
                [B,j]=Bamat(xeta(i,:),eltObj_E);
                Baux = B(1,:) + B(2,:) + B(4,:);
                [X]=Mapx(xeta(i,:),eltObj_E);
                [N]=Na(xeta(i,:),eltObj_P);
                Ceel=Ceel+(2*pi*X(1))*Wl(i)*j*(Baux'*alpha*N);
            end
    end
    
else
    error('Element not yet implemented');
end

end

