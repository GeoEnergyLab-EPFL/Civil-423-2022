function [S_el]=SetElementStressField(eltObj,s)
 
   if isa(eltObj,'ElementTri3')
        xeta=[1./6 1./6;
            2./3 1./6;
            1./6 2./3.;
            ];

        Wl=[1./6 1./6 1./6];
        
        [B,j]=Bamat(xeta(1,:),eltObj);
        
        switch eltObj.type
            case '2D'
                S_el=Wl(1)*j*B'*s;
                for i=2:length(Wl)
                    [B,j]=Bamat(xeta(i,:),eltObj);
                    S_el=S_el + Wl(i)*j*B'*s;
                end
            case 'Axis'
                [X]=Mapx(xeta(1,:),eltObj);
                S_el=(2*pi*X(1))*Wl(1)*j*B'*s;
                for i=2:length(Wl)
                    [B,j]=Bamat(xeta(i,:),eltObj);
                    [X]=Mapx(xeta(i,:),eltObj);
                    S_el=S_el + (2*pi*X(1))*Wl(i)*j*B'*s;
                end
        end

   elseif isa(eltObj,'ElementTri6')
        % Quadratic Triangle
        
        xeta = [1/6,0.7886751346;
             0.6220084679,0.2113248654;
             0.4465819874*10^(-1),0.7886751346;
             1/6,0.2113248654];

        Wl= [0.5283121635*10^(-1);
                 0.1971687836;
                 0.5283121635*10^(-1);
                 0.1971687836];
        
        [B,j]=Bamat(xeta(1,:),eltObj);
        
        switch eltObj.type
            case '2D'
                S_el=Wl(1)*j*B'*s;
                for i=2:length(Wl)
                    [B,j]=Bamat(xeta(i,:),eltObj);
                    S_el=S_el + Wl(i)*j*B'*s;
                end
            case 'Axis'
                [X]=Mapx(xeta(1,:),eltObj);
                S_el=(2*pi*X(1))*Wl(1)*j*B'*s;
                for i=2:length(Wl)
                    [B,j]=Bamat(xeta(i,:),eltObj);
                    [X]=Mapx(xeta(i,:),eltObj);
                    S_el=S_el + (2*pi*X(1))*Wl(i)*j*B'*s;
                end
        end
   else
        error('Element not yet implemented');
    end
    
end