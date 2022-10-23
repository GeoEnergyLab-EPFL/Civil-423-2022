%  Class containing all necessary functions for FE 
%  Triangular element - linear interpolation - 3 nodes
classdef ElementTri3

    properties
        
        xae;   % matrix containing the coordinates of the nodes - size 
        type; % '2D' or 'Axis'
        
        nnodes=3;
        
    end
    
    methods
        % constructor
        function obj=ElementTri3(xae,type) %ncon_nodes,
            % todo : enforce a check 
            if (size(xae,1)~=3) 
                error('Error in element coordinates !');
            end
            obj.xae=xae;
            obj.type = type; 
            
        end
        
        % B matrix - for Strain
        function [B,j]=Bamat(xil,obj)
            
            DNaDxi=[-1 1 0;...%w.r. to xi
                  -1 0 1;];%w.r. to eta
            
            J=DNaDxi*(obj.xae); % Jacobian matrix            
            DxiDx = inv(J);  % inverse
            j=det(J);
            
            DNaDx = DxiDx*(DNaDxi); 
            
            switch obj.type
                
                case '2D'
                    % exx, eyy, 2 exy
                    
                    % 3 by 6 mat
                    B = [DNaDx(1,1)  0 DNaDx(1,2) 0 DNaDx(1,3) 0  ;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3)
                        ];
                    
                case 'Axis'
                 % err, ezz, 2 erz, ett

                    xaux=Mapx(xil,obj);
                    [Nax]=Na(xil,obj);
                    
                    % 4 by 6 mat
                    B = [DNaDx(1,1)  0 DNaDx(1,2) 0 DNaDx(1,3) 0  ;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3);...
                        Nax(1)/xaux(1)  0.  Nax(2)/xaux(1) 0. Nax(3)/xaux(1) 0. ];                    
 
            end
            
        end

         % x from xil
         function [Xaux]=Mapx(xil,obj)
            
             Na_xi=[1-xil(1)-xil(2) xil(1) xil(2)] ;
             Xaux=[Na_xi*(obj.xae(:,1)) Na_xi*(obj.xae(:,2))];  

         end
         
         % shape function - for mass matrix etc.
        function [Na]=Na(xil,obj)
            
            Na=[1-xil(1)-xil(2) xil(1) xil(2)] ;
            
        end
        
        % for Laplacian 
         function [DNaDx,j]=GradN(xil,obj)
            
            DNaDxi=[-1 1 0;...%w.r. to xi
                  -1 0 1;];%w.r. to eta
            
            J=DNaDxi*(obj.xae); % Jacobian matrix            
            DxiDx = inv(J);  % inverse
            j=det(J);
            
            DNaDx = DxiDx*(DNaDxi); 
            
         end
        
        % Det jacobian
        function [j]=Jacobian(xil,obj)
            
            DNaDxi=[-1 1 0;...%w.r. to xi
                  -1 0 1;];%w.r. to eta
            
            J=DNaDxi*(obj.xae); % Jacobian matrix            
            DxiDx = inv(J);  % inverse
            j=det(J);
            
        end
                
        % Element volume
        function [vol]=Volume(obj)
            [vol]=Jacobian([0. 0. ],obj)/2.;
        end
           
        % Element side length
        function [l]=sidelength(ind,obj)
            [l]=norm(obj.xae(ind(1),:)-obj.xae(ind(2),:));
        end
        
    end
    
    
end
