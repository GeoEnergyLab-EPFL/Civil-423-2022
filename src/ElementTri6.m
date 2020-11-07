%  Class containing all necessary functions for FE
%  Triangular element - linear interpolation - 3 nodes
classdef ElementTri6
    
    properties
        
        xae;   % matrix containing the coordinates of the nodes - size
        type; % '2D' or 'Axis'
        
        nnodes=6;
        
    end
    
    methods
        % constructor
        function obj=ElementTri6(xae,type) %ncon_nodes,
            % todo : enforce a check
            if (size(xae,1)~=6)
                error('Error in element coordinates !');
            end
            obj.xae=xae;
            obj.type = type;
            
        end
        
        % B matrix - for Strain
        function [B,j]=Bamat(xil,obj)
            
            DNaDxi = [4*xil(1)-1, 0,    -3+4*xil(1)+4*xil(2), 4*xil(2),...
                -4*xil(2), 4*(1-2*xil(1)-xil(2));...
                0,4*xil(2)-1,-3+4*xil(2)+4*xil(1),...
                4*xil(1),4*(1-2*xil(2)-xil(1)),-4*xil(1);];
            
            J=DNaDxi*(obj.xae); % Jacobian matrix
            DxiDx = inv(J);  % inverse
            j=det(J);
            
            DNaDx = DxiDx*(DNaDxi);
            
            switch obj.type
                
                case '2D'
                    % exx, eyy, 2 exy
                    
                    % 3 by 12 mat
                    B = [DNaDx(1,1)  0 DNaDx(1,2) 0 DNaDx(1,3) 0 ...
                        DNaDx(1,4)  0 DNaDx(1,5)  0 DNaDx(1,6)  0;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) ...
                        0   DNaDx(2,4) 0 DNaDx(2,5) 0 DNaDx(2,6);...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) ...
                        DNaDx(2,3) DNaDx(1,3) DNaDx(2,4) DNaDx(1,4) ...
                        DNaDx(2,5) DNaDx(1,5) DNaDx(2,6) DNaDx(1,6)
                        ];
                    
                case 'Axis'
                    xaux=Mapx(xil,obj);
                    [Nax]=Na(xil,obj);
                    
                    % 4 by 12 mat
                    B = [DNaDx(1,1)  0 DNaDx(1,2) 0 DNaDx(1,3) 0 ...
                        DNaDx(1,4)  0 DNaDx(1,5) 0 DNaDx(1,6) 0;...
                        0   DNaDx(2,1) 0 DNaDx(2,2) 0 DNaDx(2,3) ...
                        0   DNaDx(2,4) 0 DNaDx(2,5) 0 DNaDx(2,6) ;...
                        DNaDx(2,1) DNaDx(1,1) DNaDx(2,2) DNaDx(1,2) DNaDx(2,3) DNaDx(1,3) ...
                        DNaDx(2,4) DNaDx(1,4) DNaDx(2,5) DNaDx(1,5) DNaDx(2,6) DNaDx(1,6) ;...
                        Nax(1)/xaux(1)  0.  Nax(2)/xaux(1) 0 Nax(3)/xaux(1) 0. ...
                        Nax(4)/xaux(1)  0.  Nax(5)/xaux(1) 0 Nax(6)/xaux(1) 0. ];
            end
            
        end
        
        % x from xil
        function [Xaux]=Mapx(xil,obj)
            
            Na_xi = [xil(1)*(xil(1)*2-1), xil(2)*(xil(2)*2-1),...
                (1-xil(1)-xil(2))*(2*(1-xil(1)-xil(2))-1),  4*xil(1)*xil(2),...
                4*xil(2)*(1-xil(1)-xil(2)),4*xil(1)*(1-xil(1)-xil(2))];
            % % %
            Xaux=[Na_xi*(obj.xae(:,1)) Na_xi*(obj.xae(:,2))];
            
        end
        
        % shape function - for mass matrix etc.
        function [Na]=Na(xil,obj)
            
            Na = [xil(1)*(xil(1)*2-1),xil(2)*(xil(2)*2-1),...
                (1-xil(1)-xil(2))*(2*(1-xil(1)-xil(2))-1), ...
                4*xil(1)*xil(2),4*xil(2)*(1-xil(1)-xil(2)),...
                4*xil(1)*(1-xil(1)-xil(2))];
            
        end
        
        % for Laplacian
        function [DNaDx,j]=GradN(xil,obj)
            
             DNaDxi = [4*xil(1)-1, 0,    -3+4*xil(1)+4*xil(2), 4*xil(2),...
                -4*xil(2), 4*(1-2*xil(1)-xil(2));...
                0,4*xil(2)-1,-3+4*xil(2)+4*xil(1),...
                4*xil(1),4*(1-2*xil(2)-xil(1)),-4*xil(1);];
            
            J=DNaDxi*(obj.xae); % Jacobian matrix
            DxiDx = inv(J);  % inverse
            j=det(J);
            
            DNaDx = DxiDx*(DNaDxi);
            
        end
        
        % Det jacobian
        function [j]=Jacobian(xil,obj)
            
            DNaDxi = [4*xil(1)-1, 0,    -3+4*xil(1)+4*xil(2), 4*xil(2),...
                -4*xil(2), 4*(1-2*xil(1)-xil(2));...
                0,4*xil(2)-1,-3+4*xil(2)+4*xil(1),...
                4*xil(1),4*(1-2*xil(2)-xil(1)),-4*xil(1);];
            
            J=DNaDxi*(obj.xae); % Jacobian matrix
            DxiDx = inv(J);  % inverse
            j=det(J);
            
        end
        
    end
    
    
end
