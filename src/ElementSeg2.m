classdef ElementSeg2
    
    properties
        
        xae;  % --
        type;
        
        nnodes=2;
        
    end
    
    
    methods
        
        function  obj = ElementSeg2(xae,type)
            
            if length(xae)==2
                obj.xae=xae;
            else
                error('error in coordinates of 1D edge' )
            end
            
            if (strcmp(type,'Axis') || strcmp(type,'2D') )
                obj.type=type;
                
            else
                error(' Not proper type ')
            end
            
        end
        
        function [B,j]=Bamat(xil,obj)
            switch obj.type
                
                case 'Axis'  %% One dimensional One dofs axi-symmetry elasticity e_r e_t (e_z=0)
                    
                    DNaDxi=[-0.5  0.5];
                    DxDxi=DNaDxi*(obj.xae);
                    j=DxDxi;
                    
                    DNaDx=(1./j)*DNaDxi;
                    
                    raux=[ 0.5*(1-xil) 0.5*(1+xil)]*(obj.xae);
                    
                    B=[ DNaDx(1) 0 DNaDx(2) 0 ;...
                        0.5*(1-xil)/raux 0 0.5*(1+xil)/raux 0 ; ...
                        0 0 0 0 ];
                    
                otherwise    % 1D stuff... to be checked
                    
                    DNaDxi=[-0.5  0.5];
                    DxDxi=DNaDxi*(obj.xae);
                    j=DxDxi;
                    
                    DNaDx=(1./j)*DNaDxi;
                    
                    B=[ DNaDx(1)  DNaDx(2) ;];
                    
            end
            
        end
        
        
        function [j]=Jacobian(xil,obj)
            
            DNaDxi=[-0.5  0.5];
            DxDxi=DNaDxi*(obj.xae);
            j=DxDxi;
            
        end
        
        
        %----------------------------------------------------------------
        function [Xaux]=Mapx(xil,obj)
            
            Na_xi=[0.5*(1-xil)  0.5*(1+xil)];
            Xaux=Na_xi*(obj.xae);
            
        end
        
        %----------------------------------------------------------------
        function [Na]=Na(xil,obj)
            
            Na=[0.5*(1-xil)  0.5*(1+xil)];
            
        end
        
    end
    
    
    
end
