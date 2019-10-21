classdef ElementSeg3
    
    properties
        
        xae;  % --
        type;
        
        nnodes=3;
        
    end
    
    
    methods
        
        function  obj = ElementSeg3(xae,type)
            
            if length(xae)==3 % must actually be 3*1
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
           %%% would have to be coded up !
            
        end
        
        
        function [j]=Jacobian(xil,obj)
            
            DNaDxi=[xil-0.5,-2*xil,xil+0.5];
            
            DxDxi=DNaDxi*(obj.xae);
            j=DxDxi;
            
        end
        
        
        %----------------------------------------------------------------
        function [Xaux]=Mapx(xil,obj)
            
            Na_xi=[.5*xil*(xil-1) (1+xil)*(1-xil) .5*xil*(xil+1)];
            Xaux=Na_xi*(obj.xae);
            
        end
        
        %----------------------------------------------------------------
        function [Na]=Na(xil,obj)
            
            Na=[.5*xil*(xil-1) (1+xil)*(1-xil) .5*xil*(xil+1)] ;
            
        end
        
    end
    
    
    
end
