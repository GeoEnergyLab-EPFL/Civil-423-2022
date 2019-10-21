function [fel]=ElementNeumann(eltObj,t)
% this function perform the 1D integration
% along an edge for a Neumann BC
%    \int_\Gamma_e     N t  dx
%
% currently only 2D problems are tackled by the library, so this
% is a 1D integration
%
%  t :: scalar component assumed constant over the element here ....for now
%  
% note also for axis problem - only segment along e_r are properly
% integrated.... with this procedure => needs to be generalized ....

if (~isa(eltObj,'ElementSeg3') && (~isa(eltObj,'ElementSeg2') ) )
    error('Note this routine performs a 1D integration !');
    
end

if isa(eltObj,'ElementSeg2')
    
    switch eltObj.type
        
        case '2D'
            
            xeta = [0.];
            Wl =[2];
            
            [j]=Jacobian(xeta(1),eltObj);
            
            [N]=Na(xeta(1),eltObj);
            fel=Wl(1)*j*(N'*t);
            
            
        case 'Axis'
            
            xeta = [-1./sqrt(3) 1./sqrt(3)];
            Wl= [1 1];
            
            [j]=Jacobian(xeta(1),eltObj);
            [N]=Na(xeta(1),eltObj);
            [Xaux]=Mapx(xeta(1),eltObj);
            
            fel=2*pi*Xaux*Wl(1)*j*(N'*t);
            for i=2:length(Wl)
                [N]=Na(xeta(i),eltObj);
                [Xaux]=Mapx(xeta(i),eltObj);
                fel=fel+2*pi*Xaux*Wl(i)*j*(N'*t);
            end
            
    end
    
    
elseif isa(eltObj,'ElementSeg3')
    
    switch eltObj.type
        
        case '2D'
            
            xeta = [-sqrt(3/5.) 0.  sqrt(3/5.)];
            Wl= [5/9 8/9  5/9];
            
            [j]=Jacobian(xeta(1),eltObj);
            [N]=Na(xeta(1),eltObj);
            fel=Wl(1)*j*(N'*t);
            for i=2:length(Wl)
                [N]=Na(xeta(i),eltObj);
                fel=fel+Wl(i)*j*(N'*t);
            end
            
        case 'Axis'
            xeta=[-0.861136 -0.339981      0.339981 0.861136];
             
            Wl= [0.347855 0.652145  0.652145   0.347855 ];
            
            [j]=Jacobian(xeta(1),eltObj);
            [N]=Na(xeta(1),eltObj);
            [Xaux]=Mapx(xeta(1),eltObj);
            fel=2*pi*Xaux*Wl(1)*j*(N'*t);
            
            for i=2:length(Wl)
                [N]=Na(xeta(i),eltObj);
                [Xaux]=Mapx(xeta(i),eltObj);
                fel=fel+2*pi*Xaux*Wl(i)*j*(N'*t);
            end
            
    end
    
end


end