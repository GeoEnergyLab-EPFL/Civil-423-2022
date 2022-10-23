function [Sel]=GetElementStress(eltObj,D,usol)
% Function to calculate the mean of the principal stresses of the gauss
% point
%   Sel = D B usol
% inputs:
% mesh :: eltObj containing the local element
% D :: the constitutive matrix
% usol : the nodal values of the displacement field of the element
%
% outputs:
%   Sel :: The local principal stresses (x,y,tau_xy) for '2D' and
%   (r,z,tau_rz,theta) for 'Axis'

   if isa(eltObj,'ElementTri3')
        % Constant strain Triangle
        
        %  Single Gauss point (Constant gradient for linear triangle)
        xeta=[1./3 1./3];
        % get B matrix
        [B,~]=Bamat(xeta,eltObj);
        % compute the solution
        Sel= D*B*usol;
        % Note: As there is only one gauss point the value is automatically
        % corresponding to the mean value.
   elseif isa(eltObj,'ElementTri6')
        % Linear strain Triangle
        
        %  three Gauss points (linear gradient for quadratic triangle)
        xeta=[1./6 1./6;
            2./3 1./6;
            1./6 2./3];
        
        % get B matrix
        [B,~]=Bamat(xeta(1,:),eltObj);
        % calculate solution
        Sel = (D*B*usol)';
        
        %-- loop over the gauss points
        for i = 2:length(xeta(:,1))
            % get B matrix
            [B,~]=Bamat(xeta(i,:),eltObj);
            % calculate solution
            Sel=[Sel;
                (D*B*usol)'];
        end
        %-- calculate mean value
        Sel = mean(Sel);
   else
        error('Element not yet implemented');
    end
    
end
