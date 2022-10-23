function [ur,uz]=BoussinesqSolution_PointForce(r,z,E,nu)
% Solution for the displacement due to an unit point force on the half-space !! 
%
% this solution should be valid in the far field (away from the foundation)

R=sqrt(r.^2+z.^2);

ur=((1+nu)/(2*pi*E)*(1./R)).*((r.^2).*abs(z)./(R.^3) -(1-2*nu)*(1.-abs(z)./R));

uz =((1+nu)/(2*pi*E)*(1./R)).*(2*(1-nu)+(z.^2)./(R.^2)); 


end
