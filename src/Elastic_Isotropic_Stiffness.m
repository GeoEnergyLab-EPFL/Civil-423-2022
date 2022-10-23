function L_elas=Elastic_Isotropic_Stiffness(k,g,Geo)
% k :: bulk modulus
% g :: shear modulus
% Geo :: 'PlaneStrain' or 'Axisymmetry'
% create Elastic Stifnness matrix from bulk and shear moduli and type of
% problem

La=k+(4./3.)*g;
Lb=k-(2./3.)*g;

switch Geo
    case 'PlaneStrain'
        
        L_elas=[La Lb 0. ;...
            Lb La 0.;...
            0 0  g];  % g here not 2*g
        % implement PlaneStress here
        
    case 'Axisymmetry'
        
        L_elas=[La Lb 0. Lb ;...
            Lb La 0. Lb;...
            0 0  g 0.;... % g here not 2*g
            Lb Lb 0. La];
        
        %     Lax=[L(1,1) L(1,2) L(1,4) L(1,3) ;...
        %                                 L(2,1) L(2,2) L(2,4)    L(2,3) ; ...
        %                                 L(4,1)  L(4,2)  L(4,4)/2.  L(4,3)  ;...
        %                                 L(3,1)  L(3,2)  L(3,4)      L(3,3);];
        
    otherwise
        
        L_elas=[La Lb Lb 0 0 0 ;...
            Lb La Lb 0 0 0; ...
            Lb Lb La 0 0 0 ; ...
            0 0 0 2*g 0 0 ;...
            0 0 0 0 2*g 0 ;...
            0 0 0 0 0 2*g];
        
end
end