%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%               Plain strain poroelastic closure              %%%%%%%
%%%%%%%                                                             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all
close all

%% Material properties
%%All stiffnesses are given in [MPa] and correspond to a ohio Sandstone.

k=8.4e3;            % elastic drained bulk modulus [MPa]
g=6.8e3;            % shear modulus [MPa]
b=0.707692;         % Biot coefficient 
M=9.18478e3;        % Biot Modulus
k_u=k+b^2*M;        % undrained bulk modulus [MPa]
perm = 0.137549e-3;  % permeability value adjusted such that c= k /(mu S)=1 
                    % , where S is storage coefficient
B =(k_u-k)/(b*k_u);
mu_f=1;             % fluid viscosity [Pa s]
rho=1;              % density
kappa = perm / mu_f;% conductivity 

nu_u=(1-2*(g)/(3*k_u))/(2+2*(g)/(3*k_u));   % undrained poisson's ratio
E_u = (9*(k_u))/(1+3*(k_u)/(g));            % undrained youngs modulus           
params_u = [E_u, nu_u];

nu=(1-2*(g)/(3*k))/(2+2*(g)/(3*k));         % poisson's ratio
E=(9*(k))/(1+3*(k)/(g));                    % youngs modulus
params = [E, nu];

c = 2 * kappa * g * (1-nu)*(nu_u-nu) / ...
    (b^2 * (1-2*nu)^2 * (1-nu_u));

eta = b*(1-2*nu_u)/(2*(1-nu_u)); % poroelastic eta parameter

beta = 1e-5; % fluid compressibility

%% Calculate initial fracture variables
Vgiven = [0.01]%[0.01 0.02];
nEl = [10,100,1000]%,2500,5000];
KIc = 1.5;
Kp = 4*(2/pi())^(1/2)* KIc;
Ep = params(1) / (1-params(2)^2);
Ep_u = params_u(1) / (1-params_u(2)^2); 

%% Prepare error values
RMSE_op = zeros(length(nEl),length(Vgiven));
RMSE_po = zeros(length(nEl),length(Vgiven));
RMSE_sx = zeros(length(nEl),length(Vgiven));
fig_ind = 1;
%% Loop over the number of elements
for iterV = 1:length(Vgiven)
    Vo = Vgiven(iterV);

    Radius = 2/pi()^(2/3)*(Ep * Vo / Kp)^(2/3);

    p_f = pi()^(1/3)/8 * (Kp^4 / (Vo * Ep))^(1/3);


    %% Mesh generation
    %Radius = 1.;    % Radius of fracture
    l_area = floor(10*Radius);    % length of observed area
    d_area = floor(10*Radius);    % Depth of observed area
    
    for iterN = 1:length(nEl)
        clear nfr external_approach fracture_zone meshE meshP
        nfr = nEl(iterN);

        external_approach = [flip(Radius.*logspace(0,log10(2),nfr));
                    -d_area*ones(1,nfr)]';
        fracture_zone = flip([(1-flip(logspace(0,log10(2),nfr)-1)).*Radius;
                    -d_area*ones(1,nfr)],2)';

        node = unique([ 0., 0. ;
            l_area, 0. ;
            l_area, -d_area;
            external_approach;
            Radius, -d_area;
            fracture_zone],'rows','stable');
        node(end,1) = 0.;

        edge = [];
        for e=1:length(node)
            edge = [ edge ; e e+1 ];
        end
        edge(end,2) = 1;

        %-- call mesh-gen.
        [meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = ...
            refine2(node,edge,[],[],4) ;
        [meshL.nodes,meshL.edge, meshL.connectivity,meshL.id] = smooth2(meshL.nodes,...
            meshL.edge, meshL.connectivity,meshL.id); % Refine the mesh
        %-- Transform into quadratic mesh
        [meshQ.nodes,meshQ.connectivity] = Tri3ToTri6(meshL);
        meshQ.edge = meshL.edge;
        meshQ.id = meshL.id;

        % Decide which mesh to use for elasticity (meshE) and pressure (meshP)
        meshE = meshQ;
        meshP = meshL;

        % Plot the mesh and higlight the boundary.
        figure(fig_ind);
        plotmesh(meshE.nodes,meshE.connectivity(:,1:3),[.2 .2 .2],'w')
        patch('faces',edge(:,1:2),'vertices',node, ...
            'facecolor','w', ...
            'edgecolor',[.1,.1,.1], ...
            'linewidth',1.5) ;
        fig_ind = fig_ind + 1;

        %% Boundary conditions

        %----- Boundaries
        % find the corresponding nodes
        bottom = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)>=Radius);
        bottom_P = find(meshP.nodes(:,2)==-d_area & meshP.nodes(:,1)>=Radius);
        left = find(meshE.nodes(:,1)==0. & meshE.nodes(:,2)>-d_area );
        top = find(meshE.nodes(:,2)==0);
        right = find(meshE.nodes(:,1)==l_area);

        %----- Displacements
        % adjust the fixed dofs and nodes
        fixed_dof_left = [2*(left-1)+1];
        fixed_dof_bottom = [2*bottom];
        fixed_nodes = unique([bottom;left]);
        fixed_dof = unique([fixed_dof_left;fixed_dof_bottom]);

        %----- Pressure
        fracture_E = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)<Radius);
        fracture_P = find(meshP.nodes(:,2)==-d_area & meshP.nodes(:,1)<=Radius);

        %----- Plotting the boundarys onto the initial plot of the mesh
        % Boundaries where displacements are fixed will be red, such that are fixed
        % in pressure blue and the ones fixed in tractions green
        hold on;
        plot(meshE.nodes(bottom,1),meshE.nodes(bottom,2),'or')
        hold on;
        plot(meshE.nodes(left,1),meshE.nodes(left,2),'or')
        hold on;
        plot(meshP.nodes(fracture_P,1),meshP.nodes(fracture_P,2),'*b')
        hold on;
        plot(meshE.nodes(right,1),meshE.nodes(right,2),'sg')
        hold on;
        plot(meshE.nodes(top,1),meshE.nodes(top,2),'sg')
        hold off;
        %% Initial condition
        %%As an initial condition we will set a unifrom pressure field with
        %%deviatoric stresses where s_x = -20 and s_y = -40 but without shear
        %%stresses. Additionally we will define an initial pore pressure field
        %%which is uniform as well.

        %----- Initial stress field
        % Set the initial stressfield in MPa
        %
        Po = 30; %mean compressive 
        So = 0; % deviatoric compressive
        Sig_unif = -[(Po-So);(Po+So);0];  %  [MPa]
        % Assemble the stress vector at the nodes
        Sig_o = SetStressField(meshE,'2D',...
            {[1:length(meshE.nodes(:,1))],Sig_unif});

        %----- Initial pore pressure field
        % Define the fixed and free nodes
        eq_free_p = setdiff([1:length(meshP.nodes(:,1))],fracture_P)';
        eq_fixed_p = setdiff([1:length(meshP.nodes(:,1))],eq_free_p)';

        % Define the  initial pore pressure field
        p_o = 5; %[MPa]
        p_f = Po + p_f;
        p_fixed = sparse(length(meshP.nodes(:,1)),1);
        p_fixed(eq_free_p) = p_o;
        p_fixed(eq_fixed_p) = p_f;
        %% Set boundary tractions
        %%We will need to set the corresponding surface tractions on the "free"
        %%(i.e. no displacment fixed) edges on the top and right. We will perform
        %%this by calculating the equivalent surface tractions from the initial
        %%stress field and apply them. Note that you cansimply check if they
        %%correspond to the stresses given in the section before.

        %----- Get the corresponding tractions
        f_y_top = sum(Sig_o(top*2))/l_area;
        f_x_right = sum(Sig_o(right*2-1))/d_area;
        f_y_fracture = p_f;

        %----- Assemble the force vector
        f_s=AssembleTractionsOverLine(meshE,right,[f_x_right,0],'2D') + ...
            AssembleTractionsOverLine(meshE,top,[0,f_y_top],'2D') + ...
            AssembleTractionsOverLine(meshE,fracture_E,[0,f_y_fracture],'2D') - ...
            Sig_o;

        %% Assembling the different matrices and vectors

        % Stiffness matrix
        [K]=AssembleStiffnessMatrix(meshE,params,'2D');

        % mass matrix term
        [Storage]=AssembleMassMatrix(meshP,1/M,'2D');

        % laplacian term
        [C]=AssembleConductivityMatrix(meshP,kappa,'2D');

        % Coupling term 
        [Ce]=AssembleCouplingMatrix(meshE,meshP,b,'2D');

        %% Solving directly for the undrained solution at t=0 ! 
        %%First we solve the undrained response. Where we already set the boundary
        %%condition within the borehole to zero pressure and zero stress.
        dt= 0.;         % time step

        AA=Storage+dt*C;   % lower right part of the system matrix

        %----- Getting the different number of DOF and assembling the total matrix
        % getting DOF
        ntot_P = length(C(:,1));
        ntot_E = length(K(:,1));
        ntot=ntot_E+ntot_P;

        % Assemble the matrix
        TotMat=sparse(ntot,ntot);
        TotMat(1:ntot_E,1:ntot_E)=K;
        TotMat(ntot_E+1:ntot,ntot_E+1:ntot)=-AA;
        TotMat(ntot_E+1:ntot,1:ntot_E)=-Ce';
        TotMat(1:ntot_E,ntot_E+1:ntot)=-Ce;
        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %----- Preparing the system
        % Get free equations.
        eq_free_u = setdiff([1:ntot_E],fixed_dof)';
        eq_fixed_u = setdiff([1:ntot_E],eq_free_u)';
        eq_free= [eq_free_u;eq_free_p+ntot_E];
        eq_fixed = setdiff([1:ntot],eq_free)';

        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %----- Setting the force vector to the initial conditions
        % Note that you need to account for the initial pore pressure field which
        % will influence the initial force vector.
        Ftot=sparse(ntot,1);
        Ftot(1:ntot_E) = f_s-Ce*p_fixed;
        Ftot(ntot_E+1:ntot) = -AA*p_fixed;
        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %----- Solving the system
        Undrained_sol=sparse(ntot,1);
        Undrained_sol(eq_free) = TotMat(eq_free,eq_free)\Ftot(eq_free);

        %----- Seperate displacements from pressure
        U_u=sparse(ntot_E,1);
        U_u(eq_free_u) = Undrained_sol(eq_free_u);
        pp_u = sparse(ntot_P,1);
        pp_u(eq_free_p) = Undrained_sol(ntot_E+eq_free_p);
        pp_u(eq_fixed_p) = p_fixed(eq_fixed_p);
        %% Calculating needed quantities
        % we need to reshape the mesh
        udisp = [U_u(1:2:end) U_u(2:2:end)];

        % Calculate the projected stress
        clear S
        S = ProjectStress(meshE,'2D',[E_u nu_u],U_u);

        %% Plotting the opening profile along the bottom line
        ind_tip = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1)==Radius);

        op_prof = (meshE.nodes(fracture_E,:)+...
                    udisp(fracture_E,:))';
        op_prof(2,:) = op_prof(2,:)+d_area;

        figure(fig_ind)
        scatter(op_prof(1,:)./Radius,op_prof(2,:)/((2/(E_u))*(p_f-Po)),'Displayname',...
            'Numerical solution (course code)')
        hold on;
        scatter(op_prof(1,:)./Radius,...
            (Radius^2-op_prof(1,:).^2).^(1/2),...
            'Displayname','Analytical solution','LineWidth',3)
        hold off;
        title('Fracutre opening at the undrained state')
        legend
        xlabel(' x [m]');
        ylabel(' w_u [m]');
        fig_ind = fig_ind + 1;

        %% Calculating the mean square error of the opening
        op_an = (Radius^2-op_prof(1,:).^2).^(1/2);
        RMSE_op(iterN,iterV) = sqrt(mean((op_prof(2,:)/((2/(E_u))*(p_f-Po)) - op_an).^2))

        %% Compare undrained pore pressure singularity

        %----- pressure undrained response \theta=0
        % define the analytical solution
        clear Analyticp_u
        Analyticp_u = @(x) -B*KIc./((x-Radius)*pi()*2).^(1/2)+p_o;

        % extract and re-arrange the numerical solution
        clear bottom_p p_bottom ia
        bottom_p = find(meshP.nodes(:,2)==-d_area & meshP.nodes(:,1) > Radius);
        [p_bottom,ia]=sort(meshP.nodes(bottom_p,1));

        % plotting the pore pressure
        figure(fig_ind)
        plot(p_bottom./Radius, pp_u(bottom_p(ia))/p_o,'-*b')
        hold on
        plot(p_bottom./Radius,Analyticp_u(p_bottom)/p_o,'-r')
        title('Pore pressure along x=-d_{area} - undrained response')
        xlabel(' x [m]');
        ylabel(' p_u [MPa]');
        hold off;
        fig_ind = fig_ind + 1;

        % plotting the error on the pore pressure along the axis.
        figure(fig_ind)
        semilogy(p_bottom./Radius,abs(pp_u(bottom_p(ia))-Analyticp_u(p_bottom))./...
            abs(Analyticp_u(p_bottom)),'b','LineWidth',2.5,...
            'Displayname','Error')
        title('Relative error on p_u along x=-d_{area} - undrained response')
        hold off;
        legend
        xlabel(' x [m] ');
        ylabel(' E(p_u)_{rel} ');
        fig_ind = fig_ind + 1;

        %% Calculating the mean square error of the opening
        clear po_an
        po_an = Analyticp_u(p_bottom)./p_o;
        RMSE_po(iterN,iterV) = sqrt(mean((pp_u(bottom_p(ia))./p_o - po_an).^2))

        %% Compare undrained stress singularity ahead of crack tip

        %----- pressure undrained response \theta=0
        % define the analytical solution
        clear AnalyticSxx_u
        AnalyticSxx_u = @(x) KIc./((x-Radius)*pi()*2).^(1/2);

        % extract and re-arrange the numerical solution
        clear bottom_s S_bottom is
        bottom_S = find(meshE.nodes(:,2)==-d_area & meshE.nodes(:,1) > Radius);
        [S_bottom,is]=sort(meshE.nodes(bottom_S,1));

        % plotting the pore pressure
        figure(fig_ind)
        plot(S_bottom./Radius, S(bottom_S(is),2)./Po,'-*b')
        hold on
        plot(S_bottom./Radius,AnalyticSxx_u(S_bottom)./Po,'-r')
        hold off;
        title('Horizontal stress along x=-d_{area} - undrained response')
        xlabel(' x [m]');
        ylabel(' S_{xx} [MPa]');
        fig_ind = fig_ind + 1;

        % plotting the error on the pore pressure along the axis.
        figure(fig_ind)
        semilogy(S_bottom./Radius,abs(S(bottom_S(is),1)-AnalyticSxx_u(S_bottom))./...
            abs(Analyticp_u(S_bottom)),'b','LineWidth',2.5,...
            'Displayname','Error')
        title('Relative error on S_{xx} along x=-d_{area} - undrained response')
        legend
        xlabel(' x [m] ');
        ylabel(' E(S_{xx})_{rel} ');
        fig_ind = fig_ind + 1;

        %% Calculating the mean square error of the opening
        clear sy_an
        sy_an = Analyticp_u(S_bottom)./Po;
        RMSE_sx(iterN,iterV) = sqrt(mean((S(bottom_S(is),2)./Po - sy_an).^2))
    end
end