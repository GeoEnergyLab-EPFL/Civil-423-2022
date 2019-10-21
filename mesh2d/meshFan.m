function [base] = meshFan(varargin)
%meshFan Function generating a fan within the mesh at a spc-
%ified location

%---------------------------------------------- extract args
    opts = [] ;
    base = [];

    if (nargin>=+1), base = varargin{1}; end
    if (nargin>=+2), opts = varargin{2}; end
    
%---------------------------------------------- check ops
    if ~isfield(opts,'center')
        error('meshFan: No center for fan provided') ;
    end
    
    if ~isfield(opts,'n_fan')
        opts.nfan = 10;
    end
    
    if ~isfield(opts,'l_max')
        opts.l_max = 0.1;
    end
    
    if ~isfield(opts,'angle')
        opts.angle = [0,2*pi];
    end
    
    if ~isfield(opts,'r_fan')
        opts.angle = 1.;
    end
    
    if ~isfield(opts,'boundary')
        opts.boundary = false;
    end
    
    if ~opts.boundary
        [~,on] = inpolygon(opts.center(1),opts.center(2),...
            base.node(:,1)',base.node(:,2)');
        if on
            opts.boundary = true;
        end
    end
    
%-------------------------------------------- Calculations
   
%-- Change parameters if the center is on the boundary
if opts.boundary
    n_corr = 1;
    n_liopts.n_fan = opts.n_fan + 2;
    eval = 2:n_liopts.n_fan-1;
elseif opts.angle == [0,2*pi]
    n_corr = 1;
    n_liopts.n_fan = opts.n_fan+n_corr;
    eval = 2:n_liopts.n_fan;
else
    n_corr = 0;
    n_liopts.n_fan = opts.n_fan+n_corr;
    eval = 1:n_liopts.n_fan;
end

%-- Calculate the new nodes
endpoints = zeros(opts.n_fan,2);
fan_angle = linspace(opts.angle(1),opts.angle(2),n_liopts.n_fan);
N_p = ceil(opts.r_fan/opts.l_max);
radi = [0:1/N_p:1].^2;
Nodes = [];
for i = eval
    endpoints(i-n_corr,:) = opts.center + opts.r_fan.*[cos(fan_angle(i)),...
        sin(fan_angle(i))];
    x= (1-radi).*opts.center(1)+radi.*endpoints(i-n_corr,1);
    y= (1-radi).*opts.center(2)+radi.*endpoints(i-n_corr,2);
    %-- Check if the nodes are inside the boundary
    [in] = inpolygon(x,y,...
            base.node(:,1)',base.node(:,2)');
    out = find(in == 0);
    newNodes = [x',y'];
    %-- Remove nodes outside boundary and add new points
    newNodes(out,:) = [];
    if isempty(newNodes)
        endpoints(i-n_corr,1) = NaN;
    else
        Nodes = [Nodes;
            newNodes]; 
        endpoints(i-n_corr,:) = newNodes(end,:);
    end
end

%------------------------------------------- Assemble edges

endpoints(isnan(endpoints(:,1)),:) = [];
%-- find indexes of endpoints
endpoints_ind = zeros(length(endpoints(:,1)),1);
for i = 1:length(endpoints(:,1))
    endpoints_ind(i) = find(Nodes(:,1) == endpoints(i,1)...
        & Nodes(:,2) == endpoints(i,2));
    %-- generate edge connectivity
    if i == 1
        edges = [1:endpoints_ind(i)-1;2:endpoints_ind(i)]';
    else
        edges = [edges;
            [endpoints_ind(i-1)+1:endpoints_ind(i)-1;...
            endpoints_ind(i-1)+2:endpoints_ind(i)]'];
    end
    
end

%------------------------------------------- Assemble system
base.edge = [base.edge;edges+length(base.node(:,1))];
base.node = [base.node;Nodes];


end

