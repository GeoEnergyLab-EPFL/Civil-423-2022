function [kr,pm]=ComputeRelPerm(mesh,p_res,eps)

kr=ones(length(mesh.connectivity(:,1)),1);
pm=ones(length(mesh.connectivity(:,1)),1);

for e=1:length(mesh.connectivity(:,1))
    
    n_e = mesh.connectivity(e,:);
    
    p_loc=p_res(n_e);
    
    p_e=mean(p_loc); % taking the mean....
    
    
    kr(e)=0.5*(1.+tanh(p_e/eps))+1e-10; %+1.e-9
    pm(e)=p_e;
    
%    disp([p_e kr(e)]);
end


end
