function [kr]=ComputeRelPerm(mesh,p_res,eps)

kr=ones(length(mesh.connectivity(:,1)),1);


for e=1:length(mesh.connectivity(:,1))
    
    n_e = mesh.connectivity(e,:);
    
    p_loc=p_res(n_e);
    
    p_e=mean(p_loc);
    
    kr(e)= (1.+tanh(p_e/eps))/2.; 
    
end


end
