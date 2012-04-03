function gx = g_nsess1(x,P,u,in)


gx = zeros(in.dim_output,1);
nsess = in.nsess;
nx_ps = size(x,1)/in.nsess;
nu_ps = size(u,1)/in.nsess;

for i=1:nsess

    [gi] = feval(...
        in.g_fname,... % the function to be called for session i
        x(nx_ps*(i-1)+1:nx_ps*i),... % the indices of the hidden states concerned by the evolution function for session i
        P(in.ind_phi(i,:)),... % the indices of the evolution parameters for session i
        u(nu_ps*(i-1)+1:nu_ps*i),... % the indices of the data for session i
        []); 
    
    gx(in.dim_output_ps*(i-1)+1:in.dim_output_ps*i) = gi;
    
end


