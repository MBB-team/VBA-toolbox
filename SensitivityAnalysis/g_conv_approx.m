function [gx,dgdx,dgdp] = g_conv_approx(~,P,u,in)

nt = round(size(u,1)./in.dim.nu);

gx = P(end)*ones(nt,1);
dgdx = [];
dgdp = zeros(size(P,1),nt);


for i=1:in.dim.nu
    
    idx = (i-1)*3+1:3*i;
    Pi = P(idx);
    [kernel,dkdp] = kernel_sinexp(Pi(1),Pi(2),Pi(3),(0:in.dim.n_t-1)*in.deltat);
    ui = u((i-1)*nt+1:i*nt)';
    
    gx = gx + myconv(ui,kernel)';
    
    dgdp(idx,:) = [myconv(ui, dkdp(1,:)); myconv(ui, dkdp(2,:)) ; myconv(ui, dkdp(3,:))];
        

end
dgdp(end,:) = 1;

function c = myconv(u,k)
    c=conv(u,k,'full');
    c(end-length(k)+2:end) = [];
end

end


