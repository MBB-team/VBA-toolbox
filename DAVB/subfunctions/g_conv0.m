function [gx,dgdx,dgdp] = g_conv0(x,P,u,in)

gx = P(end);
nt = size(u,1)./in.dim.nu;
dgdx = [];
dgdp = zeros(size(P,1),nt);
for i=1:in.dim.nu
    K = zeros(nt,nt);
    Pi = [P((i-1)*in.dim.n_t+1:i*in.dim.n_t);zeros(nt-in.dim.n_t,1)];
    ui = u((i-1)*nt+1:i*nt);
    for j=1:nt
        K(:,j) = circshift(Pi,j-1);
        K(1:j-1,j) = 0;
        if j<=in.dim.n_t
            dgdp((i-1)*in.dim.n_t+j,:) = circshift(ui',[0,j-1]);
            dgdp((i-1)*in.dim.n_t+j,1:j-1) = 0;
        end
    end

    gx = gx + K*ui;
end
dgdp(end,:) = 1;
