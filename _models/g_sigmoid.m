function [gx,dG_dX,dG_dPhi] = g_sigmoid(Xt,Phi,ut,inG)
% partially observable sigmoid mapping

n = size(Xt,1);
try
    ind1 = inG.ind;
catch
    ind1 = 1:n;
end
G = eye(n);
G = G(ind1,:);

try
    G0 = G;
    G = Phi(1).*G;
end

if size(Phi,1) >=2
    [Sx,dsdx,dsdp] = sigm(Xt,inG,Phi(2:end));
else
    [Sx,dsdx] = VBA_sigmoid(Xt,inG);
end

gx              = G*Sx(:);
dG_dX           = [G*diag(dsdx(:))]';

dG_dPhi = zeros(size(Phi,1),length(ind1));
if size(Phi,1) >=1
    dG_dPhi(1,:)   = [G0*Sx(:)]';
end
if size(Phi,1) >=2
    dG_dPhi(2:end,:)     = dsdp(:,ind1);
end


% [dsdp2] = numericDiff(@sigm,3,Xt,inG,Phi(2:end));
% 
% 
% hf=figure,
% subplot(2,2,1),imagesc(dsdp),colorbar,title('dsdp')
% subplot(2,2,2),imagesc(dsdp2),colorbar,title('dsdp2')
% subplot(2,2,3),imagesc(dsdp-dsdp2),colorbar,title('dsdp-dsdp2')
% subplot(2,2,4),plot(dsdp(:),dsdp2(:),'.')
% 
% grid on
% disp('pause')
% pause
% close(hf)


