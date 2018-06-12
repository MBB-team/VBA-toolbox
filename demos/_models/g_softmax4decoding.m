function [gx,dgdx,dgdp] = g_softmax4decoding(Xt,P,ut,in)


eta=P;
Xt = Xt+eta;
beta=100;

sx = length(Xt);
dgdx=zeros(sx,length(Xt));
dgdp=zeros(length(P),length(Xt));

% if size(Xt,1) >1
%     
%     % choice
%     ex = exp(beta*Xt);
%     gx=ex/sum(ex);
%     
%     % derivative wrt state
%     dgdx(:,:)=diag(beta*gx);
%     dgdx(:,:) = dgdx(:,:) - beta*(gx*gx');
%     
%     %derivative wrt parameters
%     dgdp = beta*(Xt.*gx - gx*sum(gx.*Xt))';
%     
% else
    % choice
    ex = exp(-beta.*Xt);
    gx=1./(1+ex);
    
    % derivative wrt state
    dgdx=diag(beta.*ex ./(ex+1).^2); 
    
    %derivative wrt parameters
    dgdp = diag(beta.*ex./(ex+1).^2);   
dgdp = diag(dgdp);
% end