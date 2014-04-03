function [gx,dgdx,dgdp] = g_demo_extended(Xt,P,ut,in) 
% canonical observation function for extended DCM.
% The output is the concatenation of:
% - classical HRF applied on neural states
% - multinomial prediction from hidden extended states

%- initializations
idX1 = 1:in.n5(end);
nreg=length(in.n5);
nresp=length(in.r);


gxr=zeros(nresp,1);
dgdxr=zeros(length(Xt),nresp);
dgdpr=zeros(length(P),nresp);

%== compute gradients
[gxn,dgdxn,dgdpn] = g_HRF3(Xt(idX1),P,ut,in) ;
if nresp>0
[gxr,dgdxr(in.r,:),dgdpr] = g_softmax4decoding(Xt(in.r),P,ut,in);
else
    gxr=[];
    dgdxr=[];
    dgdpr=[];
end
%== concatenate BOLD and behavioral responses
%- state
gx = [gxn ; gxr] ;

%- Jacobian
dgdx = zeros(length(Xt),length(gx));
dgdx(idX1,1:nreg) = dgdxn ;
dgdx(:,nreg+(1:nresp)) = dgdxr ;

%- gradient wrt parameters
dgdp = zeros(length(P),nreg+nresp);
dgdp(1:size(dgdpn,1),1:nreg) = dgdpn ;
dgdp(:,nreg+(1:nresp)) = dgdpr ;


end