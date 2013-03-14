function [gx,dgdx,dgdp] = g_demo_multisource(Xt,P,ut,in) 

%- initializations
idX1 = 1:in.n5(end);
nreg=length(in.n5);
nresp=length(in.r);

gxr=zeros(nresp,1);
dgdxr=zeros(length(Xt),nresp);
dgdpr=zeros(length(P),nresp);

%== compute gradients
[gxn,dgdxn,dgdpn] = g_HRF3(Xt(idX1),P,ut,in) ;
[gxr(1),dgdxr(in.r(1),1),dgdpr(:,1)] = g_softmax4decoding(Xt(in.r(1)),P,ut,in);
[gxr(2),dgdxr(in.r(2),2),dgdpr(:,2)] = g_softmax4decoding(Xt(in.r(2)),P,ut,in);

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