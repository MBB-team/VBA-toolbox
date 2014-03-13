function [gx,dgdx,dgdp] = g_AVL(x,Phi,u,in)
% reaction time OTO observation function
% function [gx,dgdx,dgdp] = g_AVL(x,Phi,u,in)
% This function computes the observation function (reaction time) for the
% audio-visual associative learning task.
% IN:
%   - x: the 5x1 vector of states. The third first entries are sufficient
%   statistics of the updated belief. The two last ones are summary
%   statistics of the convergence rate of the VB updates.
%   - Phi: 2x1 vector of observation parameters. Phi(1) is an intercept
%   constant and Phi(2) is a scaling factor.
%   - u: useless.
%   - in: structure for flagging to different cost functions
% OUT:
%   - gx: expected reaction time
%   - dgdx: derivative wrt states

try
    flag = in.flag;
catch
    flag = 'linear';
end

switch flag
    
    case 'linear'
        
        dgdx = zeros(size(x,1),1);
        dgdp = zeros(size(Phi,1),1);
        if ~isinf(x(end))
            gx = exp(-Phi(1)).*(Phi(1) + x(end) -2.*Phi(2));
            dgdx(end) = exp(-Phi(1));
            dgdp(1) = exp(-Phi(1)) - gx;
            dgdp(2) = -2*exp(-Phi(1));
        else
            gx = 0;
        end

        
    case 'polynomial'
        
%         tic
%         opt.args = {exp(x(end)),Phi,@pol};
%         [mu,curv,out] = optimCost(@ODT,0,opt);
%         gx = exp(mu);
%         toc
        
%         tic
%         options.GnMaxIter = 256;
%         options.args = {exp(x(end)),Phi,@pol};
%         [mu2,curv2,out2] = GaussNewton(@ODT,0,options);
%         toc
% %         gx = exp(mu);
%         

        fname = @pol;
%         gx = get_gx(x,Phi,fname);
% 
% 
% 
%         dF = exp(x(end));
%         alpha = exp(Phi(1));
%         [dc,d2c,d3c,dcdp] = pol(gx,Phi(2:end));
%         nt = (alpha*dc + 2*d2c.^2).^-1;
        
%         dcdp2 = numericDiff(@pol,2,gx,Phi(2:end))
%         dcdp
%         pause
        
        
%         dgdx = zeros(size(x,1),1);
%         dgdx(end) = dc*nt./gx;
%         
%         dgdp = zeros(size(Phi,1),1);
%         dgdp(1) = dc*(alpha.^-1-gx).*nt./gx;
%         dgdp(2:end) = -d2c*dcdp.*nt./gx;
        
                
        [dgdx,gx] = numericDiff(@get_gx,1,x,Phi,fname);
        
        [dgdx2] = numericDiff_old(@get_gx,1,x,Phi,fname);
        
%         dgdx
%         
        [dgdp,gx] = numericDiff(@get_gx,2,x,Phi,fname);
        [dgdp2] = numericDiff_old(@get_gx,2,x,Phi,fname);
        
        
        dgdx
        dgdx2
        
        dgdp
        dgdp2
        pause
        
%         dgdp
%         
%         gx
%         pause
        
        
        
% 
% 
%         
%         A = log(alpha.*dF);
%         exp(mu)
% %         exp(mu2)
%         tau = A./alpha - 2*log(feval(@pol,exp(mu),Phi(2:end)))./alpha
%         pause
        
%         dgdx = NaN;
%         dgdp = NaN;


    case 'postponed'
        
        opt.args = {exp(x(end)),Phi,@post};
        [mu,curv,out] = optimCost(@ODT,0,opt);
        gx = exp(mu);
        
        
        
        dF = exp(x(end));
        alpha = exp(Phi(1));
        [dc,d2c,d3c,dcdp] = post(exp(mu),Phi(2:end));
        
        nt = (alpha*dc + 2*d2c).^-1;
        
        dgdx = zeros(size(x,1),1);
        dgdx(end) = 0.5*gx.*nt;
        
        dgdp = zeros(size(Phi,1),1);
        dgdp(1) = 0.5*gx.*(alpha.^-1-gx).*nt;
        dgdp(2:end) = -gx.*dcdp.*nt;


        
        
%         
%         options.GnMaxIter = 256;
%         options.args = {exp(x(end)),Phi,@post};
%         [mu,curv,out] = GaussNewton(@ODT,0,options);
%         gx = exp(mu);
        
%         dF = exp(x(end));
%         alpha = exp(Phi(1));
%         A = log(alpha.*dF);
%         exp(mu)
%         tau = A./alpha - 2*log(feval(@post,exp(mu),Phi(2:end)))./alpha
%         pause
        
        
end

function gx = get_gx(x,Phi,fname)
opt.args = {exp(x(end)),Phi,fname};
[mu,curv,out] = optimCost(@ODT,0,opt);
gx = exp(mu);

function [I,S,dx] = ODT(x,dF,Phi,fname)
alpha = exp(Phi(1));
A = log(alpha.*dF);
t = exp(x);
[y,dy,d2y] = feval(fname,t,Phi(2:end));
I = -(alpha.*t + 2*log(y) -A).^2;
dI = -(2*alpha.^2*t - 2*alpha*A + 4*alpha*log(y) + 4*alpha*t*dy./y ...
    + 8*log(y).*dy./y) ...
    *t;
d2I = dI - (...
    2*alpha.^2 + 8*alpha*dy./y ...
    + 4*alpha*t*(d2y*y-dy.^2)./y.^2 ...
    + 8*(dy.^2./y.^2 + log(y).*(d2y*y-dy.^2)./y.^2) ) ...
    *t;
S = -d2I.^-1;
dx = S*dI;

function [y,dy,d2y,dydp] = pol(x,P)
y = exp(P(1)).*exp(P(2)).*x.^(exp(P(2))-1);
dy = exp(P(1)).*exp(P(2)).*(exp(P(2))-1).*x.^(exp(P(2))-2);
d2y = exp(P(1)).*exp(P(2)).*(exp(P(2))-1).*(exp(P(2))-2).*x.^(exp(P(2))-3);
dydp = [y;y*(1+log(x).*exp(P(2)))];

function [y,dy,d2y,dydp] = post(x,P)
in.G0 = 1;
in.beta = 1;
[y,dy,dydp] = sigm(x,in,P);
d2y = 1;







