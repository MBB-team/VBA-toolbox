function [EUi,Ui,Z,d2Uidx2,d2UidP2,Vy] = VBA_evalAL2(Xt,P,beta,ut,yt,options)
% smart wrapper for unnormalized log-likelihoods
% function [Eui,Ui,Zi,d2Uidx2,d2UidP2] = VBA_evalAL2(Xt,P,beta,ut,yt,options)
%
% This function evaluates the re-normalized log-likelihood at the observed
% data point. This is done by computing the log- partition function, using
% numerical integration. The ensuing gradients and Hessians are derived
% from gradients of the unnormalized log-likelihood.
%
% IN:
%   - Xt: the hidden states
%   - P: the parameters
%   - ut: the inputs to the system
%   - yt: the observations, at which the arbitrary likelihood is evaluated
%   - options: the options structure (it contains the name of the function
%   that evaluates the unnormalized log-likelihood).

% 1- evaluate Hessians at y
[Ui,dUidx,dUidp,d2Uidx2,d2UidP2] = evalUNL(Xt,P,ut,yt,options);
if options.dim.n == 0
    d2Uidx2 = zeros(options.dim.n,options.dim.n);
end
if options.dim.n_phi == 0
    d2UidP2 = zeros(options.dim.n_phi,options.dim.n_phi);
end

% 2- evaluate partition function
sigma = options.UNL_width/sqrt(beta); % to be changed for adaptive evaluation {1e1}
dy = sigma./options.UNL_ng;
gy = dy:dy:sigma;
gy = [yt-fliplr(gy),yt,yt+gy];
ng = length(gy);
Ux = zeros(ng,1);
for j=1:ng
    Ux(j) = feval(options.g_fname,Xt,P,ut,gy(j),options.inG);
end
mU = max(beta*Ux);
eU = exp(beta*Ux-mU);
Z = log(sum(eU*dy))+mU;
EUi = dy.*exp(mU-Z).*Ux'*eU;
Ey = dy.*exp(mU-Z).*gy*eU;
Vy = dy.*exp(mU-Z).*(gy-Ey).^2*eU;
% figure,plot(Ux)
% keyboard

function [Ux,dUdx,dUdp,d2Udx2,d2Udp2] = evalUNL(Xt,P,ut,yt,options)
deriv = [1 1 1 1 1];
nout = options.g_nout;
g_fname = options.g_fname;
in = options.inG;
dim = options.dim;
switch nout
    case 5
        [Ux,dUdx,dUdp,d2Udx2,d2Udp2] = feval(g_fname,Xt,P,ut,yt,in);
        if isempty(dUdx)
            deriv(1) = 0;
        end
        if isempty(dUdp)
            deriv(2) = 0;
        end
        if isempty(d2Udx2)
            deriv(3) = 0;
        end
        if isempty(d2Udp2)
            deriv(4) = 0;
        end
    case 4
        [Ux,dUdx,dUdp,d2Udx2] = feval(g_fname,Xt,P,ut,yt,in);
        deriv(4) = 0;
        if isempty(dUdx)
            deriv(1) = 0;
        end
        if isempty(dUdp)
            deriv(2) = 0;
        end
        if isempty(d2Udx2)
            deriv(3) = 0;
        end
    case 3
        [Ux,dUdx,dUdp] = feval(g_fname,Xt,P,ut,yt,in);
        deriv(3:4) = 0;
        if isempty(dUdx)
            deriv(1) = 0;
        end
        if isempty(dUdp)
            deriv(2) = 0;
        end
    case 2
        [Ux,dUdx] = feval(g_fname,Xt,P,ut,yt,in);
        deriv(2:4) = 0;
        if isempty(dUdx)
            deriv(1) = 0;
        end
    case 1
        [Ux] = feval(g_fname,Xt,P,ut,yt,in);
        deriv(1:4) = 0;
end
if ~deriv(1)
    if dim.n==0
        dUdx = zeros(1,dim.n);
    else
        if deriv(3)
            dUdx = VBA_numericDiff(g_fname,1,Xt,P,ut,yt,in);
        else
            [d2Udx2,dUdx] = VBA_numericDiff(@numericDiff,3,g_fname,1,Xt,P,ut,yt,in);
            deriv(3) = 1;
        end
    end
end
if ~deriv(3)
    if dim.n==0
        d2Udx2 = zeros(dim.n,dim.n);
    else
        d2Udx2 = VBA_numericDiff(@getDU,3,g_fname,2,Xt,P,ut,yt,in);
    end
end
if ~deriv(2)
    if dim.n_phi==0
        dUdx = zeros(1,dim.n);
    else
        if deriv(4)
            dUdp = VBA_numericDiff(g_fname,2,Xt,P,ut,yt,in);
        else
            [d2Udp2,dUdp] = VBA_numericDiff(@numericDiff,4,g_fname,2,Xt,P,ut,yt,in);
            deriv(4) = 1;
        end
    end
end
if ~deriv(4)
    if dim.n_phi==0
        d2Udx2 = zeros(dim.n_phi,dim.n_phi);
    else
        d2Udp2 = VBA_numericDiff(@getDU,4,g_fname,3,Xt,P,ut,yt,in);
    end
end

function dU = getDU(g_fname,ind,Xt,P,ut,yt,in)
[Ux,dUdx,dUdp] = feval(g_fname,Xt,P,ut,yt,in);
if ind==2
    dU = dUdX;
elseif ind==3
    dU = dUdp;
end
