function [fx] = f_AVL(x,Theta,u,in)
% OTO: associative learning evolution function
% function [fx] = f_AVL(x,Theta,u,in)
% This function implements the VB recognition of audio-visual cue/outcomes
% pairs (it learns about the cue/outcome association). The recognition
% process expects:
%   - sensory cue (u). This is a vector of 'sensory data', presented to the
%   subject, that is compared with two categories, which are defined
%   through their 'average cue' (given in in.mu(:,1:2)). This comparison
%   yields a posterior belief about which category the sensory cue belongs
%   to. NB: the probability of the sensory cue belonging to the first
%   category is the first entry in the state vector, x(1).
%   - parameters of the perceptual model (Theta). These are the variability
%   of sensory cues within each category (Theta(1) is the class
%   log-precision).
% IN:
%   - x: only the 2d and 3d entries of the x vector are used. They contain
%   the prior sufficient statistics of the unknown variables to be
%   recognised, i.e. the visual outcome label (face/house) and its
%   frequency (prior probability, i.e. auditory-visual association).
%   - Theta: 2 x 1 vector containing the required parameters of the
%   recognition model: Theta(1) is the log-precision of the sensory signals
%   and Theta(2) is the log-volatility of the association frequency.
%   - u: scalar containing the visual outcome "signal", which is then
%   compared to "template houses" and "template face", in order to be
%   classified as such.
%   - in: optional input structure. Can contain the following fields:
%       .flag: flag for the perceptual model (flag=1: static,
%       flag=2:dynamic, flag=3: volatile)
%       .n: the max number of VB updates for recognition {32}
%       .mu: 2x1 vector containing the "template house/face" to which the
%       visual outcome signal (u) is compared {[-1;1]}.
% OUT:
%   - fx: in.n x 1 vector containing the series of VB updates of the
%   sufficient statistics of the recognised variables.

% initialize current belief with top-down prior predictions from highest
% hierarchical level
switch in.flag
    case 1 % static model
        prior = x(1:3);
    case 2 % dynamic model
        prior = x(1:3);
        prior(3) = prior(3) + exp(Theta(2));
    case 3 % volatile model
        prior = x(1:5);
        prior(5) = prior(5) + exp(Theta(2));
end

% VB update of the sufficient statistics (in.n iterations)
[pxi,Ex,Vx] = VBA_AVL(prior,Theta,u(in.uu),in);

switch in.flag
    case {1,2}
        fx = [pxi(end);Ex(1,end);Vx(1,end);pxi(1)];
    case 3
        fx = [pxi(end);Ex(1,end);Vx(1,end);Ex(2,end);Vx(2,end);pxi(1)];
end

% posterior = [pxi(end);Ex(end);Vx(end)];
% if posterior(1) > 0.99
%     posterior(1) = 0.99;
% end
% if posterior(1) <0.01
%     posterior(1) = 0.01;
% end
% % posterior(1) = sigm(0.5*posterior(1)+0.25,struct('INV',1))
% % d2f = numericDiff(@numericDiff,3,@FE,1,posterior,u,prior,Theta,in);
% [F,d2f,d2f0] = FE(posterior,u,prior,Theta,in);
% % [hf] = VBA_displayGrads(d2f,d2f0,'Jacobian','FE','g');
% % pause
% % [F,d2f,d2f0] = FE(posterior,u,prior,Theta,in);
% % [hf] = VBA_displayGrads(d2f,pinv(d2f),'Jacobian','FE','g');
% J = - pinv(d2f)*d2f0';
% dfdx = [    J                       , zeros(3,1)
%             [0,pxi(1)*(1-pxi(1)),0] , 0  ];
% dfdx = dfdx';


function [pxi,Ex,Vx] = VBA_AVL(prior,Theta,u,in)

% allocate dummy variables
opt.GnTolFun = 1e-4;
opt.iffig = 0;
opt.oneIteration = 1;
opt.GnMaxIter = 8;
mu=in.mu;
n=in.n;
tdf=in.tdf;
% initialize VB sufficient statistics
pxi = zeros(1,n+1);
pxi(1) = VBA_finiteBinomial (1./(1+exp(-prior(2))));
switch in.flag
    case {1,2}
        Ex = zeros(1,n+1);
        Vx = zeros(1,n+1);
        Ex(1) = prior(2);
        Vx(1) = prior(3);
    case 3
        Ex = zeros(2,n+1);
        Vx = zeros(2,n+1);
        Ex(:,1) = [prior(2);prior(4)];
        Vx(:,1) = [prior(3)+ exp(prior(4));prior(5)];
end

% allocate likelihood variables
delta1 = sum((u - mu(:,1)).^2);
delta2 = sum((u - mu(:,2)).^2);
iva = exp(Theta(1));
l2pi = log(2*pi);

% iterative VB algorithm from prediction error (lowest hierarchical level)
F = zeros(n+1,1);
F(1) = freeEnergy(pxi(1),Ex(:,1),Vx(:,1),u,prior,Theta,in);
for i = 1:n
    % update xi
    [sx] = VBA_finiteBinomial (1./(1+exp(-Ex(i))));
    lsx = log(sx);
    Elsx = lsx + 0.5.*Vx(i).^2.*(sx.^2-sx);
    p(1) = -0.5.*(iva.*delta1 + Elsx +l2pi - Theta(1));
    p(2) = -0.5.*(iva.*delta2 + Elsx-Ex(i) +l2pi - Theta(1));
    p = exp(p-max(p));
    pxi(i+1) = VBA_finiteBinomial (p(1)./sum(p));
    switch in.flag
        case {1,2}
            % update x1
            opt.args = {prior(2),prior(3),pxi(i+1)};
            [out1,out2] = VBA_GaussNewton('expBinom',Ex(i),opt);
            if ~isempty(out1)
                Ex(i+1) = out1;
                Vx(i+1) = out2;
            else % VB convergence
                Ex(i+1:end) = Ex(i);
                Vx(i+1:end) = Vx(i);
                pxi(i+1:end) = pxi(i+1);
                break
            end
        case 3
            % update x1
            ee2 = exp(Ex(2,i));
            tmp = prior(3) + ee2;
            mft = 0;%max([0.5.*Vx(2,i).*ee2.*(ee2-prior(3))./tmp.^3,0]);
            P0 = 1./tmp + mft;
            opt.args = {prior(2),1./P0,pxi(i+1)};
            [out1,out2] = VBA_GaussNewton('expBinom',Ex(1,i),opt);
            if ~isempty(out1)
                Ex(1,i+1) = out1;
                Vx(1,i+1) = out2;
            else % VB convergence
                Ex(1,i+1:end) = Ex(1,i);
                Vx(1,i+1:end) = Vx(1,i);
                Ex(2,i+1:end) = Ex(2,i);
                Vx(2,i+1:end) = Vx(2,i);
                pxi(i+1:end) = pxi(i+1);
                break
            end
            % update x2
            a = prior(3);
            b = (prior(2)-Ex(1,i+1)).^2 + Vx(1,i+1);
            opt.args = {prior(4),prior(5),a,b};
            o.conv = 0;
            init = Ex(2,i);
            while ~o.conv
                [out1,out2,o] = VBA_GaussNewton('VarVolatility',init,opt);
                init = out1;
            end
            if ~isempty(out1)
                Ex(2,i+1) = out1;
                Vx(2,i+1) = out2;
            else % VB convergence
                Ex(2,i+1:end) = Ex(2,i);
                Vx(2,i+1:end) = Vx(2,i);
                Ex(1,i+1:end) = Ex(1,i+1);
                Vx(1,i+1:end) = Vx(1,i+1);
                pxi(i+1:end) = pxi(i+1);
                break
            end

    end
    % check VB convergence
    F(i+1) = freeEnergy(pxi(i+1),Ex(:,i+1),Vx(:,i+1),u,prior,Theta,in);
    dF = F(i+1) - F(i);
    if abs(dF) <= tdf
        switch in.flag
            case {1,2}
                Ex(i+1:end) = Ex(i+1);
                Vx(i+1:end) = Vx(i+1);
                pxi(i+1:end) = pxi(i+1);
            case 3
                Ex(2,i+1:end) = Ex(2,i+1);
                Vx(2,i+1:end) = Vx(2,i+1);
                Ex(1,i+1:end) = Ex(1,i+1);
                Vx(1,i+1:end) = Vx(1,i+1);
                pxi(i+1:end) = pxi(i+1);
        end
        break
    end
end

function [F,d2f,d2f0] = FE(posterior,u,prior,Theta,in)
% posterior(1) = sigm(posterior(1)); % for numerical derivation stability
[F,d2f,d2f0] = freeEnergy(...
    posterior(1),posterior(2),posterior(3),u,prior,Theta,in);

function [F,d2f,d2f0] = freeEnergy(pxi,Ex,Vx,u,prior,Theta,in)
try;mu=in.mu;catch;mu=[1+eps,-1];end
delta1 = sum((u - mu(:,1)).^2);
delta2 = sum((u - mu(:,2)).^2);
iva = exp(Theta(1));
Sqp = -pxi.*log(pxi) - (1-pxi).*log(1-pxi);
l2pi = log(2*pi);
[sx] = VBA_finiteBinomial (1./(1+exp(-Ex(1,:))));
sx = sx(:)';
lsx = log(sx);
if size(Ex,1) < 2
    Ex0 = prior(2);
    Vx0 = prior(3);
    F = -0.5*iva.*(pxi.*delta1+(1-pxi).*delta2) ...
        + 0.5*Theta(1) - 0.5*l2pi ...
        + Ex.*(pxi-1) + lsx + 0.5.*(sx.^2-sx).*Vx ...
        + 0.5*(-Vx0.^-1.*((Ex-Ex0).^2+Vx) -l2pi -log(Vx0)) ...
        + Sqp + 0.5.*log(Vx) + 0.5.*(l2pi+1);

else
    Ex0 = prior([2;4]);
    Vx0 = prior([3;5]);
    ee = exp(Ex(2,:));
    V = 1./(Vx0(1)+ee);
    F = -0.5*iva.*(pxi.*delta1+(1-pxi).*delta2) ...
        + 0.5*Theta(1) - 0.5*l2pi ...
        + Ex(1,:).*(pxi-1) + lsx + 0.5.*(sx.^2-sx).*Vx(1,:) ...
        - 0.5*l2pi - 0.5*(log(V) + 0.5*Vx(2,:).*Vx0(1).*ee.*V.^2) ...
        - 0.5*(V)... + 0.5*Vx(2,:).*ee.*(ee-Vx0(1)).*V.^3)...
        .*((Ex(1,:)-Ex0(1)).^2+Vx(1,:)) ...
        -0.5*l2pi - 0.5*log(Vx0(2)) - 0.5*Vx0(2).^-1 ...
        .*((Ex(2,:)-Ex0(2)).^2+Vx(2,:)) ...
        + Sqp + 0.5.*log(Vx(1,:)) + 0.5.*log(Vx(2,:)) + l2pi+1;
end

if nargout < 2 ; return; end

ds = sx*(1-sx);
d2s = sx*(1-sx)*(1-2*sx);
d3s = sx*(1-sx)*(1-2*sx)*(1-3*sx);
d2f = [ -1/(pxi*(1-pxi))    ,   1                   ,   0
        1                   , -ds-0.5*Vx*d3s-Vx0    , -0.5*d2s
        0                   , -0.5*d2s              , -0.5/Vx^2 ];

d2f0 = [ 0    ,   0      ,   0
         0    , 1/Vx0    , (Ex-Ex0)/Vx0
         0    , 0        , 0.5/Vx0^2     ];
        



