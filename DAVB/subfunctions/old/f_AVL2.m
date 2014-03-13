function [fx] = f_AVL2(x,Theta,u,in)
% OTO associative learning evolution function
% function [fx] = f_AVL2(x,Theta,u,in)
% This function implements the VB recognition of audio-visual cue/outcomes
% pairs (it learns about the cue/outcome association).
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
%       .flag: flag for the percfeptual model (flag=1: static,
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
in.G0 = 1;
in.beta = 1;

args.u = u;
args.prior = prior;
args.Theta = Theta;
args.in = in;
fname = @gradF;
dt = 2e-1;
timeGrid = 0:dt:1e3;
% prior(1) = 0.5;
x0 = [sigm(prior(1),struct('INV',1));prior(2);log(prior(3));0];
[x,out] = integrateDynSys(fname,timeGrid,x0,struct('verbose',0),args);

pxi = sigm(x(1,:));
dpdt = diff(pxi)./dt;
Ex = x(2,:);
Vx = exp(x(3,:));

hf = figure;
subplot(2,2,1)
plot(abs(dpdt));
hold on
plot([1,length(dpdt)],[exp(-2),exp(-2)],'g');
subplot(2,2,2)
plot(x(4,:))
subplot(2,2,3)
plot(x')

ind = find(abs(dpdt)<=exp(-2), 1 );
tau = timeGrid(ind+1)

fx = [pxi(end);Ex(1,end);Vx(1,end);tau];
pause


function [dF] = gradF(t,x,args)
dF = numericDiff(@Fwrapper,1,x(1:3),args);
dF = dF(:);
dF(4) = Fwrapper(x(1:3),args);

function F = Fwrapper(x,args)
F = freeEnergy(x(1),x(2),x(3),args.u,args.prior,args.Theta,args.in);

function F = freeEnergy(pxi,Ex,Vx,u,prior,Theta,in)
try;mu=in.mu;catch;mu=[1+eps;-1+eps];end
delta1 = u - mu(1);
delta2 = u - mu(2);
pxi = sigm(pxi);
Vx = exp(Vx);
iva = exp(Theta(1));
Sqp = zeros(size(pxi));
ind = find(pxi>=1e-3&pxi<=1-1e-3);
n = length(pxi);
Sqp(setdiff(1:n,ind)) = -pxi.*log(pxi) - (1-pxi).*log(1-pxi);
l2pi = log(2*pi);
[sx] = sigm(Ex(1,:),in);
sx = sx(:)';
lsx = log(sx);
if size(Ex,1) < 2
    Ex0 = prior(2);
    Vx0 = prior(3);
    F = -0.5*iva.*(pxi.*delta1.^2+(1-pxi).*delta2.^2) ...
        + 0.5*Theta(1) - 0.5*l2pi ...
        + Ex.*(pxi-1) + lsx + 0.5.*(sx.^2-sx).*Vx ...
        + 0.5*(-Vx0.^-1.*((Ex-Ex0).^2+Vx) -l2pi -log(Vx0)) ...
        + Sqp + 0.5.*log(Vx) + 0.5.*(l2pi+1);
    
else
    Ex0 = prior([2;4]);
    Vx0 = prior([3;5]);
    ee = exp(Ex(2,:));
    V = 1./(Vx0(1)+ee);
    F = -0.5*iva.*(pxi.*delta1.^2+(1-pxi).*delta2.^2) ...
        + 0.5*Theta(1) - 0.5*l2pi ...
        + Ex(1,:).*(pxi-1) + lsx + 0.5.*(sx.^2-sx).*Vx(1,:) ... 
        - 0.5*l2pi - 0.5*(log(V) + 0.5*Vx(2,:).*Vx0(1).*ee.*V.^2) ...
        - 0.5*(V)... + 0.5*Vx(2,:).*ee.*(ee-Vx0(1)).*V.^3)...
        .*((Ex(1,:)-Ex0(1)).^2+Vx(1,:)) ...
        -0.5*l2pi - 0.5*log(Vx0(2)) - 0.5*Vx0(2).^-1 ...
        .*((Ex(2,:)-Ex0(2)).^2+Vx(2,:)) ...
        + Sqp + 0.5.*log(Vx(1,:)) + 0.5.*log(Vx(2,:)) + l2pi+1;
end

