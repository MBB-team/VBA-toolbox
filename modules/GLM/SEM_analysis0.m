function [estimates,posterior,out] = SEM_analysis0(Y,A,B,C,flag)
% structural equation modelling (bayesian stats)
% FORMAT: function [estimates,posterior,out] = SEM_analysis0(Y,A,B,C)
% This function assumes a GLM of the following form:
% Y' = A*Y' + sum_i B{i}*Y(i)*Y' + C + e
% where e are normal iid residuals.
% IN:
%   - Y: nxk data matrix, where n are the number of samples per data
%   dimension (and there are k dimensions)
%   - A: kxk "connectivity" matrix
%   - B: kx1 cell array of kxk modulatory matrices
%   - C: kx1 vector of constant terms
%   - flag: flag for enforcing hard constraints (on all path coeff):
%       flag='unbounded': no constraint {default}
%       flag='exp': positivity constraint
%       flag='sigm': constraint on the [0,1] interval
% NB: only non-zero entries of matrices A, B and C are considered, i.e. SEM
% fit proceeds by estimating unknown path weights that correspond to these
% non-zero entries.
% OUT:
%   - estimates: a structure with the following fields:
%       .mean: a structure containing the parameter estimates of the path
%       weigths that correspond to non-zero entries of A, B and C matrices.
%       .var: a structure containing the posterior variance over path
%       weights.
%   - posterior/out: VBA posterior structures.

try
    switch flag
        case 'unbounded'
        case 'exp'
        case 'sigm'
        otherwise
            disp('SEM_analysis0: unknown ''flag'' input. Defaulting to ''unbounded''.')
            flag='unbounded';
    end
catch
    disp('SEM_analysis0: unspecified ''flag'' input. Defaulting to ''unbounded''.')
    flag='unbounded';
end


% initialize dummy variables
[n,k] = size(Y);
in.A = A;
in.indInA = find(A~=0);
in.indA = 1:length(in.indInA);
in.B = B;
np = length(in.indInA);
for i=1:k
    in.indInB{i} = find(B{i}~=0);
    in.indB{i} = np+1:np+length(in.indInB{i});
    np = np+length(in.indInB{i});
end
in.C = C;
in.indInC = find(C~=0);
in.indC = np+1:np+length(in.indInC);
in.flag = flag;

% deal with NaN and/or missing data
inan = find(isnan(Y));
Y(inan) = 0;
if ~isempty(inan)
    options.isYout = zeros(size(Y));
    options.isYout(inan) = 1;
    options.isYout = options.isYout';
end

% call VBA inversion routine
y = Y';
u = Y';
g_fname = @g_SEM;
dim.n_phi = np + length(in.indInC);
dim.n = 0;
dim.n_theta = 0;
options.inG = in;
options.priors.muPhi = zeros(dim.n_phi,1);
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.DisplayWin = 0;
options.verbose = 1;
[posterior,out] = VBA_hyperparameters(y,u,[],g_fname,dim,options);

% wrap-up SEM estimates
mP = posterior.muPhi;
vP = diag(posterior.SigmaPhi);
if dim.n_phi>0
    switch flag
        case 'unbounded'
        case 'exp'
            mP = exp(mP);
            vP = vP.*exp(2*mP); % Laplace approx
        case 'sigm'
            mP = VBA_sigmoid(mP);
            vP = vP.*(mP.*(1-mP)).^2; % Laplace approx
    end
end
estimates.mean.A = in.A;
estimates.mean.A(in.indInA) = mP(in.indA);
estimates.var.A = in.A;
estimates.var.A(in.indInA) = vP(in.indA);
estimates.mean.C = in.C;
estimates.mean.C(in.indInC) = mP(in.indC);
estimates.var.C = in.C;
estimates.var.C(in.indInC) = vP(in.indC);
for i=1:k
    estimates.mean.B{i} = in.B{i};
    estimates.mean.B{i}(in.indInB{i}) = mP(in.indB{i});
    estimates.var.B{i} = in.B{i};
    estimates.var.B{i}(in.indInB{i}) = vP(in.indB{i});
end


function [gx] = g_SEM(x,P,u,in)
switch in.flag
    case 'unbounded'
    case 'exp'
        P = exp(P);
    case 'sigm'
        P = VBA_sigmoid(P);
end
A = in.A;
A(in.indInA) = P(in.indA);
C = in.C;
C(in.indInC) = P(in.indC);
gx = A*u + C;
for i=1:size(u,1)
    B = in.B{i};
    B(in.indInB{i}) = P(in.indB{i});
    gx = gx + B*u*u(i);
end