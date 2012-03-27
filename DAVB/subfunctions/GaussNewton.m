function [opt,sigma,out] = GaussNewton(fname,init,options)
% Gauss-Newton maximization scheme
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

% Get default options
try
    args = options.args;
catch
    nargs = nargin(fname)-1;
    args = cell(nargs,1);
end
try
    maxIter = options.GnMaxIter;
catch
    maxIter = 32;
end
try
    rdI = options.GnTolFun;
catch
    rdI = 1e-12;
end
try
    oneIteration = ~~options.oneIteration;
catch
    oneIteration = 0;
end
    
out.conv = 0;
out.nReg = 0;

% Initialization
PreviousMu = init;
[I,PreviousSigma,deltaMu] = feval(fname,PreviousMu,args{:});
[deltaMu,flag] = VBA_checkGN(PreviousSigma,deltaMu);
out.nReg = out.nReg + flag;
PreviousI = I;
out.I = I;

% Main Gauss-Newton scheme
opt = init;
sigma = PreviousSigma;
stop = 0;
it = 0;
while ~stop
    it = it+1;
    % make a move
    mu = PreviousMu + deltaMu;
    % get energy function (as well as next move)
    try
        [I,NextSigma,NextdeltaMu] = feval(fname,mu,args{:});
        [NextdeltaMu,flag] = VBA_checkGN(NextSigma,NextdeltaMu);
        out.nReg = out.nReg + flag;
    catch
        I = -Inf;
    end
    % calculate relative energy function improvement
    deltaI = I-PreviousI;
    % check whether to stop, to accept move or to halve step
    if it <= maxIter && abs(deltaI./PreviousI)>rdI
        if deltaI<0     % halve step size
            deltaMu = 0.5*deltaMu;
        else            % accept move
            % 1- propose a new move according to new local quadratic approx
            deltaMu = NextdeltaMu;
            % 2- update sufficient stats
            PreviousMu = mu;
            PreviousI = I;
            % 3- build posterior and sufficient statistics structures
            opt = mu;
            sigma = NextSigma;
            stop = oneIteration;
            out.I = [out.I,I];
        end
    elseif abs(deltaI./PreviousI)<rdI
        % stop Gauss-Newton search
        stop = 1;
        out.conv = 1;
    elseif it > maxIter
        stop = 1;
    end
end

% Optional output
out.it = it;


