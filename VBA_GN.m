function [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,flag)
% regularized Gauss-Newton optimization of variational energies
% function [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,flag)
%
% This function computes a gauss-Newton scheme on the variational energy of
% either hidden states, initial conditions, observation parameters or
% evolution parameters. The algorithm is the following  modified
% Gauss-Newton scheme:
%
%   Initialize lastMU, lastI, DMU
%   While dI > epsilon
%       MU = lastMU + DMU
%       [I,nextDMU] = GaussNewton(MU)
%       dI = I- lastI
%       if dI < 0
%           DMU = DMU./2
%       else
%           DMU = nextDMU
%           lastMU = MU
%           lastI = I
%           store(MU)
%       end
%   end
%
% where GaussNewton() is a function that implements the standard
% Gauss-Newton scheme for optimizing the variational energy and store() is
% an operation that stores the 'accepted' Newton update and the inverse
% curvature of the variational energy (Laplace approximation).
% NB: if options.gradF = 1, then this routine ensures that the free energy
% (as opposed to the variational energy) increases between two
% iterations...

switch flag
    case 'X'
        if numel(options.sources)>1 || options.sources(1).type==2
            error('*** Stochastic multichannel VB is not yet supported !');
        end
        indIn = options.params2update.x;
        PreviousMu = posterior.muX;
        if ~options.binomial
            fname = @VBA_IX_lagged;
        else
            fname = @VBA_IX_lagged_binomial;
        end
        s1 = 'I(<X>) =';
        s2 = '<dX>';
    case 'X0'
        indIn = options.params2update.x0;
        PreviousMu = posterior.muX0(indIn);
        fname = @VBA_IX0;
        s1 = 'I(<X0>) =';
        s2 = '<dX0>';
    case 'Phi'
        indIn = options.params2update.phi;
        PreviousMu = posterior.muPhi(indIn);
        if options.UNL % to be rationalized...
            fname = @VBA_Iphi_UNL;
        else
            if  numel(options.sources)>1 || options.sources(1).type==2
                fname = @VBA_Iphi_extended;
            else
                if options.nmog > 1
                    if options.extended
                        error('*** Splitted multichannel VB is not yet supported !');
                    end
                    fname = @VBA_Iphi_split;
                elseif options.binomial
                    fname = @VBA_Iphi_binomial;
                else
                    fname = @VBA_Iphi;
                end
            end
        end
        s1 = 'I(<Phi>) =';
        s2 = '<dPhi>';
    case 'Theta'
        indIn = options.params2update.theta;
        PreviousMu = posterior.muTheta(indIn);
        fname = @VBA_Itheta;
        s1 = 'I(<Theta>) =';
        s2 = '<dTheta>';
end

if isempty(indIn)
    return
end

% Get variational energy (I) and propose move (deltaMu)
try
   [I,Sigma,deltaMu,suffStat2] = feval(fname,PreviousMu,y,posterior,suffStat,dim,u,options);
    PreviousI = I;
catch
    VBA_disp(['Warning: could not evaluate variational energy on ',flag,'!'],options)
    return
end

% Plot current mode
if options.GnFigs
    try suffStat.haf = suffStat2.haf; end
    str = [s1,num2str(PreviousI,'%4.3e')];
    hf = figure('visible','off','color',[1,1,1]);
    pos = get(hf,'position');
    set(hf,'position',pos-[pos(3)./2 0 0 0],'visible','on')
    ha = axes('parent',hf);
    plot(ha,deltaMu')
    title(ha,[s2,' ; ',str])
    drawnow
end

% Regularized Gauss-Newton VB-Laplace update
it = 0;
stop = it>=options.GnMaxIter;
conv = 0;
posterior0 = posterior;
posterior = updatePosterior(posterior,PreviousMu,Sigma,indIn,flag);
while ~stop
    it = it+1;
    % make a move
    mu = PreviousMu + deltaMu;
    try
        % get next move and energy step
        [I,Sigma,NextdeltaMu,suffStat2] = feval(fname,mu,y,posterior,suffStat,dim,u,options);
        % get increment in variational/free energy
        [rdf,deltaI,F] = getCostIncrement(I,PreviousI,mu,Sigma,suffStat2,options,posterior,flag);
    catch
        rdf = -1;
        deltaI = -Inf;
    end
    % display move
    if options.GnFigs
        try, clf(hf); catch, hf = figure; end
        ha = axes('parent',hf);
        plot(ha,deltaMu')
        str = [s1,num2str(I,'%4.3e'),' ,dI/I =',num2str(rdf,'%4.3e'),' ,it #',num2str(it)];
    end
    VBA_pause(options)  % check 'pause' button
    % accept move or halve step?
    if deltaI<0     % halve step size
        deltaMu = 0.5*deltaMu;
        try,title(ha,[s2,': halve step ; ',str]);end
    else            % accept move
        % 1- propose a new move according to new local quadratic approx
        deltaMu = NextdeltaMu;
        % 2- update sufficient stats
        PreviousMu = mu;
        PreviousI = I;
        % 3- update posterior, model evidence and sufficient statistics
        posterior = updatePosterior(posterior,mu,Sigma,indIn,flag);
        if ~options.gradF
            [F] = VBA_FreeEnergy(posterior,suffStat2,options);
        end
        suffStat2.F = [suffStat2.F,F];
        suffStat = suffStat2;
        % 4- update display
        switch flag
            case {'X','X0'}
                VBA_updateDisplay(posterior,suffStat,options,y,[],'X')
            case 'Phi'
                VBA_updateDisplay(posterior,suffStat,options,y,[],'phi')
            case 'Theta'
                VBA_updateDisplay(posterior,suffStat,options,y,[],'theta')
        end
        try,title(ha,[s2,': accept move ; ',str]);end
        conv = 1;
    end
    % check convergence criterion
    if abs(rdf)<=options.GnTolFun || it==options.GnMaxIter
        stop = 1;
        try close(hf); end
        try close(suffStat.haf); end
    end
    drawnow
end
if ~conv
    suffStat.F = [suffStat.F,suffStat.F(end)];
    posterior = posterior0;
end


function posterior = updatePosterior(posterior,mu,Sigma,indIn,flag)
switch flag
    case 'X'
        posterior.muX = mu;
        posterior.SigmaX = Sigma;
    case 'X0'
        posterior.muX0(indIn) = mu;
        posterior.SigmaX0(indIn,indIn) = Sigma;
    case 'Phi'
        posterior.muPhi(indIn) = mu;
        posterior.SigmaPhi(indIn,indIn) = Sigma;
    case 'Theta'
        posterior.muTheta(indIn) = mu;
        posterior.SigmaTheta(indIn,indIn) = Sigma;
end



function [rdf,deltaI,F] = getCostIncrement(I,PreviousI,mu,Sigma,suffStat2,options,posterior,flag)
if options.gradF
    if ~isinf(I)
        switch flag
            case 'X'
                indIn = [];
                previousMu = posterior.muX;
            case 'X0'
                indIn = options.params2update.x0;
                previousMu = posterior.muX0(indIn);
            case 'Phi'
                indIn = options.params2update.phi;
                previousMu = posterior.muPhi(indIn);
            case 'Theta'
                indIn = options.params2update.theta;
                previousMu = posterior.muTheta(indIn);
        end
        posterior = updatePosterior(posterior,mu,Sigma,indIn,flag);
        [F] = VBA_FreeEnergy(posterior,suffStat2,options);
        deltaI = F - suffStat2.F(end);
        rdf = sum((previousMu(:)-mu(:)).^2./(mu(:).^2+eps));
    else
        F = -Inf;
        deltaI = -Inf;
        rdf = -Inf;
    end
else
    % calculate relative variational energy improvement
    deltaI = I - PreviousI;
    rdf = deltaI./abs(PreviousI);
    F = [];
end
