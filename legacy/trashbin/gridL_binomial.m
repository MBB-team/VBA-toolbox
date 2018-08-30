function [logP,bestPhi] = gridL_binomial(y,u,f_fname,g_fname,dim,options)

% nuemrically evaluates the posterior of (binomial) likelihood parameters

[options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options);
iQ = VBA_inv(options.priors.SigmaPhi);

gridPhi1 = -3:0.2e-1:3;
gridPhi2 = -0.5:1e-2:0.5;
n1 = length(gridPhi1);
n2 = length(gridPhi2);
logP = zeros(n1,n2);

if options.verbose
    fprintf(1,'Numerical evaluation of the log-posterior ...')
    fprintf(1,'%6.2f %%',0)
end

for i=1:n1
    if options.verbose
        fprintf(1,'\b\b\b\b\b\b\b\b')
        fprintf(1,'%6.2f %%',100*i/n1)
    end
    for j=1:n2
        Phi = [gridPhi1(i);gridPhi2(j)];
        logL = 0;
        for t=1:dim.n_t
            gx(:,t) = VBA_evalFun('g',[],Phi,u(:,t),options,dim);
            gx(:,t) = checkGX(gx(:,t));
            logL =  logL ...
                + y(:,t)'*log(gx(:,t)) + (1-y(:,t))'*log(1-gx(:,t));
        end
        dPhi = Phi - options.priors.muPhi;
        logP(i,j) = logL -0.5*dPhi'*iQ*dPhi;
    end
end

[~, iy, ix] = VBA_maxMat(logP);

bestPhi = [gridPhi1(iy),gridPhi2(ix)];

if options.verbose
    fprintf(1,'\b\b\b\b\b\b\b\b')
    fprintf(' OK.')
    fprintf('\n')
end

function x = checkGX(x)
lim = 1e-8;
x(x<=lim) = lim;
x(x>=1-lim) = 1-lim;
