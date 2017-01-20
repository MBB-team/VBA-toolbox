function [out] = mediation_contrast(Y,X,M,c,options)
% performs contrast-based mediation analyses
% function [out] = mediation_contrast(Y,X,M)
% Let us consider the following regression equations:
% M = X*a + e_M
% Y = M*b + X*c + e_Y
% where X is the design matrix (containing the independent variables), M is
% the mediator variable and Y is the dependent variable.
% Let c be a contrast of interest, i.e. we want to test whether the effect
% of c*b onto Y is mediated by M. This is the statistical test performed by
% this function.
% IN:
%   - Y: nx1 dependent variable
%   - X: nxk design matrix (including independent variables)
%   - M: the mediator variable
%   - c: kxm contrast matrix (default is eye(k) -> omnibus F-test)
%   - options: an optional structure containing the following fields:
%       .alpha: significance level {0.05}
%       .verbose: verbose mode {1}
% OUT:
%   - out: output structure containing the following fields:
%       .p: mediation p-value (conjunctive testing approach)
%       .p0: p-value of the marginal contrast effect X->Y (F-test)
%       .pa: p-value of the contrast effect X->M (F-test)
%       .pb: p-value of the marginal effect M->Y (F-test)


n = size(Y,1);
if size(X,1)~= n || size(M,1)~= n
    disp('Error: input variables must have the same size!')
    out = [];
    return
end

try;c;catch;c=eye(size(X,2));end
try;alpha=options.alpha;catch;alpha=0.05;end
try;verbose=options.verbose;catch;verbose=1;end

% Y = Xa + e --> test for H0: c(:,i)*a=0
[out.p0,s0,df0,all0] = GLM_contrast(X,Y,c,'F',0);

% M = Xa + e --> test for H0: c*a=0
[out.pa,sa,dfa,alla] = GLM_contrast(X,M,c,'F',0);

% Y = Xa + Mb + e --> test for H0: b=0
[out.pb,sb,dfb,allb] = GLM_contrast([X,M],Y,[zeros(size(X,2),1);1],'F',0);

% loop over contrast components
pa = NaN(size(c,2),1);
pc = NaN(size(c,2),1);
for i=1:size(c,2)
    % M = Xa + e --> test for H0: c(:,i)*a=0
    [pa(i)] = GLM_contrast(X,M,c(:,i),'F',0);
    % conjunctive test
    pc(i) = max([pa(i),out.pb]);
end
% global null test
out.p = max([min(pc).*size(c,2),1]);

% summary
if verbose
    if out.p < alpha
        if out.p0 < alpha
            os = 'partial mediation';
        else
            os = 'full mediation';
        end
    else
        os = 'no mediation';
    end
    disp(' ')
    disp(['Date: ',datestr(clock)])
    disp(' ')
    disp(['-- Results (',os,', at alpha=',num2str(alpha),') --'])
    disp(['0) Regression of Y on X: R2=',num2str(round(1e3*all0.R2_a)/10),'% (p=',num2str(out.p0),')'])
    disp(['1) Regression of M on X: R2=',num2str(round(1e3*alla.R2_a)/10),'% (p=',num2str(out.pa),')'])
    disp(['2) Regression of Y on M (with X): R2=',num2str(round(1e3*allb.R2_a)/10),'% (p=',num2str(out.pb),')'])
    disp(['3) Mediation (conjunctive testing): p=',num2str(out.p),')'])
end


