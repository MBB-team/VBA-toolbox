function [mu,curv,out] = VBA_optimCost(costFcn,init,options)
% generic call function for function f(x,args) maximization
% function [opt,curv,out] = VBA_optimCost(costFcn,init,options)
% IN:
%   - costFtn: the name/handle of the function that computes the cost,
%   taking x and some arguments 'args'. Its I/O format should be:
%   I = costFcn(x,args{1},args{2},...,args{n})
%   - init: the initialization of the main argument x
%   - options: a structure that may contains the following fields:
%       .args: a cell array containing the arguments of costFcn
%       .GnMaxIter: the max number of iterations for the main regularized
%       Gauss-MNewton optimization scheme
%       .GnTolFun: the min relative increase in costFcn
%       .minimize: {0} if 1, optimCost MINIMIZES the function
% OUT:
%   - mu: the optimal x
%   - curv: the curvature of costFcn evaluated at mu
%   - out: a structure containing:
%       .it: the number of iterations of the regularized Gauss-Newton
%       optimization scheme
%       .conv: =1 if optimization has converged, =0 if not.
%       .nReg: the number of times Levenberg-Marquardt regularization has
%       been applied.
%       .I: the series of cost function evaluations along the iterative
%       optimization process
%       .elapsedTime = what do you think (in sec)?


try
    minimize = options.minimize;
catch
    minimize = 0;
end
opt.args = {costFcn;minimize};
try
    opt.args = cat(1,opt.args,options.args(:));
end
try
    opt.GnMaxIter = options.GnMaxIter;
catch
    opt.GnMaxIter = 256;
end
try
    opt.GnTolFun = options.GnTolFun;
catch
    opt.GnTolFun = 1e-4;
end
try
    opt.verbose = options.verbose;
catch
    opt.verbose = 0;
end
tic
[mu,sigma,out] = VBA_GaussNewton(@objFun,init,opt);
curv = -pinv(sigma);
out.elapsedTime = toc;
    


function [I,S,dx] = objFun(x,myFun,minimize,varargin)
args = varargin;
I = myfun(x,myFun,minimize,args);
dIdx = VBA_numericDiff(@myfun,1,x,myFun,minimize,args);
d2Idx2 = VBA_numericDiff(@VBA_numericDiff,3,@myfun,1,x,myFun,minimize,args);
S = -VBA_inv(full(d2Idx2));
dx = S*dIdx;
function I = myfun(x,myFun,minimize,args)
Args = cat(1,x,args(:));
I = (-1)^minimize.*feval(myFun,Args{:});


