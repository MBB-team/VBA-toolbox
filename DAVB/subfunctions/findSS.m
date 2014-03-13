function [xf] = findSS(f_fname,x0,P,u,in)

in.f_fname = f_fname;

%   - options: a structure that may contains the following fields:
%       .args: a cell array containing the arguments of costFcn
%       .GnMaxIter: the max number of iterations for the main regularized
%       Gauss-Newton optimization scheme
%       .GnTolFun: the min relative increase in costFcn
%       .minimize: {0} if 1, optimCost MINIMIZES the function

options = struct(...
    'GnMaxIter',4e2,...
    'GnTolFun',1e-4,...
    'minimize',1,...
    'verbose',0);
options.args = {P;u;in};

[xf,curv,out] = optimCost(@squaredEv,x0,options);



function f2 = squaredEv(x,P,u,in)
fx = feval(in.f_fname,x,P,u,in);
f2 = fx-x;
f2 = f2'*f2;


