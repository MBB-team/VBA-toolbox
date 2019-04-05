function [fx,J,dfdp] = f_embed0(Xt,Theta,ut,in)
% evolution function for dynamical system's delay embedding
% function [fx,J,dfdp] = f_embed(Xt,Theta,ut,in)

% First form delayed state vector and input
iX = (1:in.dim.n) + (in.options.delays(:)'.*in.dim.n);
dX = Xt(iX,:);

% Evaluate evolution function at the delayed state vector:
[opt,dim] = getOptions4EvalFun(in);
[fx0,J0,dfdp0] = VBA_evalFun('f',dX,Theta,ut,opt,dim);

% Construct full embedding flow and gradients:
Xe = Xt(1:in.dim.n_embed);
fx = [fx0;Xe];
J = zeros(in.dim.n_embed+in.dim.n);
J(iX,1:in.dim.n) = J0;
J(1:in.dim.n_embed,in.dim.n+1:in.dim.n+in.dim.n_embed) = ...
    eye(in.dim.n_embed);
dfdp = [dfdp0,zeros(in.dim.n_theta,in.dim.n_embed)];
       
        
function [opt,dim] = getOptions4EvalFun(in)
opt = in.options;
opt.f_fname = in.f_fname;
opt.f_nout = in.f_nout;
opt.checkGrads = 0;
dim.n = in.dim.n;
dim.n_theta = in.dim.n_theta;