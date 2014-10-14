function [gx,G,dgdp] = g_embed(Xt,Phi,ut,in)
% observation function for dynamical system's dimension embedding
% function [gx,G,dgdp] = f_embed(Xt,Phi,ut,in)

% First evaluate standard evolution function
Xt0 = Xt(1:in.dim.n);
[opt,dim] = getOptions4EvalFun(in);
[gx,G0,dgdp] = VBA_evalFun('g',Xt0,Phi,ut,opt,dim,0);
G = [ G0 ; zeros(in.dim.n_embed,size(G0,2)) ];


function [opt,dim] = getOptions4EvalFun(in)
opt = in.options;
opt.g_fname = in.g_fname;
opt.g_nout = in.g_nout;
opt.checkGrads = 0;
dim = in.dim;