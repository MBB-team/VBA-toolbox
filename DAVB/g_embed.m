function [gx,G,dgdp] = g_embed(Xt,Phi,ut,in)
% observation function for dynamical system's dimension embedding
% function [gx,G,dgdp] = g_embed(Xt,Phi,ut,in)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

% First evaluate standard evolution function
Xt0 = Xt(1:in.dim.n);
[gx,G0,dgdp] = feval(in.g_fname,Xt0,Phi,ut,in.options.inG);
G = [ G0 ; zeros(in.dim.n_embed,size(G0,2)) ];
