function [gx] = g_AVL2(x,Phi,u,in)
% reaction times OTO observation function
% function [gx] = g_AVL2(x,Phi,u,in)
% This function computes the observation function (reaction time) for the
% audio-visual associative learning task. 
% IN:
%   - x: the 5x1 vector of states. The third first entries are sufficient
%   statistics of the updated belief. The two last ones are summary
%   statistics of the convergence rate of the VB updates.
%   - Phi: 2x1 vector of observation parameters. Phi(1) is an intercept
%   constant and Phi(2) is a scaling factor.
%   - u: useless.
%   - in: structure for flagging to different cost functions
% OUT:
%   - gx: expected reaction time
%   - dgdx: derivative wrt states

   
gx = x(end);