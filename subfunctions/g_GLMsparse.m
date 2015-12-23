function [gx] = g_GLMsparse(x,P,u,in)
sP = sparseTransform(P,in.sparseP);
[gx] = g_GLM(x,sP,u,in);
