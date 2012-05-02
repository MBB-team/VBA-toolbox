function [g,dgdx,dgdP] = g_GLM(x,P,u,in)

if isempty(u) || isequal(u,0)
    g = in.X*P;
    dgdx = [];
    dgdP = in.X';
else
    g = in.X(u,:)*P;
    dgdx = [];
    dgdP = in.X(u,:)';
end
