function [gx] = g_LinDecomp(x,P,u,in)
X = NaN(in.size.X(1),in.size.X(2));
Y = NaN(in.size.Y(1),in.size.Y(2));
for i=1:in.n
    X(:,i) = P(in.ind(i).X);
    Y(i,:) = P(in.ind(i).Y)';
end
gx = VBA_vec(X*Y) + P(in.ind0);