% demo Laplacian
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

hf = figure;
ha = axes('parent',hf,'nextplot','add');
N = 128;
dx = 1;
for n=2:N
    G = diag(ones(N,1),0)-diag(ones(N-1,1),1);
    L = -G'*G;
    L(N,N) = -1;
    L = L./(dx/n);
    ev = eig(L);
    plot(ha,ev)
    plot(ha,ev,'.')
    pause
end

% hf = figure;
% ha = axes('parent',hf,'nextplot','add');
% N = 64;
% for n=2:N
%     L = diag(ones(n,1),0)-diag(ones(n-1,1),1)./5;
%     ev = eig(L'*L);
%     plot(ha,ev)
%     plot(ha,ev,'.')
% end