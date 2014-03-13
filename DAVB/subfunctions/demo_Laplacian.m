function demo_Laplacian

hf = figure('name','Laplacian eigenvalues');
ha = axes('parent',hf,'nextplot','add');
N = 16;
dx = 1;
col = getColors(N);
str = cell(0);
for n=1:N
    G = diag(ones(N,1),0)-diag(ones(N-1,1),1);
    L = -G'*G;
    L(N,N) = -1;
    L = L./(dx/n)^2;
    ev = flipud(eig(L));
    hp = plot(ha,ev,'color',col(n,:),'marker','.');
%     plot(ha,ev,'.','color',col(n,:))
    str{end+1} = ['dx=',num2str(dx/n)];
end

legend(str)

% hf = figure;
% ha = axes('parent',hf,'nextplot','add');
% N = 64;
% for n=2:N
%     L = diag(ones(n,1),0)-diag(ones(n-1,1),1)./5;
%     ev = eig(L'*L);
%     plot(ha,ev)
%     plot(ha,ev,'.')
% end