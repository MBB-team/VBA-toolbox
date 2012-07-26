function [Lf] = models2families(Lm,partition)


nf = length(partition); % # families
[ns,nm] = size(Lm); % # models and # subjects

pm = zeros(nm,1);
for j=1:nf
    pm(partition{j}) = 1./length(partition{j});
end
pm = pm./sum(pm);

%hf = figure;
%subplot(2,2,1,'parent',hf),bar(pm),title('priors over models')

Lm = Lm - max(Lm(:));
pmy = exp(Lm)*diag(pm);
pmy = diag(sum(pmy,2).^-1)*pmy;
%subplot(2,2,2,'parent',hf),bar(pmy),title('posterior over models')

%figure,plot(sum(pmy,2),'.')

Lf = zeros(ns,nf);
pfy = zeros(ns,nf);
for j=1:nf
    pfy(:,j) = sum(pmy(:,partition{j}),2);
    Lf(:,j) = log(pfy(:,j));
end
%subplot(2,2,3,'parent',hf),bar(pfy),title('posterior over families')
%subplot(2,2,4,'parent',hf),bar(Lf),title('log-posterior over families')