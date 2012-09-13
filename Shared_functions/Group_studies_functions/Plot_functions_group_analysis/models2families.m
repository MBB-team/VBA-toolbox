function [Lf] = models2families(Lm,partition)


nf = length(partition); % # families

%--- restricting to models declared in families
nm = size(Lm,2); % initial size of log evidence matrix
i_map = zeros(1,nm);
i_models = sort(cell2mat(partition),'ascend'); % indices of the models considered
nm = length(i_models); % numbers of models considered

for i = 1 : length(i_models)
    i_map(i_models(i)) = i; %mapping initial model indices to their new ones
end
for j=1:nf
    partition{j} = i_map(partition{j});% changing partition indices accordingly
end
Lm = Lm(:,i_models); % reducing matrix of logev


%--- Computing flat priors on each partition
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