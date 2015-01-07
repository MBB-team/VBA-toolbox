function s=priorPrettifyer(labeler,muPriors)

names = fieldnames(labeler);
s = struct ;
% s=[];
for i=1:numel(names)
    s.(names{i}) = muPriors(labeler.(names{i}),:)';
%     s=[s sprintf('%s:\t\t%s \n',names{i},num2str(muPriors(labeler.(names{i}),:)'))];
end