function [eq,summary] = compare_struct(s1,s2)

fNames1 = fieldnames(s1);
fNames2 = fieldnames(s2);


summary.fieldSame = isempty(setdiff(fNames1,fNames2));

fieldSame = intersect(fNames1,fNames2);
for iField = 1:numel(fieldSame)
    fieldName = fieldSame{iField};
    summary.(fieldName) = isequaln(s1.(fieldName),s2.(fieldName));
end


eq = all(cell2mat(struct2cell(summary))) ;