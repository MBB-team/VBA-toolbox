function accuracy=balanced_accuracy(pred,data,isOut)

if ~exist('isOut','var')
    isOut = zeros(size(pred));
end

pred=vec(pred(isOut==0));
data=vec(data(isOut==0));


bg = pred > 0.5;                   % binarized model predictions
tp = sum(   data   .*   bg    );   % true positives
fp = sum( (1-data) .*   bg    );   % false positives
fn = sum(   data   .* (1-bg)  );   % false positives
tn = sum( (1-data) .* (1-bg)  );   %true negatives
P = tp + fn;
N = tn + fp;

accuracy = 0.5*(tp./P + tn./N);

end