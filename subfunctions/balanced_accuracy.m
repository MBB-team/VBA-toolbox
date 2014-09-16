function acc=balanced_accuracy(pred,data)
pred=vec(pred);
data=vec(data);
    pos_idx = find(pred>=.5);
    neg_idx = find(pred<=.5);    
    acc = mean( [mean(data(pos_idx)) , mean(1-data(neg_idx))] );
end