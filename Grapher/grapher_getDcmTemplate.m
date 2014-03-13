function nodes = grapher_getDcmTemplate(pos,labels)


n_nodes = numel(labels) ;

for i=1:n_nodes
    nodes(i).n = n_nodes;
    nodes(i).id = i;
    if ~isempty(pos)
        nodes(i).pos = pos(i,:);
    else
        thet = (2*pi) * i/n_nodes ;
        nodes(i).pos = 300*[cos(thet) sin(thet)] ;
    end
    if ~isempty(labels)
        nodes(i).lbl = labels{i};
    else
        nodes(i).lbl = num2str(i) ;
    end
end

end