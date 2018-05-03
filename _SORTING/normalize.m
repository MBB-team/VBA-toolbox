function yn = normalize(y)

yn = zeros(size(y));
for i=1:size(y,2)
    y0 = y(:,i);
    my = mean(y0);
    sy = std(y0);
    yn(:,i) = (y0 - my)./sy;
end
