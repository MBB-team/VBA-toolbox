function y = approxOnGrid(x,gy)
% approximates entry x on grid gy

sx = size(x);
x = x(:);
gy = gy(:);

d = zeros(length(x),length(gy));
if length(x) >= length(gy)
    for i=1:length(gy)
        d(:,i) = abs(x-gy(i));
    end
else
    for i=1:length(x)
        d(i,:) = abs(x(i)-gy)';
    end
end
[C,I] = min(d,[],2);
y = gy(I);
y = reshape(y,sx);

