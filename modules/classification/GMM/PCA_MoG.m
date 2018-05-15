function [u,ss,K,v]=PCA_MoG(y,varThresh,verbose)

if ~exist('varThresh','var')
    varThresh = 0.95;
end
if ~exist('verbose','var')
    verbose = 1;
end

n = size(y,2);
p = size(y,1);

[u,s,v] = svd(y',0);
ss = diag(s);

va = cumsum(ss.^2./sum(ss.^2));

ind = find(va>=varThresh);
K = min(ind);
u = u(:,1:K);

if ~verbose
    return
end

figure,subplot(2,1,1),plot(ss,'bo');
title('singular values');


subplot(2,1,2);plot(va,'bo');
hold on,plot([1 length(ss)],repmat(varThresh,1,2),'r');
title('percentage of explained variance')

disp([num2str(K),' significant components.']);

figure,imagesc(u(:,1:K)),colorbar,
title(['Data projection on the ',num2str(K),' significant principal axis']);


my = mean(y,2);
stdy = std(y(:));

f = figure;

for i = 1:K
    a = u(:,i);
    [n,x] = hist(a);
    figure(f),subplot(K,1,i), bar(x,n)
    title(['histogram of the projections on the ',num2str(i),'th PA'])
    
end

% sqrt(ss)

n = size(y,2);
xii = ones(n,1);

lab{1} = '1st PA';
for k = 2:K
    lab{k} = [num2str(k),'d PA'];
end
if p == 2
    pl = figure;hold on
    title('Data space')
    for i = 1:n
        hh(i) = plot(y(1,i),y(2,i),'k.');
    end
    xlabel('1st dim')
    ylabel('2d dim')
    for k = 1:K
        b = v(:,k)*sqrt(ss(k));
        h = quiver(my(1),my(2),b(1),b(2));
        for i=1:length(h)
            set(h(i),'color',[1 0 0],'LineWidth',2);
        end
        t=text(my(1)+b(1),my(2)+b(2),lab{k});
        set(t,'color',[1 0 0])
    end
elseif p >= 3
    pl = figure;hold on
    title('Data space (3 first dimensions)')
    for i = 1:n
        hh(i) = plot3(y(1,i),y(2,i),y(3,i),'k.');
    end
    xlabel('1st dim')
    ylabel('2d dim')
    zlabel('3d dim')
    for k = 1:K
        b = v(:,k)*sqrt(ss(k));
        h = quiver3(my(1),my(2),my(3),b(1),b(2),b(3));
        for i=1:length(h)
            set(h(i),'color',[1 0 0],'LineWidth',2);
        end
        t=text(my(1)+b(1),my(2)+b(2),my(3)+b(3),lab{k});
        set(t,'color',[1 0 0])
    end
end


if K == 1
elseif K == 2
    pl = figure;hold on
    title([num2str(K),'-D principal component space'])
    for i = 1:n
        hh(i) = plot(u(i,1),u(i,2),'k.');
    end
    xlabel('first PA')
    ylabel('second PA')
elseif K >= 3
    pl = figure;hold on
    title([num2str(K),'-D principal component space (3 first dimensions)'])
    for i = 1:n
        hh(i) = plot3(u(i,1),u(i,2),u(i,3),'k.');
    end
    xlabel('first PA')
    ylabel('second PA')
    zlabel('third PA')
end





