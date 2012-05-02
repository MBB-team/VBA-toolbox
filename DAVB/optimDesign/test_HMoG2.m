function test_HMoG2
% test the influence of number of models on Chernoff bound (with shifted mean)

alpha = 1;
dx = 0:1:8;

gri = -10:0.01:10;
col = repmat('rgbcmy',1,1e2);

n = 1;
In = eye(n);

mus{1} = zeros(n,1)-4;
Qs{1} = alpha.*In;
tmp = mus{1}-gri;
p1 = exp(-0.5*tmp.^2./Qs{1});
p1 = p1./sum(p1);

hf = figure('color',[1 1 1]);
pos = get(hf,'position');
set(hf,'position',[pos(1)-pos(3)/2 pos(2) 2*pos(3) pos(4)])
ha = subplot(1,2,1,'parent',hf);
set(ha,...
    'nextplot','add',...
    'xlim',[min(gri) max(gri)],...
    'fontsize',14);

str = cell(0);
xt = cell(0);
DJS = zeros(length(dx),1);
b = DJS;
pe = DJS;

for i = 1:length(dx)

    
    mus{i+1} = mus{1}+dx(i);
    Qs{i+1} = Qs{1};
    
    tmp = mus{i+1}-gri;
    p2 = exp(-0.5*tmp.^2./Qs{i+1});
    p2 = p2./sum(p2);
    
    plot(ha,gri,p2,'color',col(i+1));
    str{i} = ['m_2: mu = ',num2str( mus{i+1})];
    legend(str)
    
    ps = ones(i+1,1)./(i+1);
    [DJS(i),b(i)] = JensenShannon(mus,Qs,ps);
    [pe(i)] = ProbError(mus,Qs,ps);
    
    
end
plot(ha,gri,p1,'k--');
str{end+1} = ['m_1: mu = ',num2str(mus{1})];
legend(ha,str)
grid(ha,'on')
xlabel(ha,'y')
ylabel(ha,'p(y|m)')

ha2 = subplot(1,2,2,'parent',hf);
set(ha2,...
    'nextplot','add',...
    'xlim',[0.5,length(dx)+.5],...
    'xtick',1:length(dx),...
    'xticklabel',2:length(dx)+1,...
    'fontsize',14)
plot(ha2,pe,'g')
plot(ha2,-DJS,'b')
plot(ha2,0.5*b,'r')
plot(ha2,b.^2./(4*[1:length(dx)]'),'r--')
plot(ha2,-DJS,'b.')
plot(ha2,0.5*b,'r.')
plot(ha2,b.^2./(4*[1:length(dx)]'),'r.')
plot(ha2,pe,'g.')
grid(ha2,'on')
xlabel(ha2,'number of models')
ylabel(ha2,'p(e=1|u)')
legend(ha2,{'p(e=1|u)','r(u)','b(u)/2','b(u)^2/4(M-1)'})
getSubplots

