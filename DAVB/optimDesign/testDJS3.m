function testDJS3
% test the influence of data dimension on Chernoff bound (with variance scaling)


alpha1 = 1;
alpha2 = 4;
ps = 0.5*ones(2,1);


gridn = 1:4;

DJS = zeros(length(gridn),1);
b = DJS;
pe = DJS;

for i = 1:length(gridn)
    
    In = eye(i);
    mus{1} = zeros(i,1)-1;
    Qs{1} = alpha1.*In;
    mus{2} = mus{1};
    Qs{2} = alpha2.*In;
    
    [DJS(i),b(i)] = JensenShannon(mus,Qs,ps);
    [pe(i)] = ProbErrorND(mus,Qs,ps);

end


hf = figure('color',[1 1 1]);
pos = get(hf,'position');
set(hf,'position',[pos(1)-pos(3)/2 pos(2) 2*pos(3) pos(4)])
ha2 = subplot(1,2,2,'parent',hf);
set(ha2,...
    'nextplot','add',...
    'xlim',[0.5,length(gridn)+.5],...
    'xtick',1:length(gridn),...
    'xticklabel',gridn,...
    'fontsize',14)
plot(ha2,pe,'g')
plot(ha2,-DJS,'b')
plot(ha2,0.5*b,'r')
plot(ha2,0.25*b.^2,'r--')
plot(ha2,-DJS,'b.')
plot(ha2,0.5*b,'r.')
plot(ha2,0.25*b.^2,'r.')
plot(ha2,pe,'g.')
grid(ha2,'on')
xlabel(ha2,'data dimension')
ylabel(ha2,'p(e=1|u)')
legend(ha2,{'p(e=1|u)','r(u)','b(u)/2','b(u)^2/4'})
getSubplots

