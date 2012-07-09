load sujet2p
INVERSION2p = INVERSION;
load sujet2Q
INVERSION2Q = INVERSION;


Nsubjects = length(INVERSION);

LogEv = zeros(2,Nsubjects);
for i = 1 : Nsubjects
    
    LogEv(1,i)=INVERSION2Q{i}.results.out.F;
    LogEv(2,i)=INVERSION2p{i}.results.out.F;
    
end

figure
bar(LogEv(1,:)-LogEv(2,:))
title('Difference of log evidence (M1-M2)')
xlabel('Subject inder')


% From the inversion, model 2 seems to be favored.