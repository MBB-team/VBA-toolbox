N=100;
%y2 = bernoulli(.25,N);
y1 = NaN(1,N);
level=3;
HS=cat(3,[1,0;0,1],[0,1;1,0]);
inF.player=2;
inF.lev=level;
inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1; % Nb para evol
inF.indParobs=2; %Nb para obs
inF.dummyPar=[1;0;0];
NtotPar=inF.indParev +inF.indParobs;

phi=[0;1];
theta=-2;
%initialize hidden states
matx1= zeros(sizeXrec(level,NtotPar),N+1);
[matx1(:,1), indlev1] =RecToMfunction(matx1(:,1),theta, [],inF);
gx=NaN(1,N);
%indlev1=defIndlev(level,NtotPar);
inG=inF;
inG.indlev=indlev1;
inG.npara=NtotPar;
%%

for i=1:N
    gx(i)= ObsRecGen(matx1(:,i)',phi, [],inG) ; 
    y1(i)= gx(i)>.5;
    matx1(:,i+1)= RecToMfunction(matx1(:,i),theta, [y2(i);y1(i)],inF);
    
end