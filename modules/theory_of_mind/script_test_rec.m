%script test recursive
%%%%% quantitative differences: because of the way the level is updated
%%%%% (should be analytically equivalent but not numerically). Small
%%%%% deiations for 2-ToM (but accumulates) larger for 3-ToM
x=zeros(2,101); xtest=x;
y1=VBA_random ('Bernoulli', 0.7, 100, 1);
inF.player=1;
inF.lev=0;
Par=0; %theta
for i=1:100
    x(:,i+1)=RecToMfunction(x(:,i),Par,y1(i),inF);
    xtest(:,i+1)=evolution0bisND(xtest(:,i),Par,y1(i),inF);
end
%%
N=1000
xtest=zeros(8,N+1); x=zeros(9,N+1);
y1=VBA_random('Bernoulli', 0.7, N, 1);
y2=VBA_random('Bernoulli', 0.5, N, 1);
HS=cat(3,[1,0;0,1],[0,1;1,0]);
inF.player=1;
inF.lev=1;
inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1;
inF.indParobs=1;
inF.dummyPar=[1;0];
in={HS,1};
Par=0; %theta
for i=1:N
    [x(:,i+1),indlev]=RecToMfunction(x(:,i),Par,[y2(i);y1(i)],inF);
    xtest(:,i+1)=evolution1bisND(xtest(:,i),Par,[y2(i);y1(i)],in);
end

%%%%compare
indcorr1=[1:2,6,8,7,9,3:4];
xcomp=x(indcorr1,:);
find(abs(xcomp-xtest)>1e-15)
%%
N=50;

y1=VBA_random('Bernoulli',0.7,100,1);
y2=VBA_random('Bernoulli',0.5,100,1);
%%
xtest=zeros(23,N+1); x=zeros(26,N+1);
HS=cat(3,[1,0;0,1],[0,1;1,0]);
inF.player=1;
inF.lev=2;
inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1;
inF.indParobs=1;
inF.dummyPar=[1;0];
in={HS,1};
Par=0; %theta
for i=1:N
    [x(:,i+1),indlev]=RecToMfunction(x(:,i),Par,[y2(i);y1(i)],inF);
    xtest(:,i+1)=evolution2bisND(xtest(:,i),Par,[y2(i);y1(i)],in);
end
indcorr2=[2:5,10+indcorr1,20,21,7,9,8,10,23,25,24,26,1];
xcomp=x(indcorr2,:);
xcomp(end,:)=-xcomp(end,:);
%%%%compare
find(abs(xcomp-xtest)>1e-5)
%%
N=50;
xtest=zeros(53,N+1); x=zeros(60,N+1);
y1=VBA_random('Bernoulli',0.7,100,1);
y2=VBA_random('Bernoulli',0.5,100,1);
HS=cat(3,[1,0;0,1],[0,1;1,0]);
inF.player=1;
inF.lev=3;
inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1;
inF.indParobs=1;
inF.dummyPar=[1;0];
in={HS,1};
Par=0; %theta
for i=1:3
    [x(:,i+1),indlev]=RecToMfunction(x(:,i),Par,[y2(i);y1(i)],inF);
    xtest(:,i+1)=evolution3bisND(xtest(:,i),Par,[y2(i);y1(i)],in);
end
 indcorr3=[3:6,11+indcorr1,21,22,27+indcorr2,54,55,8,10,9,11,24,26,25,27,57,59,58,60,1,2];
 xcomp=x(indcorr3,:);
 xcomp(37,:)=-xcomp(37,:); %estimation level 2 

%%%%compare
find(abs(xcomp-xtest)>1e-4)

%%
N=10;
inF.lev=7;
%xtest=zeros(53,N+1); 
x=zeros(sizeXrec(inF.lev,2),N+1);
y1=VBA_random('Bernoulli',0.7,100,1);
y2=VBA_random('Bernoulli',0.5,100,1);

HS=cat(3,[1,0;0,1],[0,1;1,0]);
inF.player=1;

inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1;
inF.indParobs=1;
inF.dummyPar=[1;0];
in={HS,1};
Par=0; %theta
for i=1:N
    [x(:,i+1),indlev]=RecToMfunction(x(:,i),Par,[y2(i);y1(i)],inF);
   % xtest(:,i+1)=evolution3bisND(xtest(:,i),Par,[y2(i);y1(i)],in);
end
%  indcorr3=[3:6,11+indcorr1,21,22,27+indcorr2,54,55,8,10,9,11,24,26,25,27,57,59,58,60,1,2];
%  xcomp=x(indcorr3,:);
%  xcomp(37,:)=-xcomp(37,:); %estimation level 2 
% 
% %%%%compare
% find(abs(xcomp-xtest)>1e-10)
%%
N=100;
 x=zeros(264,N+1);
y1=VBA_random('Bernoulli',0.7,100,1);
y2=VBA_random('Bernoulli',0.5,100,1);
HS=cat(3,[1,0;0,1],[0,1;1,0]);
inF.player=1;
inF.lev=5;
inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1;
inF.indParobs=1;
inF.dummyPar=[1;0];
in={HS,1};
Par=0; %theta
for i=1:N
    [x(:,i+1),indlev]=RecToMfunction(x(:,i),Par,[y2(i);y1(i)],inF);
end
