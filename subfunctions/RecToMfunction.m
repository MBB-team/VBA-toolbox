function [fx,indlev] =RecToMfunction(x,P,u,inF)
% Marie Devaine wrote this function in November 2015
% Recursive function to do the evolution update of any ToM level (here in
% the context of a game) but the form is quite general.
%The update of the hidden states is done by calling recursively the
%function
%The update of level estimation and Parameters estimations are performed at
%each step by the function
%update of the Parameters is performed by the updatePar function
d=0.1;
b=-0.4;
c=0.7;
level=inF.lev;
fx = zeros(size(x));
if level==0
    fx=evolution0bisND(x,P,u,inF);
    indlev=[];
else
    fobs=inF.fobs;
    NParev=inF.indParev; %=1
    NParobs=inF.indParobs;
    NtotPar=inF.indParev+inF.indParobs;
    indlev = defIndlev(level, NtotPar);
    if ~isempty(u)
        y2t=u(1);
        if level>1
            vecV=zeros(1,level);% V factor used in the update of the level of adversary
            f=zeros(1,level);
            for j=1:level
                f(j)=x(indlev(j).f);
                df=x(indlev(j).df);
                for i=1:NtotPar
                    Sig=exp(x(indlev(j).Par(2*i)));
                    vecV(j)=vecV(j)+Sig*df(i)^2;
                end
            end
            
            %update estimate level adversary
            vecp=zeros(1,level);
            vecp(1:(level-1))=sigmoid(x(1:(level-1)));
            vecp(end)=max(.0001,1-sum(vecp));
            Nvecp= vecp.*sigmoid((f+b*vecV.^c)./sqrt(1+d*vecV)).^y2t.*(sigmoid((-f+b*vecV.^c)./sqrt(1+d*vecV))).^(1-y2t);
            Nvecp= Nvecp./sum(Nvecp);
            fx(1:(level-1))=invsigmoid(Nvecp(1:(end-1))); %just keep the first
            
        else
            Nvecp=1;
        end
    else %initialisation
        Nvecp= 1./level.*ones(1,level);
        if level>1
            fx(1:(level-1)) =invsigmoid(Nvecp(1:(end-1)));
        end
    end
    %update estimate Parameter and hidden states
    
    for j=1:level
        in.k=Nvecp(j); %estimate level opponent
        in.f=x(indlev(j).f);
        in.df=x(indlev(j).df);
        in.dummyPar=inF.dummyPar;%1 by default otherwise do not use volatility to update
        
        if ~isempty(u) %else just initialize accordingly f and df
        %%%%%%%%%%% Update parameters
        fx(indlev(j).Par) = updatePar(x(indlev(j).Par),P,y2t,in);
        u_t=u([2;1]);
        else
           fx(indlev(j).Par)=x(indlev(j).Par);  
           u_t=[];
        end
        inL=inF;
        inL.lev=j-1;
        inL.player=3-inF.player;
        
        %%%%%%%%%%% Update parameters hidden states
        Parjev = fx(indlev(j).Par(mod(1:2*NParev,2)==1)); %take mu and not sigma
        Parobs = fx(indlev(j).Par(find(mod((1:2*NParobs),2)==1)+2*NParev)); %take mu and not sigma
        
        [X, indlevj]  = RecToMfunction(x(indlev(j).X),Parjev,u_t,inL);
        fx(indlev(j).X)=X;
        %computing derivatives
        inG=inF;
        inG.npara=NtotPar;
        inG.lev=j-1;
        inG.k=in.k;
        inG.player= 3-inF.player;
        inG.indlev=indlevj;
        
        f_new=invsigmoid(fobs(X,Parobs,u_t,inG));
        fx(indlev(j).f)=f_new;
        eps=zeros(1,NParev);
        f_neweps=zeros(NParev,1);
        df_newev=zeros(NParev,1);
        %%%%% the for loop is useless in the usual case where NParev==1;
        for l=1:NParev
            eps(l)=1e-4*Parjev(l);
            if abs(eps(l))<1e-4
                eps(l)=1e-4;
            end
            Parjevl=Parjev;
            Parjevl(l)=Parjevl(l)+eps(l);
            Xeps = RecToMfunction(x(indlev(j).X),Parjevl,u_t,inL);
            f_neweps(l)=invsigmoid(fobs(Xeps,Parobs,u_t,inG));% index !
            df_newev(l)=(f_neweps(l)-f_new)/eps(l);
        end
        fx(indlev(j).df(1:NParev))=df_newev;
        eps=zeros(1,NParobs);
        f_neweps=zeros(NParobs,1);
        df_newobs=zeros(NParobs,1);
        %    df_newobs= -f_new;
        %%%%% in general use below
        for l=1:NParobs
            eps(l)=1e-4*Parobs(l);
            if abs(eps(l))<1e-4
                eps(l)=1e-4;
            end
            Parobsl=Parobs;
            Parobsl(l)=Parobsl(l)+eps(l);
            f_neweps(l)=invsigmoid(fobs(X,Parobsl,u_t,inG));
            df_newobs(l)=(f_neweps(l)-f_new)/eps(l);
        end
        fx(indlev(j).df(NParev+(1:NParobs)))=df_newobs;
    end
end
end


function fPar_k=updatePar(Par_k,theta,u_t,in) %
fPar_k=zeros(size(Par_k));
p=in.k; % weight for levels
f=in.f;
Proba=sigmoid(f);
df=in.df;
V=(exp(Par_k(mod(1:end,2)==0))+exp(theta).*in.dummyPar).^-1;
Sig=1./(V+p*Proba*(1-Proba).*(df.^2));
fPar_k(mod(1:end,2)==0)=log(Sig);
fPar_k(mod(1:end,2)==1)=Par_k(mod(1:end,2)==1)+p.*Sig.*(u_t-Proba).*df;
end
