function findlev= defIndlev(level, NtotPar)
% defIndlev returns a structure containing the indexes corresponding to the
% hidden states of a ToM agent of ToM=level (level) and a number of
% parameters (NtotPar) of the generative model of the opponent (i.e.
% volatility, temperature, bias...)
if level==0
    findlev=[];
else
curind=level-1;
for i=1:level %works only for i>0 
%careful here the order of hidden states is different from the equivalent
%evolution function for level>1 
% [estimated levels;simulated hidden state level 0; mu, Sigma para level 0; f0;df0;...;
% simulated hidden state level N-1; mu, Sigma para level N-1;
% fN-1;dfN-1]
    findlev(i).k=i-1;
    sizeXi=sizeXrec(i-1,NtotPar);
    findlev(i).X=curind+(1:sizeXi);%index simulated hidden states for level i-1
    curind=curind+sizeXi;
    findlev(i).f=curind+1;
    curind=curind+1;
    findlev(i).df=curind+(1:NtotPar);
    curind=curind+NtotPar;
    findlev(i).Par=curind+(1:2*NtotPar);%index parameters to update
    curind=curind+2*NtotPar;
end
end

end