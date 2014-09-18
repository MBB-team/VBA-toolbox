function [inF, n_r] = extend_dcm(inF,hA,hB,hC,hD,dim,sources)

%%- get dimensions
n=dim.n;
n_u=dim.n_u;
n_r = get_dims(hA,hB,hC,hD) ;

sourcesDim = cellfun(@numel,{sources(2:end).out});
for i=1:sum(sourcesDim)
    respSource(i) = sum(i>cumsum(sourcesDim))+1 ;
end


%%
[inFtemp] = prepare_dcm(hA,hB,hC,hD,respSource);

%% = save extended structure
%- raw matrices
inF.hA = inFtemp.A ;
inF.hB = inFtemp.B ;
inF.hC = inFtemp.C ;
inF.hD = inFtemp.D ;

%- parameter indices
inF.indhA = inFtemp.indA + inF.indself ;

for i=1:n_u
    inF.indhB{i} = inFtemp.indB{i} + inF.indself ;
end
inF.indhC = inFtemp.indC + inF.indself ;
for i=1:n
    inF.indhD{i} = inFtemp.indD{i} + inF.indself ;
end
inFtemp.indself = inFtemp.indself ;
inF.indhself = inFtemp.indself + inF.indself ;

%- indicators
inF.dhA = inFtemp.dA;
inF.dhB = inFtemp.dB;
inF.dhC = inFtemp.dC;
inF.dhD = inFtemp.dD;

%- save dimensions
% dim.n_r=n_r;

end


%% ========== subfunctions =============

%- compute number of predicted responses
function nr=get_dims(hA,hB,hC,hD)
nr=0;
nr=max(nr,size(hA,1));
for i=1:length(hB)
    nr=max(nr,size(hB{i},1));
end
nr=max(nr,size(hC,1));
for i=1:length(hD)
    nr=max(nr,size(hD{i},1));
end
end



