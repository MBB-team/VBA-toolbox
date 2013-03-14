function [options] = prepare_extendedDCM(options,hA,hB,hC,hD)
% precalculates intermediary variables for the extension of DCM for
% behavioraldata, update the options accordingly
% function [options] = prepare_fullDCM(options,hA,hB,hC,hD)
% IN:
%   - options: 
%       structure containing DCM specification (see prepare_fullDCM)
%   - hX :
%       binary matrices 
% OUT:
%   - options: 
%       incomplete optimal structure for VB inversion of the
%       specified model (this does not include priors)...

%%
options.extended = true;

%% recompute inF for extended matrices
[n,nu,nr] = get_dims(options.inF,hA,hB,hC,hD) ;
    A = options.inF.A;
    A(n+1:nr,n+1:nr) = 0;
    B = options.inF.B;
    for i=1:length(B)
        B{i}(n+1:nr,n+1:nr) = 0;
    end
    C = options.inF.C;
    C(n+1:nr,size(C,2):nu) = 0;
    D = options.inF.D;
    for i=1:length(D)
        D{i}(n+1:nr,n+1:nr) = 0;
    end
    if options.inF.fullDCM
        options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
    else
        inF = prepare_DCM(A,B,C,D);
        options.inF=inF;
    end
    

% :: indices in the trimmed and vectorized description
[indhA,indhB,indhC,indhD,indhself,n_psi,n_r] = find_hDCM(hA,hB,hC,hD);



% :: corresponding indicator matrices in vectorized system
[inH.dhU,inH.dhU2,inH.dhX,inH.dhX2,inH.dhBi] = get_dMatdvec(inH);

inF=options.inF;
inF.hA = hA;
inF.hB = hB;
inF.hC = hC;
inF.hD = hD;


end

%% subfunctions

function [n,nu,nr]=get_dims(inF,hA,hB,hC,hD) %compute number of predicted responses
    % DCM
    n=0;
    nu=0;
    n=max(n,size(inF.A,1));
    nu=max(nu,length(inF.B));
    for i=1:length(inF.B)
        n=max(n,size(inF.B{i},1));
    end
    n=max(n,size(inF.C,1));
    nu=max(nu,size(inF.C,2));
    for i=1:length(inF.D)
        n=max(n,size(inF.D{i},1));
    end
    
    % decoding
    nr=0;
    nr=max(nr,size(hA,1));
    nu=max(nu,length(hB));
    for i=1:length(hB)
        nr=max(nr,size(hB{i},1));
    end
    nr=max(nr,size(hC,1));
    nu=max(nu,size(hC,2));
    for i=1:length(hD)
        nr=max(nr,size(hD{i},1));
    end
end

function [indhU,indhU2,indhX,indhX2,indhBi,indhself,n_psi,n_r] = find_hDCM(hA,hB,hC,hD)

n_psi = 0;
% hU    
ihU = find(hU~=0);
if ~isempty(ihU)
    indhU = 1:length(ihU);
else
    indhU = [];
end


%hU2
ihU2 = [];
indhU2 = cell(length(hU2),1);
for i=1:length(hU2)
    tmp = find(hU2{i}~=0);
    if ~isempty(tmp)
        indhU2{i} = length(ihU)+length(ihU2)+1:...
            length(ihU)+length(ihU2)+length(tmp);  
    else
        indhU2{i} = [];
    end
    ihU2 = [ihU2;tmp];
    
end
%hX   
ihX = find(hX~=0);
if ~isempty(ihX)
    indhX = length(ihU)+length(ihU2)+1:...
            length(ihU)+length(ihU2)+length(ihX);
else
    indhX = [];
end
n_r=max(n_r,size(hX,1));
%hX2
ihX2 = [];
indhX2 = cell(length(hX2),1);
for i=1:length(hX2)
    tmp = find(hX2{i}~=0);
    if ~isempty(tmp)
        indhX2{i} = length(ihU)+length(ihU2)+length(ihX)+1:...
                    length(ihU)+length(ihU2)+length(ihX)+length(tmp);
    else
        ihX2{i} = [];
    end
    ihX2 = [ihX2;tmp];
    n_r=max(n_r,size(hX2{i},1));
end
%hBi
ihBi = [];
indhBi = cell(length(hBi),1);
for i=1:length(hBi)
    tmp = find(hBi{i}~=0);
    if ~isempty(tmp)
        indhBi{i} = length(ihU)+length(ihU2)+length(ihX)+length(ihX2)+1:...
                    length(ihU)+length(ihU2)+length(ihX)+length(ihX2)+length(tmp);
    else
        ihBi{i} = [];
    end
    ihBi = [ihBi;tmp];
    n_r=max(n_r,size(hBi{i},1));
end
indhself = length(ihU)+length(ihU2)+length(ihX)+length(ihX2)+length(ihBi)+1;



end

function [dhU,dhU2,dhX,dhX2,dhBi] = get_dMatdvec(inH)
%hU
if ~isempty(inH.indhU)
    dhU = dMatdvec(inH.hU);
else
    dhU = [];
end
% hU2
dhU2 = cell(length(inH.hU2),1);
for i=1:length(inH.hU2)
    if ~isempty(inH.indhU2{i})
        dhU2{i} = dMatdvec(inH.hU2{i});
    else
        dhU2{i} = [];
    end
end
%hX
if ~isempty(inH.indhX)
    dhX = dMatdvec(inH.hX);
else
    dhX = [];
end
% hX2
dhX2 = cell(length(inH.hX2),1);
for i=1:length(inH.hX2)
    if ~isempty(inH.indhX2{i})
        dhX2{i} = dMatdvec(inH.hX2{i});
    else
        dhX2{i} = [];
    end
end
% hBi
dhBi = cell(length(inH.hBi),1);
for i=1:length(inH.hBi)
    if ~isempty(inH.indhBi{i})
        dhBi{i} = dMatdvec(inH.hBi{i});
    else
        dhBi{i} = [];
    end
end


end

function C = dMatdvec(A)
A = ~~A;
ind = find(A~=0);
n = numel(A);
ni = length(ind);
C = zeros(n,ni);
for i=1:length(ind)
    C(ind(i),i) = 1;
end
end

