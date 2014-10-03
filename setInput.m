function [ u, options ] = setInput(options,varargin )
% allow to specify labels to inputs (u)
% IN
%     options: a structure
%     varargin: should be organized as: my_n1_times_nTrial_maxtrix1, name1, my_nk_times_nTrial_maxtrixk, namek...
% OUT
%     u : a concatenate (n1+...+nk) * nTrial matrix
%     options is filled with options.in[FG].inputLabel.(each name) fields indexing the lines of u

n = numel(varargin);
assert(mod(n,2)==0,'***Parameters should be: matrix1, ''name1'', matrix2, ''name2'',...');
[~, nTrial]=size(varargin{1});
for i=1:2:n
    assert(isnumeric(varargin{i}), '***Parameters should be: matrix1, ''name1'', matrix2, ''name2'',...')
    assert(ischar(varargin{i+1}), '***Parameters should be: matrix1, ''name1'', matrix2, ''name2'',...')
end


u =[];
cpt = 0;
for i=1:2:n
    
    newU = varargin{i};
    [nL, nT] = size(newU);
    assert(nT == nTrial, 'all u should have the same number of columns (i.e. nTrials)')
    u = [u; varargin{i}];
    inputLabel.(varargin{i+1})=(cpt+1):(cpt+nL);
    cpt=cpt+nL;
end

options.inF.inputLabel=inputLabel;
options.inG.inputLabel=inputLabel;
end

