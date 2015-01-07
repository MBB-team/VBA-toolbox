function [ u, dim, options ] = setInput(options,dim,varargin )
% allow to specify labels to inputs (u)
% IN
%     options: a structure
%     varargin: should be organized as: my_n1_times_nTrial_maxtrix1, name1, my_nk_times_nTrial_maxtrixk, namek...
% OUT
%     u : a concatenate (n1+...+nk) * nTrial matrix
%     options is filled with options.in[FG].inputLabel.(each name) fields indexing the lines of u

n = numel(varargin);
assert(mod(n,2)==0,'***Parameters should be: name1, ''matrix1'', name2, ''matrix2'',...');
[~, nTrial]=size(varargin{2});
for i=1:2:n
    if islogical(varargin{i+1})
        varargin{i+1} = +varargin{i+1};
    end
    assert(isnumeric(varargin{i+1}), '***Parameters should be: name1, ''matrix1'', name2, ''matrix2'',...')
    assert(ischar(varargin{i}), '***Parameters should be: name1, ''matrix1'', name2, ''matrix2'',...')
end


u =[];
cpt = 0;
for i=1:2:n
    newU = varargin{i+1};
    [nL, nT] = size(newU);
    assert(nT == nTrial, 'all u should have the same number of columns (i.e. nTrials)')
    u = [u; varargin{i+1}];
    inputLabel.(varargin{i})=(cpt+1):(cpt+nL);
    cpt=cpt+nL;
end

options.inF.inputLabel=inputLabel;
options.inG.inputLabel=inputLabel;


dim.u = size(u,1) ;

end

