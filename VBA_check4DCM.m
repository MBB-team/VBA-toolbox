function [options] = VBA_check4DCM(options)

A = options.inF.A;
B = options.inF.B;
C = options.inF.C;
D = options.inF.D;
TR = options.inF.deltat.*options.decim;
microDT = TR./options.decim;
if isfield(options.inG,'homogeneous')
    homogeneous = options.inG.homogeneous;
else
    homogeneous = 0;
end
[odcm] = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
inF0 = odcm.inF;
inG0 = odcm.inG;

[options.inF] = fillin(options.inF,inF0);
[options.inG] = fillin(options.inG,inG0);


function [in] = fillin(in,in0)
% fill in default DCM options if missing
if iscell(in0)
    for i=1:numel(in0)
        if numel(in) < i
            in{i} = fillin(in{i},in0{i});
        else
            in{i} = in0{i};
        end
    end
elseif isstruct(in0)
    fn0 = fieldnames(in0);
    fn = fieldnames(in);
    for i=1:length(fn0)
        if ~ismember(fn0{i},fn)
            in.(fn0{i}) = in0.(fn0{i});
        else
            in.(fn0{i}) = fillin(in.(fn0{i}),in0.(fn0{i}));
        end
    end
end




