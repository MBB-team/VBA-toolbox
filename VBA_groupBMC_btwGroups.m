function [h,p] = VBA_groupBMC_btwGroups(Ls,options)
% test for between-groups difference in model frequencies
% function [h,p] = VBA_groupBMC_btwGroups(Ls,options)
% IN:
%   - Ls: {nmXns_1, nmXns_2} array of log-model evidences matrices of each group (nm models; ns_g subjects in the group).
%   - options: a structure containing the following fields:
%       .DisplayWin: flag for display window
%       .verbose: flag for summary statistics display
%       .families: a cell array of size nf, which contains the indices of
%       the models that belong to each of the nf families.
% OUT:
%   - h: test decision about a difference between the group (rejection of
%        the null hypothesis of equality)
%   - p: the posterior probability that the two groups have the same model
%        frequencies

if nargin < 2
    options = {};
end
% one frequency for all
L = [Ls{1} Ls{2}];
[posterior,out] = VBA_groupBMC(L,options);
Fe = out.F(end);

% separate frequencies
[posterior1,out1] = VBA_groupBMC(Ls{1},options);
[posterior2,out2] = VBA_groupBMC(Ls{2},options);
Fd = out1.F(end) + out2.F(end);

p = 1./(1+exp(Fd-Fe));
h = p<.05;

end