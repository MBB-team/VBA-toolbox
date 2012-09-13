% This script is a example of group level analysis.

close all
% clear all
  
Nsubject= 10;
lme  = [randn(Nsubject,1),   randn(Nsubject,1)+4/Nsubject, randn(Nsubject,1)+8/Nsubject];
% Row : subjects
% Columns : models
model_names = {'bad','ok','good'};
dbstop if error
partition = {1,[2,3]};
[Lf] = models2families(lme,partition)

group_level_analysis(Lf,'RFX')  
% group_level_analysis(lme,'FFX')  
% group_level_analysis(lme,'FFX',model_names)  
  

%%

Lm = rand(10,8);
partition = {[1,3],[5,7]};
[Lf] = models2families(Lm,partition);
  


