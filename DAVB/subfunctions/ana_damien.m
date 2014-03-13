% this script performs automatic statistical data analysis of Damien's data

clear all
close all

stefano = 0; % include RL-task?
zs = 1; % z-score?

if stefano
    ind = [5 6 7 65 107 113 115 117 118 126 129]; % indices of independent variables
    fn1 = ['C:\Users\JeanD\Documents\MATLAB\data\old',filesep,...
        'Tableau-Insight-StefanoPatientsTrackHD-2009-Motiv-DEF'];
else
    ind = [5 6 7 65 107 113 115 117 118 126]; % indices of independent variables
%     fn1 = ['C:\Users\JeanD\Documents\MATLAB\data',filesep,...
%         'Tableau-Insight-68PatientsTrackHD-2009-withIRM-apathie-mouvObsOnline-Motiv-DEF'];
    fn1 = ['C:\Users\JeanD\Documents\MATLAB\data',filesep,...
        'Tableau-Insight-68Patients-TrackHD-2009-DEF'];
end

i1 = 103; % score on-line
i2 = 48; % score off-line

[num,txt,raw] = xlsread(fn1,1,'A1:EG120');
% get column names
if zs
    names = txt(1,ind);
else
    names = cat(2,txt(1,ind),'MEAN');
end
na1 = txt(1,i1);
na2 = txt(1,i2);

% check missing data
X = num(:,ind-2);
sx = sum(X,2);
in = find(~isnan(sx));


% get data and z-score it
if zs
    X = zscore(X(in,:));
    y1 = zscore(num(in,i1-2));
    y2 = zscore(num(in,i2-2));
else
    X = [X(in,:),ones(size(in,1),1)];
    y1 = num(in,i1-2);
    y2 = num(in,i2-2);
end

% perform omnibus test on on-line score
n = size(X,2);
c = eye(n);
if zs
    c(:,end) = []; % remove mean from omnibus test
end
type = 'F';
[pv1,stat,df,all1] = GLM_contrast(X,y1,c,type,1,names);
set(gcf,'name',na1{1})


% augment model with on-line score and perform analysis of off-line score
X = [X,y1];
names{end+1} = na1{1};
n = size(X,2);
c = eye(n);
c(:,end) = []; % remove mean from omnibus test
type = 'F';
[pv2,stat,df,all2] = GLM_contrast(X,y2,c,type,1,names);
set(gcf,'name',na2{1})


