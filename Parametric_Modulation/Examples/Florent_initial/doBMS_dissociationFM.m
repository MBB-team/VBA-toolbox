% script to compute a bayesian model selection for 
% 	- the DISSOCIATION TASK
% 	- linear & hyperbolic models
% , by: 
%   - 1) calling the function to compute the bayesian model inversion and 
%        the log-evidence of each model and each subject
%   - 2) calling the function to save the log-evidence matrix to be use by 
%        the spm function
%   - 3) writing the spm batch
%   - 4) running the batch into spm
%   - 5) plotting the exceedence probability for each model and using the
%        right labels
%
% FOR PARALLEL COMPUTING
% Florent Meyniel 2012-05-31

clear all
% matlabpool(8)

% 1- COMPUTE MODEL LOG-EVIDENCE
% =============================
method='predifine bad model'; % 'predifine bad model' or 'bad defined on result'

% BMS: Disso
subList = [1:6,8:10,12,13,15:18]; 
model = 'linear'; % 'linear' 'hyperbolic'
type = 'FindOnOff';
[modev, ind] = computeLogEv_dissociation_v3_par(subList, type, model);

if strcmp(model, 'linear') 
    pathout = '/home/florent/ACCUMULATEUR/these/manip/ACCU-EXPLICITE-DISSOCIATION/results/DAVB/20120531/DissoLin_v3';
    mkdir(pathout)
    save(strcat(pathout, '/modev'), 'modev')
end

if strcmp(model, 'hyperbolic') 
    pathout = '/home/florent/ACCUMULATEUR/these/manip/ACCU-EXPLICITE-DISSOCIATION/results/DAVB/20120531/DissoHyp_v3';
    mkdir(pathout)
    save(strcat(pathout, '/modev'), 'modev')
end

% 2- COMPUTE BAYESIAN MODEL SELECTION
% ===================================
% save the F.mat subject x modexl matrix of log evidence

reshape_Fmat_fn(modev, pathout);

% create the matlab batch
bms_dcm.dir                        = {pathout};
bms_dcm.sess_dcm                   = {};
bms_dcm.model_sp                   = {''};
bms_dcm.load_f                     = {strcat(pathout, '/F.mat')};
bms_dcm.method                     = 'RFX';
bms_dcm.family_level.family_file   = {''};
bms_dcm.bma.bma_no                 = 0;
bms_dcm.verify_id                  = 0;

matlabbatch{1}.spm.stats.bms.bms_dcm = bms_dcm;

% save the matlab batch
save(strcat(pathout, '/batch.mat'), 'matlabbatch')

% run the matlab batch
spm('defaults', 'fmri');
spm_jobman('initcfg');
spm_jobman('run', strcat(pathout, '/batch.mat'));
spm('quit')

% 3- PLOT DATA
% ============
res   = load(strcat(pathout, '/BMS.mat'));
names = load(strcat(pathout, '/Names.mat'));
tmp = load(strcat(pathout, '/modev.mat'));
modev = tmp.modev;
nModels = length(unique(modev(:,2)));

figure(1)
subplot(1,2,1)
set(gcf, 'Color', [1 1 1])
set(gcf, 'Name', 'BMS results (with table of all models)')
models = modev(1:nModels, 3:11);
imagesc(models)
ylabel('Model #')
title('(white: allowed)')
set(gca, 'XTick', [1:9], 'XTickLabel', {'K_$', 'Je_$', 'Jr_$', 'K_F', 'Je_F', 'Jr_F', 'K_D', 'Je_D', 'Jr_D'})
colormap('bone')

subplot(1,2,2)
xp = NaN(nModels, 1);
xp(names.Names) = res.BMS.DCM.rfx.model.xp;
barh(xp)
set(gca, 'YDir', 'reverse')
axis([0 1 1 nModels])
xlabel({'Model exceedence'; 'probability'})
set(gcf, 'Color', [1 1 1])

figure(2)
subplot(1,2,1)
set(gcf, 'Color', [1 1 1])
set(gcf, 'Name', 'BMS results (with table of all models)')
ind = names.Names;
models = modev(ind, 3:11);
imagesc(models)
ylabel('Model #')
title('(white: allowed)')
set(gca, 'XTick', [1:9], 'XTickLabel', {'K_$', 'Je_$', 'Jr_$', 'K_F', 'Je_F', 'Jr_F', 'K_D', 'Je_D', 'Jr_D'})
colormap('bone')

subplot(1,2,2)
xp = res.BMS.DCM.rfx.model.xp;
barh(xp)
set(gca, 'YDir', 'reverse')
ylim([1 length(names.Names)])
xlabel({'Model exceedence'; 'probability'})
