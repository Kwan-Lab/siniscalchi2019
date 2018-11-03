%% Plot fluorescence data, for discrim tasks with no flexibility
% run this after running start_beh and start_dff_compute

clearvars;
close all;

tic;    %set clock to estimate how long this takes

%% setup path and plotting formats

flex_setPathList;

setup_figprop;  %set up default figure plotting parameters

%% which set of data to load?

data_subdir = fullfile(data_dir,'M2 omission');
[ dirs, expData ] = expData_M2_Omission(data_subdir);

%% process data files
%i=2; j=38; %cell example 1
i=5; j=38; %cell example 2
%i=1; j=9;  %cell example 3

disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);

% setup/create subdirectories to save analysis and figures
if isfield(expData(i),'onefolder')   %all the log files saved in one folder
    temp = sprintf('%04d',i);
    savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
    savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-fluo');
else
    savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
    savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-fluo');
end
if ~exist(savefluofigpath,'dir')
    mkdir(savefluofigpath);
end

% load the saved behavioral analysis (from start_beh.m)
cd(savematpath);
load('beh.mat');
load('dff.mat');

cd(savefluofigpath);

%% Plot cue-aligned dF/F for each cell
params=[];
params.trigTime = trialData.cueTimes;
params.xtitle = 'Time from sound cue (s)';
params.window = [-3:0.5:7];
params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
params.CI = 0.9;  %confidence interval
params.minNumTrial = 0; %only calc PSTH if there is this number of trials

%% plot different choices in same panels
psth_panel=[];
for k=1:3
    fieldname=[];
    if k==1 %panel 1
        fieldname{1}={'sound','upsweep','left','hit'}; col{1}='r'; linstyle{1}='-';
        fieldname{2}={'sound','downsweep','right','hit'}; col{2}='b'; linstyle{2}='-';
    elseif k==2 %panel 2
        fieldname{1}={'sound','upsweep','left','doublereward'}; col{1}='r'; linstyle{1}=':';
        fieldname{2}={'sound','downsweep','right','doublereward'}; col{2}='b'; linstyle{2}=':';
    elseif k==3 %panel 3
        fieldname{1}={'sound','upsweep','left','omitreward'}; col{1}='r'; linstyle{1}='--';
        fieldname{2}={'sound','downsweep','right','omitreward'}; col{2}='b'; linstyle{2}='--';
    end
    for kk=1:numel(fieldname)
        trialMask = getMask(trials,fieldname{kk});
        psth_panel(k).sig{kk} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
        psth_panel(k).col{kk} = col{kk};
        psth_panel(k).linstyle{kk} = linstyle{kk};
    end
end
tlabel = ['Cell ' int2str(j)];

plot_psth(psth_panel,tlabel,params.xtitle);
print(gcf,'-dpng',['cell' int2str(j) '-choice']);
saveas(gcf, ['cell' int2str(j) '-choice'], 'fig');

%% Plot different outcomes in same panels

psth_panel=[];
for k=1:3
    fieldname=[];
    if k==1 %panel 1
        fieldname{1}={'sound','upsweep','left','hit'}; col{1}='r'; linstyle{1}='-';
        fieldname{2}={'sound','upsweep','left','doublereward'}; col{2}='r'; linstyle{2}=':';
        fieldname{3}={'sound','upsweep','left','omitreward'}; col{3}='r'; linstyle{3}='--';
    elseif k==2 %panel 2
        fieldname{1}={'sound','downsweep','right','hit'}; col{1}='b'; linstyle{1}='-';
        fieldname{2}={'sound','downsweep','right','doublereward'}; col{2}='b'; linstyle{2}=':';
        fieldname{3}={'sound','downsweep','right','omitreward'}; col{3}='b'; linstyle{3}='--';
    elseif k==3 %panel 3
        fieldname{1}={'sound','hit'}; col{1}='k'; linstyle{1}='-';
        fieldname{2}={'sound','doublereward'}; col{2}='k'; linstyle{2}=':';
        fieldname{3}={'sound','omitreward'}; col{3}='k'; linstyle{3}='--';
    end
    for kk=1:numel(fieldname)
        trialMask = getMask(trials,fieldname{kk});
        psth_panel(k).sig{kk} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
        psth_panel(k).col{kk} = col{kk};
        psth_panel(k).linstyle{kk} = linstyle{kk};
    end
end
tlabel = ['Cell ' int2str(j)];

plot_psth(psth_panel,tlabel,params.xtitle);
print(gcf,'-dpng',['cell' int2str(j) '-outcome']);
saveas(gcf, ['cell' int2str(j) '-outcome'], 'fig');

% plays sound when done
load train;
sound(y,Fs);

toc