%% Summary plots from multiple analyses of behavioral performance
% run this after running start_beh

clearvars;
close all;

%% setup path and plotting formats

flex_setPathList;

setup_figprop;  %set up default figure plotting parameters

disp('Summary of behavioral data:');

%% which set of data to load?

data_subdir = fullfile(data_dir,'M2 omission');
[ dirs, expData ] = expData_M2_Omission(data_subdir);

%% collect analysis results
for i = 1:numel(expData)

    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',i);
        cd(fullfile(dirs.analysis,expData(i).sub_dir,temp));
    else
        cd(fullfile(dirs.analysis,expData(i).sub_dir));
    end
    load('beh.mat');

    choice_stat_array{i}=choice_stat;
    lick_trType_array{i}=lick_trType;
    
    edges = [0:1:30];
    histLick_missNext(:,i) = histc(numLick_missNext,edges);
    histLick_respNext(:,i) = histc(numLick_respNext,edges);
    histLick_missNext1(:,i) = histc(numLick_missNext1,edges);
    histLick_respNext1(:,i) = histc(numLick_respNext1,edges);
    histLick_missNext2(:,i) = histc(numLick_missNext2,edges);
    histLick_respNext2(:,i) = histc(numLick_respNext2,edges);

    close all;
    clearvars -except i dirs expData ...
        choice_stat_array lick_trType_array... 
        edges histLick_missNext histLick_respNext...
        histLick_missNext1 histLick_respNext1 histLick_missNext2 histLick_respNext2;
end

%% make summary plots
savebehfigpath = fullfile(dirs.summary,'figs-beh');
if ~exist(savebehfigpath,'dir')
    mkdir(savebehfigpath);
end

cd(savebehfigpath);
tlabel=strcat('Group summary, n=',int2str(numel(expData)));

%% summarize all sessions
plot_choice_stats(choice_stat_array,tlabel);

%%
plot_lickrate_overlay(lick_trType_array);

%%
plot_licknum_byTrialType(lick_trType_array);

%%
figure;
subplot(2,3,1); hold on; title('All');
bar(edges,nansum(histLick_missNext,2),'histc');
xlabel('# lick when next trial is Miss'); xlim([0 30]);
ylabel('Occurrence');
subplot(2,3,4); hold on;
bar(edges,nansum(histLick_respNext,2),'histc');
xlabel('# lick when next trial is Not Miss'); xlim([0 30]);
ylabel('Occurrence');
subplot(2,3,2); hold on; title('First half of sessions');
bar(edges,nansum(histLick_missNext1,2),'histc');
xlabel('# lick when next trial is Miss'); xlim([0 30]);
ylabel('Occurrence');
subplot(2,3,5); hold on;
bar(edges,nansum(histLick_respNext1,2),'histc');
xlabel('# lick when next trial is Not Miss'); xlim([0 30]);
ylabel('Occurrence');
subplot(2,3,3); hold on; title('Second half of sessions');
bar(edges,nansum(histLick_missNext2,2),'histc');
xlabel('# lick when next trial is Miss'); xlim([0 30]);
ylabel('Occurrence');
subplot(2,3,6); hold on;
bar(edges,nansum(histLick_respNext2,2),'histc');
xlabel('# lick when next trial is Not Miss'); xlim([0 30]);
ylabel('Occurrence');

%% plays sound when done
load train;
sound(y,Fs);