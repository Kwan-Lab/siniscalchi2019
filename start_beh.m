%% Analyze behavioral performance of flexibility task

clearvars;
close all;

tic;    %set clock to estimate how long this takes

%% setup path and plotting formats

flex_setPathList;

setup_figprop;  %set up default figure plotting parameters

%% which set of data to load?

data_subdir = fullfile(data_dir,'M2 omission');
[ dirs, expData ] = expData_M2_Omission(data_subdir);

suppressFigs = true; %do not plot anything (only save the analysis results)

%% process data files
for i = 1:numel(expData)
    disp(['--- Processing file ' int2str(i) '/' int2str(numel(expData)) '.']);
    disp([expData(i).sub_dir ' - ' expData(i).logfile]);
    
    % setup/create subdirectories to save analysis and figures
    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',expData(i).logfilenum);
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
        savebehfigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-beh');
    else
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
        savebehfigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-beh');
    end
    
    if ~exist(savematpath,'dir')
        mkdir(savematpath);
    end
    if ~exist(savebehfigpath,'dir')
        mkdir(savebehfigpath);
    end
    
    % parse the presentation log file
    [ logData ] = parseLogfile( fullfile(dirs.data,expData(i).sub_dir), expData(i).logfile );
    
    % break the parsed data into trials
    [ sessionData, trialData ] = flex_getSessionData( logData );
    
    % identify the rule blocks
    [ blocks ] = flex_getBlockData( sessionData, trialData );
    
    % identify the trial types/characteristics
    [ trials ] = flex_getTrialMasks( trialData, blocks );
    
    % identify the block types/characteristics
    [ blocks ] = flex_getMoreBlockData( blocks, trials );
    
    % save the analysis so can summarize the population results later
    save(fullfile(savematpath,'beh.mat'),...
        'logData','trialData','sessionData','blocks','trials');
    
    %% analysis and plotting of behavioral performance
    
    cd(savebehfigpath);
    tlabel=strcat('Subject=',char(logData.subject),', Time=',char(logData.dateTime(1)),'-',char(logData.dateTime(2)));
    
    % check whether reward and response times are consistent
    % (presentation scripts may use multiple pulses to deliver rewards.
    % i.e., introduce differences between response times (first pulse) and
    % reward time (event starts after last pulse)
    if ~exist('suppressFigs','var')
        check_outcomeTimes ( sessionData, trialData );
    end
    
    %% plot behavior in raw format
    time_range=[-2 6];
    
    if ~exist('suppressFigs','var')
        %        plot_session_beh_vert(trialData,trials,blocks,tlabel,time_range);
        plot_session_beh_horz(trials,blocks,tlabel);
    end
    
    
    %% plot choice behavior (makes sense if performance is stable across entire session and no flexibility is involved)
    if numel(blocks.firstTrial)==1   %if animal only perform discrim
        % plot key statistics
        choice_stat=choice_stats(trials);
        if ~exist('suppressFigs','var')
            plot_choice_stats(choice_stat,tlabel);
        end
        
        if sum(trials.omitreward)>0 %are there are omitted reward trials?
            trialType={'hit','doublereward','omitreward','err'};
        else
            trialType={'hit','err'};
        end
        edges=[-3:0.1:7];   % edges to plot the lick rate histogram
        lick_trType=get_lickrate_byTrialType(trialData,trials,trialType,edges);
        if ~exist('suppressFigs','var')
            plot_lickrate_byTrialType(lick_trType);
            plot_lickrate_overlay(lick_trType);
            plot_licknum_byTrialType(lick_trType);
        end
        
        %how well does current trial's number of licks predict whether next trial is a miss
        numLick = (trialData.numLeftLick(1:end-1)+trialData.numRightLick(1:end-1));  %number of licks in trials
        missNext = trials.miss(2:end);  %if next trial is a miss
        idxTrial = ~trials.miss(1:end-1);    %trials in which the mouse responded
        
        numLick_missNext = numLick(missNext & idxTrial);  %number of licks in trials when mouse misses next trial
        numLick_respNext = numLick(~missNext & idxTrial); %number of licks in trials when mouse responds next trial
        
        cutoff = round(3/4*sum(~trials.miss));
        numLick_missNext1 = numLick(missNext & [idxTrial(1:cutoff); false(size(idxTrial(cutoff+1:end)))]);  %number of licks in trials when mouse misses next trial
        numLick_respNext1 = numLick(~missNext & [idxTrial(1:cutoff); false(size(idxTrial(cutoff+1:end)))]); %number of licks in trials when mouse responds next trial
        numLick_missNext2 = numLick(missNext & [false(size(idxTrial(1:cutoff))); idxTrial(cutoff+1:end)]);  %number of licks in trials when mouse misses next trial
        numLick_respNext2 = numLick(~missNext & [false(size(idxTrial(1:cutoff))); idxTrial(cutoff+1:end)]); %number of licks in trials when mouse responds next trial
        
        % save the analysis so can summarize the population results later
        save(fullfile(savematpath,'beh.mat'),...
            'choice_stat','lick_trType','numLick_missNext','numLick_respNext',...
            'numLick_missNext1','numLick_respNext1','numLick_missNext2','numLick_respNext2','-append');
    end
        
    %%
    close all;
    clearvars -except i suppressFigs dirs expData;
end

% plays sound when done
load train;
sound(y,Fs);

toc