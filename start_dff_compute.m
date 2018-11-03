%% Analyze fluorescence data -- compute dF/F
% run this after running start_beh.m

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
for i = 1:numel(expData)
    disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);
    
    % setup/create subdirectories to save analysis and figures
    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',i);     %then create subfolder with names based on the file orders
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
        savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-fluo');
    else                              %otherwise, one directory = one data set, create subfolder using the same directory name
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
        savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-fluo');
    end
    if ~exist(savefluofigpath,'dir')
        mkdir(savefluofigpath);
    end
    
    % load the saved behavioral analysis (from start_beh.m)
    cd(savematpath);
    load('beh.mat');
    
    % load the fluorescence data extracted using cellROI, calculate dF/F
    cd(fullfile(dirs.data,expData(i).sub_dir))
    stackInfo = load('stackinfo.mat');
    
    %% make an initial plot of the imaging and behavior trigger times to see
    %if there are obvious problems
    cd(savefluofigpath);
    tlabel = char(stackInfo.savFile_name);
    check_imageTriggerTimes ( stackInfo.trigTime, trialData.startTimes, tlabel );
    
    %% calculate the dF/F
    neuropilSubtractFactor = 0;  %neuropil subtraction
    
    cd(fullfile(dirs.data,expData(i).sub_dir));
    [ cells ] = calc_dFF( stackInfo, trialData.startTimes, sessionData.timeLastEvent, neuropilSubtractFactor );
    
    if isfield(cells,'dFF') %if there were any fluorescence data
        % check to see that imaging and behavior timing is sync'ed properly
        cd(savefluofigpath);
        check_imageFrameTimes ( cells, stackInfo, tlabel );
        
        % save the dF/F so don't have to re-compute every time
        save(fullfile(savematpath,'dff.mat'),...
            'cells');
    end
    
    %% save the dF/F so don't have to re-compute every time
    save(fullfile(savematpath,'dff.mat'),...
        'cells');
    
    clearvars -except i dirs expData;
end

% plays sound when done
load train;
sound(y,Fs);

toc