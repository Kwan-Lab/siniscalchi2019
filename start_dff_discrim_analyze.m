%% Analyze fluorescence data, for discrim tasks with no flexibility
% run this after running start_beh and start_dff_compute

clearvars;
close all;

tic;    %set clock to estimate how long this takes

runDecode = true;  %run the decoding analysis -- takes time to compute

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
        temp = sprintf('%04d',i);
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
        savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-fluo');
    else
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
        savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-fluo');
    end
    
    % load the saved behavioral analysis (from start_beh.m)
    % load the saved dF/F (from start_dff_compute.m)
    cd(savematpath);
    load('beh.mat');
    load('dff.mat');
    
    cd(savefluofigpath);
    
    %% Calculate trial-averaged dF/F and choice selectivity
    
    params=[];
    params.trigTime = trialData.cueTimes;
    params.xtitle = 'Time from stimulus (s)';
    params.window = [-2:0.25:6.5];
    params.minNumTrial = 5;
    
    for j=1:numel(cells.dFF)
        %calculate PSTH, i.e. trial-averaged dF/F
        fieldname={'sound','upsweep','left','hit'}; trialMask = getMask(trials,fieldname);
        psth_left{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        fieldname={'sound','downsweep','right','hit'}; trialMask = getMask(trials,fieldname);
        psth_right{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
    end
    
    save(fullfile(savematpath,'dff_and_beh.mat'),...
        'psth_left','psth_right');
    
    %% Calculate trial-averaged dF/F and choice selectivity
    
    params=[];
    params.trigTime = trialData.cueTimes;
    params.xtitle = 'Time from stimulus (s)';
    params.window = [-2:0.25:6.5];  %current trial
    params.minNumTrial = 5;
    
    for j=1:numel(cells.dFF)
        %calculate PSTH, i.e. trial-averaged dF/F
        fieldname={'sound','upsweep','left','hit'}; trialMask = getMask(trials,fieldname);
        psth_left{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        fieldname={'sound','downsweep','right','hit'}; trialMask = getMask(trials,fieldname);
        psth_right{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        
        fieldname={'sound','downsweep','left','err'}; trialMask = getMask(trials,fieldname);
        psth_left_err{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        fieldname={'sound','upsweep','right','err'}; trialMask = getMask(trials,fieldname);
        psth_right_err{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        
        fieldname={'sound','upsweep','left','omitreward'}; trialMask = getMask(trials,fieldname);
        psth_left_omit{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        fieldname={'sound','downsweep','right','omitreward'}; trialMask = getMask(trials,fieldname);
        psth_right_omit{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        
        fieldname={'sound','upsweep','left','doublereward'}; trialMask = getMask(trials,fieldname);
        psth_left_double{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        fieldname={'sound','downsweep','right','doublereward'}; trialMask = getMask(trials,fieldname);
        psth_right_double{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
        
        %calculate choice selectivity = (A-B)/(A+B)
        choicesel_hit{j} = calc_selectivity(psth_left{j},psth_right{j});
        choicesel_err{j} = calc_selectivity(psth_left_err{j},psth_right_err{j});
        choicesel_omit{j} = calc_selectivity(psth_left_omit{j},psth_right_omit{j});
        choicesel_double{j} = calc_selectivity(psth_left_double{j},psth_right_double{j});
    end
    
    save(fullfile(savematpath,'dff_and_beh.mat'),...
        'choicesel_hit','choicesel_err','choicesel_omit','choicesel_double','-append');
    
    %% CHOICE AND OUTCOME: Multiple linear regression  - choice and reward and their interaction
    
    params=[];
    
    %first predictor is choice; dummy-code: left=-1, right=1, miss=NaN
    params.choiceEvent=NaN(size(trials.left));
    params.choiceEvent(trials.left) = -1;
    params.choiceEvent(trials.right) = 1;
    %second predictor is outcome; dummy-code: reward=0, omit/error=-1, double-reward=1, miss=NaN
    params.outcomeEvent=NaN(size(trials.hit));
    params.outcomeEvent(trials.hit) = 0;  %rationale: for well-learned subjects, this is par on course
    params.outcomeEvent(trials.err) = NaN;  %well-learned animal, source of error unclear
    params.outcomeEvent(trials.omitreward) = -1;
    params.outcomeEvent(trials.doublereward) = 1;
    
    params.trigTime = trialData.cueTimes;
    params.xtitle = 'Time from stimulus (s)';
    params.window = [-2:0.5:6.5];
    params.nback = 2;       %how many trials back to regress against
    params.interaction = true; %consider interaction terms (our data do not have enough trials)
    params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
    
    %only perform analysis on trials with a response (when animal is engaged)
    fieldname={'left','right'};
    trialMask = getAnyMask(trials,fieldname);
    for j=1:numel(cells.dFF)
        reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, [params.choiceEvent params.outcomeEvent], params.trigTime, trialMask, params );
    end
    
    tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};
    plot_regr(reg_cr,params.pvalThresh,tlabel,params.xtitle);
    print(gcf,'-dpng','MLR-choiceoutcome');    %png format
    saveas(gcf, 'MLR-choiceoutcome', 'fig');
    
    save(fullfile(savematpath,'dff_and_beh.mat'),...
        'reg_cr','-append');
    
    %% which cell was choice selective
    
    params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
    
    timeIdx=sum(0>reg_cr{1}.regr_time);   %find index associated with time = 0 s
    
    nCells = numel(reg_cr);
    choice_mod=false(nCells,1);
    interaction_mod=false(nCells,1);
    for j=1:nCells
        %at least 5 bins were significant for n trial
        if sum(reg_cr{j}.pval(timeIdx:end,2)<params.pvalThresh) >= 5 %2 for C(n)
            choice_mod(j)=true;
        end
        if sum(reg_cr{j}.pval(timeIdx:end,8)<params.pvalThresh) >= 5 %8 for C(n)xR(n)
            interaction_mod(j)=true;
        end
    end
    
    mod = choice_mod | interaction_mod;  %cells selective for C(n) or C(n)xR(n)
    
    %% ensemble decoding analysis
    if (runDecode)
        % create array for ensemble activity [time x cell]
        dFF_ens=[];
        for j=1:numel(cells.dFF)
            dFF_ens(:,j) = cells.dFF{j};
        end
        
        %% Decoding CHOICE, linear classifier using every cell in the ensemble
        
        % predict choice; dummy-code: left=-1, right=1, miss=NaN
        params=[];
        params.trigEvent=NaN(size(trials.left));
        params.trigEvent(trials.left) = -1;
        params.trigEvent(trials.right) = 1;
        % construct linear classifier to decode choice, using dF/F from ensemble, use signals simultaneously recorded
        params.traintest_fieldname={'hit'};     % use these trials to train classifier and then test (x-fold cross-validation)
        params.addtest_fieldname{1}={'doublereward'};    % use the same classifier to additionally test these trials
        params.addtest_fieldname{2}={'omitreward'};      % use the same classifier to additionally test these trials
        params.addtest_fieldname{3}={'err'};    % use the same classifier to additionally test these trials
        
        %--- linear classifier
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from stimulus (s)';
        params.window = [-2:0.5:6.5];
        params.frac = 0.8;      %use this fraction to construct classifier, save the rest for testing, x-fold cross-validation
        params.numRep = 100;     %number of repeats, should be >=30
        params.nBack = 1;       %look at decoding of choice for trials n=0 and 1 back
        
        lclass_ens_choice=decode_linearclassifier( dFF_ens, cells.t, trials, params );
        tlabel='Linear classifier';
        plot_decode(lclass_ens_choice,params.xtitle,tlabel);
        
        %--- same settings, but use random forest classifier
        params.numTrees = 100;  %forest with 100 trees
        
        RF_ens_choice=decode_randomforest( dFF_ens, cells.t, trials, params );
        
        tlabel='Random forest';
        plot_decode(RF_ens_choice,params.xtitle,tlabel);
        
        %--- same settings, but linear classifier based on leave-one-out sampling
        params.frac = [];
        params.numRep = [];
        lclass_ens_choiceLOO=decode_linearclassifierLOO( dFF_ens, cells.t, trials, params );
        tlabel='Linear classifier (LOO)';
        plot_decode(lclass_ens_choiceLOO,params.xtitle,tlabel);
        
        save(fullfile(savematpath,'dff_decode.mat'),...
            'lclass_ens_choice','RF_ens_choice','lclass_ens_choiceLOO');
        
        %% Decoding CHOICE, adding one cell to ensemble at a time
        % --- add cell picked randomly
        numDraw = 100;    %set to low to test the code, for actual analysis should be >=30
        cellNumList = [1:1:30];  %characterize ensembles with these numbers of cells
        
        % predict choice; dummy-code: left=-1, right=1, miss=NaN
        params=[];
        params.trigEvent=NaN(size(trials.left));
        params.trigEvent(trials.left) = -1;
        params.trigEvent(trials.right) = 1;
        % construct linear classifier to decode choice, using dF/F from ensemble, use signals simultaneously recorded
        params.traintest_fieldname={'hit'};     % use these trials to train classifier and then test (x-fold cross-validation)
        
        %--- linear/random forest classifier
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from stimulus (s)';
        params.window = [2 4];  %but only look at one time interval
        params.frac = 0.8;      %use this fraction to construct classifier, save the rest for testing, x-fold cross-validation
        params.numRep = 1;      %we are making multiple draws, that's where the bootstrap comes from, so will only do 1 5-fold cross-validation for each draw
        params.numTrees = 100;  %forest with 100 trees
        
        lclass_subset_choice = []; RF_subset_choice = [];
        for j=1:numel(cellNumList)
            disp(['Ensemble decoding with #cell = ' int2str(cellNumList(j))]);
            for k=1:numDraw
                
                idx = randperm(numel(cells.dFF));
                dFF_ens_subset = dFF_ens(:,idx(1:cellNumList(j)));  %ensemble constructed from a few randomly drawn cells
                
                % --- linear classifier
                temp = decode_linearclassifier( dFF_ens_subset, cells.t, trials, params );
                lclass_subset_choice.corrPred(j,k) = nanmean(temp.corrPred);
                lclass_subset_choice.corrPred_randsig(j,k) = nanmean(temp.corrPred_randsig);
                lclass_subset_choice.corrPred_scram(j,k) = nanmean(temp.corrPred_scram);
                
                % --- random forest
                temp = decode_randomforest( dFF_ens_subset, cells.t, trials, params );
                RF_subset_choice.corrPred(j,k) = nanmean(temp.corrPred);
                RF_subset_choice.corrPred_randsig(j,k) = nanmean(temp.corrPred_randsig);
                RF_subset_choice.corrPred_scram(j,k) = nanmean(temp.corrPred_scram);
            end
            
            lclass_subset_choice.cellNum(j) = cellNumList(j);
            RF_subset_choice.cellNum(j) = cellNumList(j);
        end
        lclass_subset_choice.fieldname = temp.fieldname;
        lclass_subset_choice.decode_time = temp.decode_time;
        lclass_subset_choice.numDraw = numDraw;
        lclass_subset_choice.numRepeat = temp.numRepeat;
        RF_subset_choice.fieldname = temp.fieldname;
        RF_subset_choice.decode_time = temp.decode_time;
        RF_subset_choice.numDraw = numDraw;
        RF_subset_choice.numRepeat = temp.numRepeat;
        
        tlabel='C(n) - linear classifier';
        plot_decode_subset(lclass_subset_choice,tlabel);
        print(gcf,'-dpng',['lclass-subset-choice']);    %png format
        saveas(gcf, ['lclass-subset-choice'], 'fig');
        
        tlabel='C(n) - random forest';
        plot_decode_subset(RF_subset_choice,tlabel);
        print(gcf,'-dpng',['RF-subset-choice']);    %png format
        saveas(gcf, ['RF-subset-choice'], 'fig');
        
        save(fullfile(savematpath,'dff_decode.mat'),...
            'lclass_subset_choice','RF_subset_choice','-append');
        
        %% Pairwise correlation, for real and pseudo-ensemble
        
        %--- quantify, all cells
        
        params=[];        
        params.numRepeat = 1000;  %number of times to scramble
        
        params.trigEvent=NaN(size(trials.left)); %dummy-code: left/single-reward=-1, right/single-reward=1, miss=NaN
        params.trigEvent(trials.left & trials.hit) = -1;
        params.trigEvent(trials.right & trials.hit) = 1;
        
        params.trigTime = trialData.cueTimes;
        params.window = [2 4];  %but only look at one time interval
        
        pcorr = pairwiseCorr( dFF_ens, cells.t, params );
        plot_pcorr(pcorr);
        print(gcf,'pcorr','-dpng');    %png format
        saveas(gcf,'pcorr', 'fig');
        
        %--- quantify, using choice-selective cells
        
        pcorr_modCells = pairwiseCorr( dFF_ens(:,mod), cells.t, params );
        plot_pcorr(pcorr_modCells);
        print(gcf,'pcorr_modCells','-dpng');    %png format
        saveas(gcf,'pcorr_modCells', 'fig');
        
        save(fullfile(savematpath,'dff_decode.mat'),...
            'pcorr','pcorr_modCells','-append');
        
    end
    
    %%
    clearvars -except i dirs expData expParam runDecode;
    close all;
end

% plays sound when done
load train;
sound(y,Fs);

toc