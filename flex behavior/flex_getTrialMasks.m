function [ trials ] = flex_getTrialMasks( trialData, blocks )
% % getTrialMasks %
%PURPOSE:   Create data structure, 'trials', containing logical masks
%           of size(nTrials,1) for task variables.
%AUTHORS:   MJ Siniscalchi 161214.
%           modified by AC Kwan 170515
%
%INPUT ARGUMENTS
%   trialData:  Structure generated by flex_getSessionData()
%   blocks:     Structure generated by flex_getBlockData()
%
%OUTPUT VARIABLES
%   trials:     Structure containing these fields, each a logical mask
%               indicating whether trial(idx) is of the corresponding subset, e.g.,
%               response==left or cue==upsweep.

%%
nTrials = numel(trialData.outcomeTimes);

%% GET CODES FROM PRESENTATION
[STIM,RESP,OUTCOME,EVENT] = flex_getPresentationCodes(blocks.presCodeSet);

%% GET MASKS FOR THOSE RESP/OUTCOME/RULE TYPES WITH CLEAR MAPPINGS
taskVar = {'cue' 'response' 'outcome' 'rule'};

for i = 1:numel(taskVar)
    clear codes;
    switch taskVar{i}
        case 'cue'
            codes.upsweep = [STIM.sound_UPSWEEP,...
                                STIM.left_UPSWEEP,...
                                STIM.right_UPSWEEP,...
                                STIM.reversal_UPSWEEP];
            codes.downsweep = [STIM.sound_DNSWEEP,...
                                STIM.left_DNSWEEP,...
                                STIM.right_DNSWEEP,...
                                STIM.reversal_DNSWEEP];
        case 'response'
            codes.left = [RESP.LEFT];
            codes.right = [RESP.RIGHT];
        case 'outcome'
            codes.hit = [OUTCOME.REWARDLEFT,...
                        OUTCOME.REWARDRIGHT];
            codes.err = [OUTCOME.NOREWARD];                        
            codes.miss = [OUTCOME.MISS];
            
            codes.doublereward = [OUTCOME.TWOREWARDLEFT,OUTCOME.TWOREWARDRIGHT];
            codes.omitreward = [OUTCOME.ZEROREWARD,OUTCOME.ZEROREWARDLEFT,OUTCOME.ZEROREWARDRIGHT];
        case 'rule'
            codes.sound   = [STIM.sound_UPSWEEP,...
                            STIM.sound_DNSWEEP];
            codes.actionL = [STIM.left_UPSWEEP,...
                            STIM.left_DNSWEEP];
            codes.actionR = [STIM.right_UPSWEEP,...
                            STIM.right_DNSWEEP];
            codes.reversal = [STIM.reversal_UPSWEEP,...
                            STIM.reversal_DNSWEEP];
            trialData.rule = trialData.cue; %Rule info is multiplexed in presentation codes for cue
    end;
    fields = fieldnames(codes);
    for j = 1:numel(fields)
        trials.(fields{j}) = ismember(trialData.(taskVar{i}),codes.(fields{j})); %Generate trial mask for each field in 'codes'
    end
end

%% GET MASKS FOR PERSEVERATIVE AND OTHER ERRORS
trials.persev = false(nTrials,1);   %perseverative choices (both hits and errors)
trials.err_p = false(nTrials,1);    %perseverative errors
trials.err_o = false(nTrials,1);    %non-perseverative errors
SR=getSRmappings(blocks.contingType);
nBlocks = numel(blocks.contingType);

%first block, classify all errors as other errors
blockMask = flex_getBlockMask(1, blocks);   %trials corresponding to block 1
trials.err_o(trials.err & blockMask) = true;

%beyond the first block
for i = 2:nBlocks
    
    %find which cue-responses consistent with prior block's rules
    tempMask=[];
    for j=1:2   % two mappings always (2 stim->resp contingencies)
        stimField=SR.mapping{i-1}{1,j};
        respField=SR.mapping{i-1}{2,j};
        tempMask(:,j) = trials.(stimField) & trials.(respField);
    end
    persevMask=tempMask(:,1) | tempMask(:,2);   %matches either of the past contingencies
    
    %find which trials belong to current block
    blockMask = flex_getBlockMask(i, blocks);
    
    trials.persev(persevMask & blockMask) = true;
    trials.err_p(persevMask & trials.err & blockMask) = true;
    trials.err_o(~persevMask & trials.err & blockMask) = true;
    
end

%% GET MASK FOR LAST 20 TRIALS PRE-SWITCH
trials.last20 = false(nTrials,1);   %perseverative choices (both hits and errors)
for i = 2:nBlocks
    trials.last20(blocks.firstTrial(i)-20:blocks.firstTrial(i)-1)=true;
end

%% GET MASK FOR N-TRIAL BACK

trials.hit_1=[false; trials.hit(1:end-1)];    %last trial was a reward?
trials.hit_2=[false; false; trials.hit(1:end-2)]; %two trials ago was a reward?

trials.err_1=[false; trials.err(1:end-1)];    
trials.err_2=[false; false; trials.err(1:end-2)]; 

trials.left_1=[false; trials.left(1:end-1)]; 
trials.left_2=[false; false; trials.left(1:end-2)];

trials.right_1=[false; trials.right(1:end-1)]; 
trials.right_2=[false; false; trials.right(1:end-2)];

trials.upsweep_1=[false; trials.upsweep(1:end-1)]; 
trials.upsweep_2=[false; false; trials.upsweep(1:end-2)];

trials.downsweep_1=[false; trials.downsweep(1:end-1)]; 
trials.downsweep_2=[false; false; trials.downsweep(1:end-2)];

%% GET MASK FOR REACTION (first lick after cue, not necessarily in response window)

trials.left_reaction = (trialData.reaction == RESP.LEFT);
trials.right_reaction = (trialData.reaction == RESP.RIGHT);

%% SPLIT SOUND CONTINGENCIES --> DOES IT PRECEDE ACTIONL OR ACTIONRIGHT?

trials.soundafAL=false(nTrials,1);   %sound contingency trials, after a prior action-left block
trials.soundafAR=false(nTrials,1);   %sound contingency trials, after a prior action-right block

for i = 2:nBlocks
    if strcmp(blocks.contingType{i},'Sound')    %if this is a sound block
        
        %find which trials belong to current block
        blockMask = false(nTrials,1);                   %Create mask for block(i)
        firstTrial = blocks.firstTrial(i);
        n = blocks.nTrials(i);
        blockMask(firstTrial:firstTrial+n-1) = true;
        
        if strcmp(blocks.contingType{i-1},'ActionL') %if it is after an action-left block
            trials.soundafAL(blockMask) = true;
        elseif strcmp(blocks.contingType{i-1},'ActionR')
            trials.soundafAR(blockMask) = true;
        end
    end
end

%% rules and contingencies masks

% Rule: sound-guided=1; non-conditional=2
trials.rule=nan(nTrials,1);
trials.rule(trials.sound)=1;
trials.rule(trials.actionL)=2;
trials.rule(trials.actionR)=2;
trials.rule(trials.reversal)=1;
trials.rule_labels = {'Sound-guided','Non-conditional'};

% Contingencies: sound=1; action-left=2; action-right=3; reversal=4
trials.conting=nan(nTrials,1);
trials.conting(trials.sound)=1;
trials.conting(trials.actionL)=2;
trials.conting(trials.actionR)=3;
trials.conting(trials.reversal)=4;
trials.conting_labels = {'Sound','ActionL','ActionR','Reversal'};

%% Check consistency among the extracted trial values
if sum(trials.hit)+sum(trials.err)+sum(trials.miss)+...
        sum(trials.doublereward)+sum(trials.omitreward) ~= numel(trials.hit)
    disp('ERROR in flex_getTrialMasks: check #1');
end
if sum(trials.upsweep)+sum(trials.downsweep) ~= numel(trials.hit)
    disp('ERROR in flex_getTrialMasks: check #2');
end
if sum(trials.hit)+sum(trials.err)+...
        sum(trials.doublereward)+sum(trials.omitreward) ~= sum(trials.left)+sum(trials.right)
    disp('ERROR in flex_getTrialMasks: check #3');
end
if sum(trials.err_p)+sum(trials.err_o) ~= sum(trials.err)
    disp('ERROR in flex_getTrialMasks: check #4');
end
if sum(trials.sound & ((trials.upsweep & trials.left) | (trials.downsweep & trials.right))) ...
        ~= sum(trials.sound & trials.hit)+sum(trials.doublereward)+sum(trials.omitreward)
    disp('ERROR in flex_getTrialMasks: check #5');
end