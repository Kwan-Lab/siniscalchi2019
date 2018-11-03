function [ sessionData, trialData ] = flex_getSessionData( logData )
% % getSessionData %
%PURPOSE: Retrieve session data for flexibility task.
%AUTHORS: MJ Siniscalchi 161212.
%         modified by AC Kwan 170515
%
%INPUT ARGUMENTS
%   logdata:    Structure obtained with a call to parseLogfile().
%
%OUTPUT VARIABLES
%   sessionData:    Structure containing these fields:
%                   {subject, dateTime, nTrials, *lickTimes, *nSwitch}.
%                   * lickTimes([1 2]):=[left right] lick times.
%   trialData:      Fields:
%                   {startTimes, cueTimes, outcomeTimes, *cue, *response, *outcome}
%               *cue, response, and outcome for each trial contains
%               correspoinding eventCode from NBS Presentation

%% What event codes were used?

%COPY FROM LOGDATA
sessionData.subject = logData.subject;
sessionData.dateTime = logData.dateTime;

%SESSION DATA <<logData.header: 'Subject' 'Trial' 'Event Type' 'Code' 'Time'>>
TYPE = logData.values{3}; %Intersectional approach necessary, because values 2,3 were reused;
CODE = logData.values{4}; %change in future Presentation scripts: with unique codes, only CODE would be needed to parse the logfile...

tempidx=(strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')); %do not consider RESPONSE or MANUAL
codeUsed = unique(CODE(tempidx));         %List the set of event codes found in this logfile

% Get the event codes that are used in this log file, and then based on the set of codes used,
% make an educated guess on the correct event code set
% The reason for this section is that the set of event code (e.g. 6=reward;
% 5=miss, etc) was overhauled completely in around August, 2017; so older
% log files are stored based on one set of event codes, and newer log files
% may be stored based on another set of event codes
sessionData.presCodeSet = NaN;   %which is the event code set used?
for l = 1:2                      %test 2 potential event code sets
    [ STIM, RESP, OUTCOME, EVENT ] = flex_getPresentationCodes(l);
    cueCodes = cell2mat(struct2cell(STIM)); %get expected list of stimulus-associated codes as vector
    outcomeCodes = cell2mat(struct2cell(OUTCOME)); %get expected outcome-associated codes as vector
    eventCodes = cell2mat(struct2cell(EVENT));     %get expected event-associated codes as vector
    
    if sum(ismember(codeUsed,[cueCodes; outcomeCodes; eventCodes])==0)==0     %are all the used code a subset of the expected codes?
        sessionData.presCodeSet = l;
    end
end
if isnan(sessionData.presCodeSet)     %did not recognize the event code set
    disp('Codes that appeared in this log file');
    codeUsed
    error('ERROR in flex_getSessionData: there are unrecognized event codes found in the log file.');
else                    %otherwise, recognize the event code set and will use it!
    [ STIM, RESP, OUTCOME, EVENT ] = flex_getPresentationCodes(sessionData.presCodeSet);
    cueCodes = cell2mat(struct2cell(STIM)); %get expected list of stimulus-associated codes as vector
    outcomeCodes = cell2mat(struct2cell(OUTCOME)); %get expected outcome-associated codes as vector
end

% We will only consider trials that are complete. If the behavior was
% aborted in the middle of the very last trial, we will not consider those
% few orphan events
lastTrialEnd=find(strcmp(TYPE,'Nothing') & CODE==EVENT.ENDEXPT,1,'last'); %get rid of CODE beyond the last completed trial
TYPE = TYPE(1:lastTrialEnd);
CODE = CODE(1:lastTrialEnd);

%% Set up the time axis, and identify lick times
time_0 = logData.values{5}(find(CODE==EVENT.STARTEXPT,1,'first'));
if isempty(time_0)
    disp('ERROR in flex_getSessionData: there are no StartExpt event codes found in the log file.');
    disp('   check the expParam.presCodeSet used for flex_getPresentationCodes');
else
    time = logData.values{5}-time_0;   %time starts at first instance of startExpt
end
time = double(time)/10000;         %time as double in seconds
sessionData.timeLastEvent = time(end);   %time of the last event logged
sessionData.lickTimes{1} = time(strcmp(TYPE,'Response') & CODE==RESP.LEFT);    %left licktimes
sessionData.lickTimes{2} = time(strcmp(TYPE,'Response') & CODE==RESP.RIGHT);   %right licktimes

%% What are the events and when did they occur?
trialData.startTimes = time(strcmp(TYPE,'Nothing') & CODE==EVENT.STARTEXPT);

trialData.cue = CODE(strcmp(TYPE,'Sound') & ismember(CODE,cueCodes));
trialData.cueTimes = time(strcmp(TYPE,'Sound') & ismember(CODE,cueCodes));

trialData.outcome =  CODE((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')) & ismember(CODE,outcomeCodes));
trialData.outcomeTimes = time((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')) & ismember(CODE,outcomeCodes));

%% Check consistency of the extracted trials and trial times

% Each cue presented should be associated with a startExpt marker (time
% when a new .tiff imaging file is triggered)
if numel(trialData.startTimes) ~= numel(trialData.cueTimes)
    %Presentation sometimes seem to skip a few of the startExpt events?!
    disp('WARNING in flex_getSessionData: missing StartExpt events');
    disp('   so will estimate startExpt times = (cue times - 1)');
    
    startexpDur = [];
    for j = 1:numel(trialData.startTimes)
        temp = abs(trialData.cueTimes - trialData.startTimes(j));  %time of cues relative to a specific startExpt time
        temp = temp(temp>0);      %cue should follow startExpt
        startexpDur(j) = min(temp);   %identify the cue that is closest to startExpt
    end
    %use mode to identify the most frequency value of startexpDur
    %expect that to be the supposed duration of the startExpt event; this
    %trick works as long as most of the time startExpt event was present,
    %and we are only missing some events
    trialData.startTimes = trialData.cueTimes - mode(startexpDur);
end

% Each cue presented should be associated with a certain outcome
if numel(trialData.cueTimes) > numel(trialData.outcomeTimes)
    %In early 2014, miss events do not have event codes
    disp('WARNING in flex_getSessionData: some events missing outcomes');
    disp('   assuming Miss trials were not assigned');
    
    %for each cue, there should be an outcome within 2 sec, otherwise we
    %will categorize that as a miss
    temp_outcomeTimes=[]; temp_outcome=[];
    for j=1:numel(trialData.cueTimes)
        
        difftime = trialData.outcomeTimes - trialData.cueTimes(j);
        idx = (difftime > 0) & (difftime < 2);
        if sum(idx) == 1     %if such outcome time exist
            temp_outcomeTimes(j,1) = trialData.outcomeTimes(idx);
            temp_outcome(j,1) = trialData.outcome(idx);
        elseif sum(idx) == 0     %if no outcome, then it is a miss trial
            temp_outcomeTimes(j,1) = trialData.cueTimes(j) + 2;
            temp_outcome(j,1) = OUTCOME.MISS;
        else
            disp('Error when trying to fix cueTimes > outcomeTimes');
        end
    end
    trialData.outcomeTimes = temp_outcomeTimes;
    trialData.outcome = uint32(temp_outcome);
end

% Each trial must have one and only one outcome
sessionData.nTrials = numel(trialData.outcomeTimes);

%% When and what type of responses?
respIdx = find(strcmp(TYPE,'Response'));
respTimes = time(respIdx);

if nanmin(trialData.outcomeTimes - trialData.cueTimes) >= 0.4
    %if all outcomes occur much later than the cue presentation time, it is
    %likely that we set a grace period (to ignore early, impulsive licks)
    sessionData.gracePeriodDur = 0.5;  %grace period of 0.5 s was employed for many early studies, especially when animals start to learn, in the lab
else
    sessionData.gracePeriodDur = 0;    %no grace period
end

trialData.reaction = nan(sessionData.nTrials,1);             %direction of the first lick after cue
trialData.reactionTimes = nan(sessionData.nTrials,1);          %time of the first lick after cue (only for trial without pre-cue licks)

trialData.response = zeros(sessionData.nTrials,1,'uint32');  %direction of the first lick during response period
trialData.responseTimes = nan(sessionData.nTrials,1);          %time of the first lick during response period

trialData.preCueLickRate = nan(sessionData.nTrials,1);

idx = find(trialData.outcome~=OUTCOME.MISS); %Idx all non-miss trials
for i = 1:numel(idx)
    %Reaction time
    temp = find(respTimes>=trialData.cueTimes(idx(i)),1,'first');
    trialData.reaction(idx(i)) = CODE(respIdx(temp));
    % if there is no lick for 0.5 s before cue onset (i.e., mouse not licking prior to cue)
    if sum((respTimes-respTimes(temp))>-0.5 & (respTimes-respTimes(temp))<0)
        trialData.reactionTimes(idx(i)) = respTimes(temp)-trialData.cueTimes(idx(i));   %relative to cue time
    end
    
    %Time of the lick that counts as animal's choice (recall there is a 0.5 sec grace period)
    temp = find(respTimes>=(trialData.cueTimes(idx(i))+sessionData.gracePeriodDur),1,'first'); %first lick within response window
    trialData.response(idx(i)) = CODE(respIdx(temp));
    trialData.responseTimes(idx(i)) = respTimes(temp);   %in absolute terms
    
    %Pre-cue lick rate (with 1 s before cue), for trial where there was an
    %animal response (otherwise unfair because animals could be disengaged
    %from task, no pre-cue licking but no responses either)
    temp = find(respTimes>(trialData.cueTimes(idx(i))-1) & respTimes<=(trialData.cueTimes(idx(i))));
    trialData.preCueLickRate(idx(i)) = numel(temp);      %number of licks in the pre-cue window
    
end

% What are the lick times associated with each trial?
for i = 1:sessionData.nTrials
    %Times of all licks from prior trial's cue to next trial's cue
    if i>1
        time1=trialData.cueTimes(i-1);
    else
        time1=0;
    end
    
    if i<sessionData.nTrials
        time2=trialData.cueTimes(i+1);
    else
        time2=time(end);
    end
    
    temp=sessionData.lickTimes{1};
    temp=temp(temp>=time1 & temp<=time2); %save only those licks within range
    trialData.leftlickTimes{i}=temp'-trialData.cueTimes(i); %make lick time relative to cue
    
    temp=sessionData.lickTimes{2};
    temp=temp(temp>=time1 & temp<=time2); %save only those licks within range
    trialData.rightlickTimes{i}=temp'-trialData.cueTimes(i); %make lick time relative to cue
end

% What are the number of licks associated with each response? 
% (within 5 s of cue, avoid counting pre-cue anticipatory licks)
for i = 1:sessionData.nTrials
    %Times of all licks within 4 s of cue
    time1=trialData.cueTimes(i);
    time2=trialData.cueTimes(i) + 5;
    
    temp=sessionData.lickTimes{1};
    trialData.numLeftLick(i,1)=sum(temp>=time1 & temp<=time2); %save only those licks within range

    temp=sessionData.lickTimes{2};
    trialData.numRightLick(i,1)=sum(temp>=time1 & temp<=time2); %save only those licks within range
end

end

