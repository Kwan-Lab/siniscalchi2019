function output = get_psth ( signal, t, eventTime, psth_label, params)
% % get_psth %
%PURPOSE:   Given a time-series signal and event times, find peri-event signal
%           by binning, and estimate CI using bootstrap
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   eventTime:      the event times
%   psth_label:     text saying what the event is about
%   params.window:  the time window around which to align signal
%   params.minNumTrial: only calculate PSTH if number of events exceed this
%                       number (too few trials then PSTH not accurate)
%
%   params.CI:                  confidence interval (CI)
%   params.numBootstrapRepeat:  number of bootstrap repeats to estimate CI
%
%   params.tSpon:       if animal has a quiescence period with no behavior,
%                       then we use this to convert dF/F to a z-score
%
%OUTPUT ARGUMENTS
%   output:        structure containing peri-event time histogram
%
% To plot the output, use plot_psth().

%%
if isfield(params,'numBootstrapRepeat')
    numRep=params.numBootstrapRepeat;
    CI=params.CI;
else
    numRep=0;
    CI=[];
end

if isfield(params,'tSpon');
    tempSig = [];
    for j = 1:size(params.tSpon,1)
        idx = find(t > params.tSpon(j,1) & t < params.tSpon(j,2));
        tempSig = [tempSig; signal(idx)];  %take all the dF/F values from the spontaneous activity period
    end
    tempSig = tempSig(~isnan(tempSig));    %get rid of time points with NaN value as signal
    
    sponMean = nanmean(tempSig);
    sponStd = nanstd(tempSig);
    
    signal = (signal - sponMean)/sponStd;   %convert the signal to a z-score
end

%set up the parameters
window=params.window;
window_dt=nanmean(diff(window));
output.t=window(1:end-1)+window_dt/2; %use center of bin as the time associated with the average signal
output.t=output.t(:);

if numel(eventTime) >= params.minNumTrial  %will do PSTH if at least this number of such event occurred
    
    %% find average (by binning/histogram)
    
    %find the segment of signal and its relative time, around each event
    tempSig=[]; tempTime=[]; tempEvent=[];
    nEvent = numel(eventTime);  %number of events that contributed to this PSTH; if event count too low, should be careful using this PSTH
    for j=1:numel(eventTime)
        relTime=t-eventTime(j);     %time relative to the event
        
        tempSig=[tempSig; signal(relTime>=window(1) & relTime<=window(end))];
        tempTime=[tempTime; relTime(relTime>=window(1) & relTime<=window(end))];
        tempEvent=[tempEvent; j*ones(sum(relTime>=window(1) & relTime<=window(end)),1)];    %the event# associated with this signal segment
        
        %mean value for entire window for each trial
        valbyTrial(j,1)=nanmean(signal(relTime>=window(1) & relTime<=window(end)));
    end
    for j=1:numel(window)-1 %for each bin of time, what is the average signal?
        psth(j,1)=nanmean(tempSig(tempTime>=window(j) & tempTime<(window(j)+window_dt)));
    end
    
    %% find CI using the bootstrap method
    
    if ~isempty(CI)
        bootSig=[];
        for i=1:numRep
            
            %draw a subset for bootstrap
            drawIndex=randsample(numel(eventTime),numel(eventTime),'true'); %each time draw a set of events with replacement
            drawSig=[]; drawTime=[];
            for k=1:numel(drawIndex)
                drawSig=[drawSig; tempSig(tempEvent==drawIndex(k))];
                drawTime=[drawTime; tempTime(tempEvent==drawIndex(k))];
            end
            
            %go through each time bin, and average
            for j=1:numel(window)-1
                bootSig(j,i)=nanmean(drawSig(drawTime>=window(j) & drawTime<(window(j)+window_dt)));
            end
        end
        
        %bootstrap mean and confidence interval
        output.bootavg=squeeze(nanmean(bootSig,2));
        output.bootlow=quantile(bootSig,0.5*(1-CI),2);
        output.boothigh=quantile(bootSig,1-0.5*(1-CI),2);
    end
    
else %if said event never occur, PSTH has no meaning
    psth = nan(size(output.t));
    valbyTrial = NaN;
    nEvent = 0;
end

output.psth_label=psth_label;
output.nEvent = nEvent;
output.valbyTrial = valbyTrial;
output.signal = psth;

end