function check_imageTriggerTimes ( imageTimes, logfileTimes, tlabel )
% % check_imageTriggerTimes %
%PURPOSE:   Plot the image file start times relative to behavioral trigger
%           times
%AUTHORS:   AC Kwan 171012
%
%INPUT ARGUMENTS
%   imageTimes:     Time stamps for start of each .tiff file
%   logfileTimes:   Time stamps for when Presentation sent TTL pulse to
%                   ScanImage to initiate a new .tiff file
%   tlabel:         Text to put in the title

%set first image frame and first logged trigger to time 0
%note: do not subtract imageTimes(1), which sometimes could have value of -1, 
%but that does not mean ScanImage started at t=-1 because stackInfo.frameTime(1)=0
if imageTimes(1) ~= -1
    imageTimes = imageTimes - imageTimes(1);
end
logfileTimes = logfileTimes - logfileTimes(1);

maxVal = max([max(logfileTimes) max(imageTimes)]);

figure; hold on;
plot([0 maxVal],[0 maxVal],'k','LineWidth',2);
if numel(imageTimes) == numel(logfileTimes)
    plot(logfileTimes/60,imageTimes/60,'r.','MarkerSize',20);
elseif numel(imageTimes) < numel(logfileTimes)
    plot(logfileTimes(1:numel(imageTimes))/60,imageTimes/60,'r.','MarkerSize',20);
    plot(logfileTimes(numel(imageTimes)+1:end)/60,zeros(size(logfileTimes(numel(imageTimes)+1:end))),'r.','MarkerSize',20);
elseif numel(imageTimes) > numel(logfileTimes)
    plot(logfileTimes/60,imageTimes(1:numel(logfileTimes))/60,'r.','MarkerSize',20);
    plot(zeros(size(imageTimes(numel(logfileTimes)+1:end))),imageTimes(numel(logfileTimes)+1:end)/60,'r.','MarkerSize',20);
end

xlim([0 maxVal/60]);
ylim([0 maxVal/60]);
axis square;
xlabel({'Time stamps for triggers according to .log file (min)';['[n = ' int2str(numel(logfileTimes)) ']']});
ylabel({'Time stamps for start of .tiff files (min)';['[n = ' int2str(numel(imageTimes)) ']']});
title({tlabel;'Value should not deviate from diagonal'},'interpreter','none');

print(gcf,'-dpng','check-triggertiming');    %png format
saveas(gcf, 'check-triggertiming', 'fig');

end