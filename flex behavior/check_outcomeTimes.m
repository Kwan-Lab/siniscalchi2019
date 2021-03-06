function check_outcomeTimes ( sessionData, trialData )
% % check_outcomeTimes %
%PURPOSE:   Check the timing for when animal made chocie, and when reward
%           is delivered. Difference in time stamps for these events should 
%           be the zero. Sometimes multiple pulses are used to calibrate
%           and deliver precise rewards. Depending on which pulse is
%           associated as the "outcome event", the time stamps could be
%           misaligned.
%AUTHORS:   AC Kwan 170519
%
%INPUT ARGUMENTS
%   sessionData:    Structure generated by flex_getSessionData().
%   trialData:      Structure generated by flex_getSessionData().
%

%% 

tlabel=[char(sessionData.subject) ' - ' char(sessionData.dateTime{1})];

figure;
subplot(2,1,1); hold on;
plot(trialData.outcomeTimes- trialData.responseTimes,'k-');
xlabel('Trial');
ylabel({'Reward time - response time (s)'});
title({tlabel;'Value should not deviate from zero - otherwise perhaps >1 pulse was used to deliver water'},'interpreter','none');
ylim([-0.2 0.2]);

print(gcf,'-dpng','check-timingoutcome');    %png format
saveas(gcf, 'check-timingoutcome', 'fig');

end