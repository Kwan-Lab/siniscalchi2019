function [ STIM, RESP, OUTCOME, EVENT ] = flex_getPresentationCodes(presCodeSet)
% % flex_getPresentationCodes %
%
%PURPOSE: To read and parse Presentation logfile for further analysis.
%AUTHORS: MJ Siniscalchi & AC Kwan, 161209.
%
%OUTPUT VARIABLES
%   presCodeSEt:    Allow toggling between different types of event code
%                   definitions
%
%OUTPUT VARIABLES
%   STIM:     fields containing stimulus-related eventCode defined in Presentation
%   RESP:     fields containing response-related eventCode defined in Presentation
%   OUTCOME:  fields containing outcome-related eventCode defined in Presentation
%   EVENT:    fields containing other event-related eventCode defined in Presentation

% for all discrim/flexibility log files prior to June 2017
if presCodeSet == 1
    STIM.sound_UPSWEEP=1;   %Rule/Sound cue combo: Sound/Upsweep
    STIM.sound_DNSWEEP=2;
    STIM.left_UPSWEEP=9;    %Rule/Sound cue combo: Action-Left/UpSweep
    STIM.left_DNSWEEP=10;
    STIM.right_UPSWEEP=11;  %Rule/Sound cue combo: Action-Right/UpSweep
    STIM.right_DNSWEEP=12;
    STIM.reversal_UPSWEEP=17;  %Rule/Sound cue combo: Action-Right/UpSweep
    STIM.reversal_DNSWEEP=16;
    
    RESP.LEFT=2;
    RESP.RIGHT=3;
    
    OUTCOME.REWARDLEFT=3;   %Following hit
    OUTCOME.REWARDRIGHT=8;
    OUTCOME.NOREWARD=4;     %aka Timeout event in Presentation
    OUTCOME.MISS=6;         %aka Pause event in Presentation
    OUTCOME.TWOREWARDLEFT=13;      %for probReward: 10% correct responses-->2x Reward
    OUTCOME.TWOREWARDRIGHT=14;
    OUTCOME.ZEROREWARD=15;         %for probReward: 10% correct responses-->0x Reward (OLD PRESENTATION CODE)
    OUTCOME.ZEROREWARDLEFT=16;     %for probReward: 10% correct responses-->0x Reward
    OUTCOME.ZEROREWARDRIGHT=17;
        
    EVENT.STARTEXPT = 5;
    EVENT.ENDEXPT = 7;      %inter-trial period
    EVENT.INTERPULSE=18;     %Presentation sends 2 pulses to valve to deliver double reward
    
elseif presCodeSet == 2
    STIM.sound_UPSWEEP=21;   %Rule/Sound cue combo: Sound/Upsweep
    STIM.sound_DNSWEEP=22;
    STIM.left_UPSWEEP=23;    %Rule/Sound cue combo: Action-Left/UpSweep
    STIM.left_DNSWEEP=24;
    STIM.right_UPSWEEP=25;  %Rule/Sound cue combo: Action-Right/UpSweep
    STIM.right_DNSWEEP=26;
    STIM.reversal_UPSWEEP=27;  %Rule/Sound cue combo: Action-Right/UpSweep
    STIM.reversal_DNSWEEP=28;
    
    RESP.LEFT=2;
    RESP.RIGHT=3;
    
    OUTCOME.REWARDLEFT=5;   %Following hit
    OUTCOME.REWARDRIGHT=6;
    OUTCOME.NOREWARD=7;     %aka Timeout event in Presentation
    OUTCOME.MISS=8;         %aka Pause event in Presentation
    
    EVENT.STARTEXPT = 4;
    EVENT.ENDEXPT = 9;      %inter-trial period
    
    %not used
    OUTCOME.TWOREWARDLEFT=999;     
    OUTCOME.TWOREWARDRIGHT=999;
    OUTCOME.ZEROREWARD=999;      
    OUTCOME.ZEROREWARDLEFT=999;       
    OUTCOME.ZEROREWARDRIGHT=999;

else
    disp('Warning: Code set is invalid for flex_getPresentationCodes');

end

