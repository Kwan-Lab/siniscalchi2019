function [ dirs, expData ] = expData_M2_Omission(data_dir)
% % expData_M2_Omission %
%
%PURPOSE: Create data structure for imaging tiff files and behavioral log files
%AUTHORS: AC Kwan, 170519.
%
%INPUT ARGUMENTS
%   data_dir:    The base directory to which the raw data are stored.  
%
%OUTPUT VARIABLES
%   dirs:        The subfolder structure within data_dir to work with
%   expData:     Info regarding each experiment

dirs.data = fullfile(data_dir,'data');
dirs.analysis = fullfile(data_dir,'analysis');
dirs.summary = fullfile(data_dir,'summary');

i=1;
        expData(i).sub_dir = '140603 L4';
        expData(i).logfile = 'L4-phase2b_auditory_discrim_IMG.log';
i=i+1;
        expData(i).sub_dir = '140608 L2';
        expData(i).logfile = 'L2-phase2b_auditory_discrim_IMG.log';
i=i+1;
        expData(i).sub_dir = '140608 M6';
        expData(i).logfile = 'M6-phase2b_auditory_discrim_IMG.log';
i=i+1;
        expData(i).sub_dir = '140613 M6';
        expData(i).logfile = 'M6-phase2b_auditory_discrim_IMG3.log';
i=i+1;
        expData(i).sub_dir = '140616 M6';
        expData(i).logfile = 'M6-phase2b_auditory_discrim_IMG4.log';

i=i+1;
        expData(i).sub_dir = '150820 M12';  
        expData(i).logfile = 'M12-phase2b_discrim_modReward.log';
i=i+1;
        expData(i).sub_dir = '150821 M20';  
        expData(i).logfile = 'M20-phase2b_discrim_modReward1.log';
i=i+1;
        expData(i).sub_dir = '150822 M12';  
        expData(i).logfile = 'M12-phase2b_discrim_modReward1.log';
i=i+1;
        expData(i).sub_dir = '150822 M13';  
        expData(i).logfile = 'M13-phase2b_discrim_modReward1.log';
i=i+1;
        expData(i).sub_dir = '150827 M12';  
        expData(i).logfile = 'M12-phase2b_discrim_varReward.log';
i=i+1;
        expData(i).sub_dir = '150828 M14';  
        expData(i).logfile = 'M14-phase2b_discrim_varReward.log';
i=i+1;
        expData(i).sub_dir = '150903 M16';  
        expData(i).logfile = 'M16-phase2b_discrim_varReward.log';
i=i+1;
        expData(i).sub_dir = '150904 M16';  
        expData(i).logfile = 'M16-phase2b_discrim_varReward1.log';
i=i+1;
        expData(i).sub_dir = '150910 M17'; 
        expData(i).logfile = 'M17-phase2b_discrim_varReward.log';
i=i+1;
        expData(i).sub_dir = '150930 M17';
        expData(i).logfile = 'M17-phase2b_discrim_varReward5.log';
i=i+1;
        expData(i).sub_dir = '150930 M22';
        expData(i).logfile = 'M22-phase2b_discrim_varReward2.log';
        
%         expData(i).sub_dir = '140602 M8'; %excluded, overall hit rate = 77% <80%
%         expData(i).logfile = 'M8-phase2b_auditory_discrim_IMG.log';
%
%         expData(i).sub_dir = '150929 M22'; %excluded, overall hit rate = 76% <80%
%         expData(i).logfile = 'M22-phase2b_discrim_varReward1.log';
%
%         expData(i).sub_dir = '150925 M17';   %fewer than 5 omit-reward + right trials
%         expData(i).logfile = 'M17-phase2b_discrim_varReward2.log';
% 
%         expData(i).sub_dir = '150820 M20';  %fewer than 5 double-reward + right trials
%         expData(i).logfile = 'M20-phase2b_discrim_modReward.log';
%
%         expData(i).sub_dir = '150821 M14';  %23/20 double-reward vs 7/6 omitted reward..?? does not seem random
%         expData(i).logfile = 'M14-phase2b_discrim_modReward1.log';
% 
%         expData(i).sub_dir = '150819 M14';  %motion correction problem -- strange dip for all cells in fluorescence at time of choice
%         expData(i).logfile = 'M14-phase2b_discrim_modReward.log';

        
