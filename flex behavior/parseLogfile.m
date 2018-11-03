function [ logData ] = parseLogfile( data_dir, logfile )
% % parseLogfile %
%
%PURPOSE: To read and parse Presentation logfile for further analysis.
%AUTHORS: MJ Siniscalchi & AC Kwan, 161209.
%
%INPUT ARGUMENTS
%   data_dir:    Path for logfile.   
%   logfile:     Filename for logfile.
%
%OUTPUT VARIABLES
%   logData:     Structure containing these fields:
%                {subject, dateTime, header, values}.

fID=fopen(fullfile(data_dir,logfile));  %use fullfile to avoid backslash/slash difference in mac/pc

header{1} = textscan(fID,'%s %*[^\n]',1);
header{2} = textscan(fID,'%s %s %c %s %s',1);
header{3} = textscan(fID,'%s %s %s %s %s %*[^\n]',1);  %Presentation column labels 

%ignore everything (e.g. lick response) that appears before the first Manual line in log file
header{4} = textscan(fID,'%s %s %s %s %s %s %s',1); 
while ~strcmp(header{4}{3},'Manual')                
    header{4} = textscan(fID,'%s %s %s %s %s %s %s',1);
end

data = textscan(fID,'%s %u %s %u %d %*[^\n]','delimiter','\t','EmptyValue',-1);
fclose(fID);

logData.subject = header{4}{1};
logData.dateTime = [header{2}{4:5}];
logData.header = [header{3}{1:5}]; %'Subject' 'Trial' 'Event Type' 'Code' 'Time'
logData.values = data;

end

