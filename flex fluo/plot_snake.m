function cellOrder = plot_snake(input,sortTime,xtitle,zscore)
% % plot_snake %
%PURPOSE:   Plot selectivity, based on PSTHs from one condition (i.e. a
%           snake plot)
%AUTHORS:   AC Kwan 170720
%
%INPUT ARGUMENTS
%   input:        Structure generated by get_psth()
%   sortTime:     Sort the cells based on signal value in this time period
%   xtitle:       Text to put as the label for x-axis.
%   zscore:       If true, then the input.signal contains z-score values,
%                 not dF/F (no need to further normalize here)
%
%OUTPUT ARGUMENTS
%   cellOrder:    Sorting order of the cells
%                 (could be fed back as input 'sortParam' for next call of
%                 this function)

if nargin < 4
    zscore = false;
end

%% Plotting parameters
%colors=cbrewer('div','RdBu',256);
%colors=flipud(colors);

colors=cbrewer('seq','OrRd',256);

%% Setup
t=input{1}.t;
nCells=numel(input);

signal=[];
for j=1:nCells
    if (zscore)  %if z-score, then use z-score
        signal(:,j)=input{j}.signal;
    else         %otherwise, normalize to range of 0 to 1
        signal(:,j)=(input{j}.signal-nanmin(input{j}.signal))./(nanmax(input{j}.signal)-nanmin(input{j}.signal));
    end
end

% do an initial scramble first, so cells from the same experiments are not
% always adjacent to each other in our matrix
scramIdx = randperm(size(signal,2));
signal = signal(:,scramIdx);

% will sort the cells using signal from the specified time range
tIdx=[max([sum(t<=sortTime(1)) 1]):sum(t<=sortTime(2))];  %index should start from at least value of 1

% sort by time of signal peak, testing set
% [~,peakBin]=nanmax(signal(tIdx,:));
% [~,cellOrder]=sort(peakBin);

% sort by center of mass (mass should be all positive)
for j=1:nCells
    mass = signal(tIdx,j) - nanmin(signal(tIdx,j));
    com(j) = sum(t(tIdx).*mass)/sum(mass);
end
[~,cellOrder]=sort(com);

% set up the y-axis range
if (zscore)
    ampl = max([abs(prctile(signal(:),95)) abs(prctile(signal(:),5))]);  %maximum deviation from zero
    maxY = ampl;   %make the scale symmetric around zero
    minY = -ampl;
else   %else dF/F had been normalized to a range from 0 to 1
    maxY=prctile(signal(:),95);
    minY=prctile(signal(:),5);
end

%% Make the pseudo-color snake plot

% PSTH is only reliable if there are at least several events.. 
% set threshold at 10
if (input{1}.nEvent >= 10)
    figure;
    
    subplot(1,3,1);
    image(t,1:nCells,signal(:,cellOrder)','CDataMapping','scaled');
    hold on; 
    plot([0 0],[0 nCells+1],'w');
    colormap(colors);
    caxis([minY maxY]);      %normalize dF/F heatmap to max of all conditions
    ylabel('Cells');
    xlabel(xtitle);
    title(input{1}.psth_label);
    
    %make a color scale bar
    subplot(3,20,50);
    image(0,linspace(minY,maxY,100),linspace(minY,maxY,100)','CDataMapping','scaled');
    colormap(colors);
    if (zscore)
        caxis([minY maxY]);
        title('z-score');
    else
        caxis([0 1]);
        title('Normalized dF/F');
    end
    set(gca,'YDir','normal');
    set(gca,'XTick',[]);
    
    print(gcf,['snake_' input{1}.psth_label],'-dpng');    %png format
    saveas(gcf, ['snake_' input{1}.psth_label], 'fig');
    print(gcf,['snake_' input{1}.psth_label],'-depsc','-painters');   %eps format
end

end
