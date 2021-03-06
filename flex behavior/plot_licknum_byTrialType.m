function plot_licknum_byTrialType(input)
% % plot_licknum_byTrialType %
%PURPOSE:   Plot lick numbers for different trial types
%AUTHORS:   AC Kwan 180120
%
%INPUT ARGUMENTS
%   input:        Structure generated by get_lickrate_byTrialType().
%   tlabel:       Text to put as title of the plot.

%%
% if called from single-session analysis, input is a struct
% if called from summary analysis, input is a cell array
% here convert everything to a cell array first
if ~iscell(input)
    temp = input;
    clear input;
    input{1} = temp;
end

% load from the cell array
trialType=input{1}.trialType;
trialLabel=[trialType trialType];

% find mean number of licks for each trial type
for j=1:numel(input)
    for l=1:numel(trialType)
        up_leftNumLickL(l,j)=nanmean(input{j}.up_leftNumLickL{l});
        up_rightNumLickR(l,j)=nanmean(input{j}.up_rightNumLickR{l});
        down_leftNumLickL(l,j)=nanmean(input{j}.down_leftNumLickL{l});
        down_rightNumLickR(l,j)=nanmean(input{j}.down_rightNumLickR{l});
        
        %combining left and right spout --> target vs non-target spout
        NumLickTarget(l,j)=nanmean([input{j}.up_leftNumLickL{l}; input{j}.down_rightNumLickR{l}]);
        NumLickNontarget(l,j)=nanmean([input{j}.up_rightNumLickR{l}; input{j}.down_leftNumLickL{l}]);
        
        % --- for the control analysis, when we adjust for lick rate
%         if l==4     %for error column, plot for subset of hit trials with low lick rates to match the omitted-reward trials
%             temp = [input{j}.up_leftNumLickL{1}; input{j}.down_rightNumLickR{1}];
%             thresh = prctile(temp,30);
%             NumLickTarget(l,j) = nanmean(temp(temp<=thresh));
%             
%             temp = input{j}.up_leftNumLickL{1};
%             thresh = prctile(temp,30);
%             up_leftNumLickL(l,j) = nanmean(temp(temp<=thresh));
%             temp = input{j}.down_rightNumLickR{1};
%             thresh = prctile(temp,30);
%             down_rightNumLickR(l,j) = nanmean(temp(temp<=thresh));
%             
%             trialLabel{4} = 'hitsubset';
%         end 
    end
end

gray=[0.7 0.7 0.7];

%% plot
ymin = 0;
ymax = 160;

figure;

for k = 1:3
    if k==1
        subplot(2,1,1); hold on;
        temp = [NumLickTarget; NumLickNontarget];  %number of licks on trials       
        temp = 100* temp ./ repmat(temp(1,:),numel(trialType)*2,1);  %normalize as fraction to the first trial type
        col = 'k';
        tlabel = 'Licks on [target spout; non-target spout]';
        
    elseif k==2
        subplot(2,1,2); hold on;
        temp = [up_leftNumLickL; down_leftNumLickL];  %number of left licks on trials
        temp = 100* temp ./ repmat(temp(1,:),numel(trialType)*2,1);  %normalize as fraction to the first trial type
        col = 'r';
        tlabel = 'Left/right licks [target spout; non-target spout]';
    elseif k==3
        subplot(2,1,2); hold on;
        temp = [down_rightNumLickR; up_rightNumLickR];  %number of right licks on trials
        temp = 100* temp ./ repmat(temp(1,:),numel(trialType)*2,1);  %normalize as fraction to the first trial type
        col = 'b';
        tlabel = 'Left/right licks  [target spout; non-target spout]';
    end
        
    for l=1:size(temp,1)
        plot(l-0.25+0.5*rand(size(temp(l,:))),temp(l,:),'^','MarkerSize',10,'LineWidth',2,'Color',gray);
        if numel(input)>1  %more than 1 data set, plot mean+-sem
            plot(l+[-0.25 0.25],nanmean(temp(l,:))*[1 1],[col '-'],'LineWidth',3);
            plot(l+[0 0],nanmean(temp(l,:))+nanstd(temp(l,:))/sqrt(numel(temp(l,:)))*[-1 1],[col '-'],'LineWidth',3);
        end
        
       % --- print the mean and std; print the statistical test results
%       [nanmean(temp(l,:)) nanstd(temp(l,:))/sqrt(numel(temp(l,:)))]
%       signrank(temp(1,:),temp(l,:))
    end
    
    xlim([0 size(temp,1)+1]);
    ylim([ymin ymax]);
    set(gca,'xtick',1:size(temp,1));
    set(gca,'xticklabel',trialLabel);
    ylabel({'Licks per trial';'relative to single-reward (%)'});
    title(tlabel);
    
end

print(gcf,'-dpng','licknum_byTrialType');    %png format
saveas(gcf, 'licknum_byTrialType', 'fig');

end


