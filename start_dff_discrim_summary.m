%% Summary plots from multiple analyses of imaging sessions
% run this after running start_dff_discrim_analyze on each session

clearvars;
close all;

%% setup path and plotting formats

flex_setPathList;

setup_figprop;  %set up default figure plotting parameters

disp('Summary of behavior + fluoresence data:');

%% which set of data to load?

data_subdir = fullfile(data_dir,'M2 omission');
[ dirs, expData ] = expData_M2_Omission(data_subdir);

%% collect analysis results
psth_left_array=[];
choicesel_hit_array=[]; choicesel_err_array=[]; choicesel_omit_array=[]; choicesel_double_array=[];
reg_cr_array=[]; 
lclass_ens_choice_array=[]; lclass_ens_choiceLOO_array=[]; lclass_subset_choice_array=[];
RF_ens_choice_array=[]; RF_subset_choice_array=[];
pcorr_array = []; pcorr_modCells_array=[];

%%
for i = 1:numel(expData)
    % setup/create subdirectories to save analysis and figures
    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',i);
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
    else
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
    end
    cd(savematpath);
    
    load('dff_and_beh.mat');
    load('dff_decode.mat');
    
    nCell(i) = numel(reg_cr);
    
    psth_left_array=[psth_left_array psth_left];
    
    choicesel_hit_array=[choicesel_hit_array choicesel_hit];
    choicesel_err_array=[choicesel_err_array choicesel_err];
    choicesel_omit_array=[choicesel_omit_array choicesel_omit];
    choicesel_double_array=[choicesel_double_array choicesel_double];
    
    reg_cr_array=[reg_cr_array reg_cr];
    
    lclass_ens_choice_array{i}=lclass_ens_choice;
    lclass_ens_choiceLOO_array{i}=lclass_ens_choiceLOO;
    lclass_subset_choice_array{i}=lclass_subset_choice;
    
    RF_ens_choice_array{i}=RF_ens_choice;
    RF_subset_choice_array{i}=RF_subset_choice;
    
    pcorr_array{i} = pcorr;
    pcorr_modCells_array{i} = pcorr_modCells;
    
end
nCells=numel(reg_cr_array);

%% make summary plots
savefluofigpath = fullfile(dirs.summary,'figs-fluo');
if ~exist(savefluofigpath,'dir')
    mkdir(savefluofigpath);
end
cd(savefluofigpath);

%% Make snake plots, all data

params.xtitle = 'Time from stimulus (s)';
plot_snake(psth_left_array,[0 6.5],params.xtitle);

%% Multiple linear regression  - choice and reward for sound trials

params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
params.xtitle = 'Time from stimulus (s)';

tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n) x R(n)','C(n-1) x R(n-1)','C(n-2) x R(n-2)'};
plot_regr(reg_cr_array,params.pvalThresh,tlabel,params.xtitle);
print(gcf,'MLR-choiceoutcome','-dpng');    %png format
saveas(gcf,'MLR-choiceoutcome', 'fig');
print(gcf,'MLR-choiceoutcome','-depsc','-painters');   %eps format

%% Choice selectivity

% find the cells with dF/F significantly modulated by C(n), according to regression analysis
params.pvalThresh = 0.01;   %p-value for coefficient be considered significant

%% comparing C(n) and C(n-1)
currtimeIdx=sum(2>reg_cr_array{1}.regr_time);   %find index associated with time = 2 s
pasttimeIdx=sum(0>reg_cr_array{1}.regr_time);   %find index associated with time = 0 s

curr_choice_mod=false(nCells,1);
past_choice_mod=false(nCells,1);
curr_choice_coeff=nan(nCells,1);
past_choice_coeff=nan(nCells,1);
for j=1:nCells
    if (reg_cr_array{j}.pval(currtimeIdx,2)<params.pvalThresh) %2 for C(n), t = 2 s
        curr_choice_mod(j)=true;
        curr_choice_coeff(j)=reg_cr_array{j}.coeff(currtimeIdx,2);
    end
    if (reg_cr_array{j}.pval(pasttimeIdx,3)<params.pvalThresh) %5 for C(n-1), t = 0 s
        past_choice_mod(j)=true;
        past_choice_coeff(j)=reg_cr_array{j}.coeff(currtimeIdx,3);
    end
end

figure;
subplot(2,3,1); hold on;
bar(1,100*sum(curr_choice_mod & ~past_choice_mod)/numel(curr_choice_mod),'k');
bar(2,100*sum(~curr_choice_mod & past_choice_mod)/numel(curr_choice_mod),'k');
bar(3,100*sum(curr_choice_mod & past_choice_mod)/numel(curr_choice_mod),'k');
ylabel('Fraction of cells (%)');
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'C(n)','C(n-1)','C(n) &'});
xlim([0 4]);
ylim([0 25]);

subplot(2,3,2); hold on;
plot(curr_choice_coeff,past_choice_coeff,'k.','MarkerSize',20);
plot([-1 1],[0 0],'k-','LineWidth',2);
plot([0 0],[-1 1],'k-','LineWidth',2);
xlabel('a_1 for C(n)');
ylabel('a_2 for C(n-1)');
axis([-0.25 0.25 -0.25 0.25]);  %this crops out 5 cells
axis square;

idx=~isnan(curr_choice_coeff) & ~isnan(past_choice_coeff);
cc=corrcoef(curr_choice_coeff(idx),past_choice_coeff(idx));
disp(['Correlation coefficient C(n), C(n-1): ' num2str(cc(1,2))]);

print(gcf,'lreg-currpastchoice','-dpng');    %png format
saveas(gcf,'lreg-currpastchoice', 'fig');
print(gcf,'lreg-currpastchoice','-depsc','-painters');   %eps format

%%
timeIdx(1)=sum(0>reg_cr_array{1}.regr_time);   %find index associated with time = 0 s

choice_mod=false(nCells,1);
choice_any_mod=false(nCells,1);
outcome_mod=false(nCells,1);
interaction_mod=false(nCells,1);
for j=1:nCells
    %at least 5 bins were significant for n trial
    if sum(reg_cr_array{j}.pval(timeIdx(1):end,2)<params.pvalThresh) >= 5 %2 for C(n)
        choice_mod(j)=true;
    end
    if sum(reg_cr_array{j}.pval(timeIdx(1):end,5)<params.pvalThresh) >= 5 %5 for R(n)
        outcome_mod(j)=true;
    end
    if sum(reg_cr_array{j}.pval(timeIdx(1):end,8)<params.pvalThresh) >= 5 %8 for C(n)xR(n)
        interaction_mod(j)=true;
    end
end

disp('Number of cells with significant p-value in multiple linear regression:');
disp(['   for C(n):' int2str(sum(choice_mod))]);
disp(['   for R(n):' int2str(sum(outcome_mod))]);
disp(['   for C(n)xR(n):' int2str(sum(interaction_mod))]);
disp(['   for C(n) or C(n)xR(n):' int2str(sum(choice_mod | interaction_mod))]);
disp(['   for (C(n) or C(n)xR(n)) AND R(n):' int2str(sum((choice_mod | interaction_mod) & outcome_mod))]);

mod = choice_mod | interaction_mod;  %for following analysis, we'll use cells selective for C(n) or C(n)xR(n)

%% if there were task-modulated cells, plot them
if sum(mod)>0
    %% Make snake plots, choice-selective cells
    params.xtitle = 'Time from stimulus (s)';
    plot_snake(psth_left_array(mod),[0 7],params.xtitle);

    %% Choice selectivity plot
    params.sortParam=[0 6.5];  %sort the cells based on amplitude of selectivity in this time period    
    params.colorRange=[-0.8 0.8];

    tlabel = ['Choice selectivity (n=' int2str(sum(mod)) ' choice-selective cells)'];
    cellOrder = plot_selectivity(choicesel_hit_array(mod),params.sortParam,tlabel,params.xtitle,params.colorRange);
    print(gcf,'choice-selectivity','-dpng');    %png format
    saveas(gcf,'choice-selectivity', 'fig');
    print(gcf,'choice-selectivity','-depsc','-painters');   %eps format
    
    %%
    compTime=[2 4];  %sort the cells based on selectivity in this time period
    minNumTrial = 1;
    tlabel = 'Choice selectivity';
    plot_psth_comp(choicesel_omit_array(mod),choicesel_double_array(mod),compTime,minNumTrial,tlabel);
    print(gcf,'choice-selectivity_omitdouble_comp','-dpng');    %png format
    saveas(gcf,'choice-selectivity_omitdouble_comp', 'fig');
    print(gcf,'choice-selectivity_omitdouble_comp','-depsc','-painters');   %eps format

end

%% Linear classifier decoding
tlabel='Linear classifier - random sub-sampling';
plot_decode(lclass_ens_choice_array,params.xtitle,tlabel);

tlabel='Linear classifier - LOO';
plot_decode(lclass_ens_choiceLOO_array,params.xtitle,tlabel);

%% Random forest
tlabel='Random forest';
plot_decode(RF_ens_choice_array,params.xtitle,tlabel);

%% Decoding accuracy for ensemble of diffrent sizes
tlabel='C(n) - linear classifier';
plot_decode_subset(lclass_subset_choice_array,tlabel);
print(gcf,'lclass-subset-choice','-dpng');    %png format
saveas(gcf,'lclass-subset-choice', 'fig');
print(gcf,'lclass-subset-choice','-depsc','-painters');   %eps format

tlabel='C(n) - random forest';
plot_decode_subset(RF_subset_choice_array,tlabel);
print(gcf,'RF-subset-choice','-dpng');    %png format
saveas(gcf,'RF-subset-choice', 'fig');
print(gcf,'RF-subset-choice','-depsc','-painters');   %eps format

%% Pairwise correlation
plot_pcorr(pcorr_array);
print(gcf,'pcorr','-dpng');    %png format
saveas(gcf,'pcorr', 'fig');

plot_pcorr(pcorr_modCells_array);
print(gcf,'pcorr_modCells','-dpng');    %png format
saveas(gcf,'pcorr_modCells', 'fig');

%% plays sound when done
load train;
sound(y,Fs);