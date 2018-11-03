function [ SR ] = getSRmappings( blockSeq )
%%% getSRmappings() %%%
%PURPOSE: Retrieve contingencies for each block of the RuleSwitching task.
%AUTHOR: MJ Siniscalchi 161215.
%--------------------------------------------------------------------------
%
%INPUT ARGUMENTS
%   blockSeq:= Cell array of strings indicating sequence of rule blocks.
%OUTPUT VARIABLES
%   SR:= Struct containing these fields, characterizing
%           stimulus-response-trialType structure of each rule block:
%
%--------------------------------------------------------------------------
nBlocks = numel(blockSeq);
type = {'Sound' 'ActionL' 'ActionR' 'Reversal'};
if sum(ismember(blockSeq,type)) < nBlocks
    disp('ERROR in function getSRTrialTypes:');
    disp('Mismatch between number of types of contingencies.')
    return
end;

SR.mapping = cell(1,nBlocks);
for i = 1:numel(type)
    switch type{i}
        case 'Sound'
            tempCell = {'upsweep' 'downsweep'; 'left' 'right'}; %should be fieldnames for getTrialMasks
        case 'ActionL'
            tempCell = {'upsweep' 'downsweep'; 'left' 'left'};
        case 'ActionR'
            tempCell = {'upsweep' 'downsweep'; 'right' 'right'};
        case 'Reversal'
            tempCell = {'upsweep' 'downsweep'; 'right' 'left'};
    end; 
    blockMask = strcmp(blockSeq,type{i});    %Blocks in the sequence matching type{i}
    SR.mapping(blockMask) = {tempCell};      %Populate with {contingency strings}
end

end