%% Command Line PLS setup for functional connectivity data
% This code was written by Laetitia Mwilambwe-Tshilobo July 20,2020 to
% run partial least-squares (PLS) analysis on connectivity data (fMRI). If
% you are unfamiliar with pldcmd outputs, a detailed description has been
% provided at the end of this script under "General Notes for interpreting
% PLS output"      

% Adapted from scripts provided by John Anderson, PhD
% Contact: laetitia.mwilambwe-tshilobo@mail.mcgill.ca
%% clear all workspace variables
clc;
clear all;
close all;
%% Setup directories and paths
% We need to set and add paths to functions to run command line PLS and to
% the directory to the data (correlation matrices, behavior).

cmdPLSDir = '/path/to/Pls/plscmd/' % Set pls directory
dataDir = '/path/to/plsConnectivityMatrices'; %Set data directory

% ----- add paths
addpath(genpath(cmdPLSDir));
addpath(genpath(dataDir)); 
%% Setup

plsThresh = 1.96 ; %bootstrap threshold (95%CI = 1.96; 99%CI = 2.98)
roiIndx = [1:400] ; %number of roi of matrices 
conditions = {'cond1'}; % total number of conditions
subjs_group1 = 181; %total number of subjects in group 1
subjs_group2 = 120; %total number of subjects in group 2
nsubj = (subjs_group1 + subjs_group2); %total number of subject in analysis 
num_perm = 1000; %number of permutations to perform for PLS
num_boot = 1000; %number of bootstraps to perform for PLS
%% Import correlation matrices
% Import data from Group 1 for all conditions. Data for each condition 
% should be formated as a 3D matrix (roiIndx,roiIndx,1:subjs_group).

% Group 1 condition 1 data 
load([dataDir,'/data/Cond1_group1/Cond1_group1.mat']); %path to data for group1 condition 1
Group1_1Z = Group1_cond1_xyz;

% Group 2 condition 1 data 
load([dataDir,'/data/Cond1_group2/Cond1_group2.mat']); %path to data for group2 condition 1
Group2_1Z = Group2_cond1_xyz;

% If your analysis involved more groups and/or conditions, input data 
% should follow a similar format as above (e.g. commented section below) 

% %Group 1 conditions 2 data
% load([dataDir,'/data/Cond2_group1/Cond2_group1.mat']); %path to data for group1 condition 2
% Group1_2Z = Group1_cond2_xyz;
% 
% %Group 2 conditions 2 data
% load([dataDir,'/data/Cond2_group2/Cond2_group2.mat']); %path to data for group2 condition 2
% Group2_2Z = Group2_cond2_xyz
%% Stack correlation matrices
% For each group, we need to vectorize the lower triangle matrix. We do 
% this by calling on the LowerTriangleIndex function to stack the lower
% triangle matrix into a row. For each group/condition, the rows of the 
% output matrix (Group#_#Stack)correspond to the connectivity data for each
% subject.

% make index for lower triangle matrix 
lowTriagDataIndx = LowerTriangleIndex(length(roiIndx));

% Group 1 condition 1
 for i = 1:subjs_group1
     g1_temp1 = Group1_1Z(:,:,i);
     Group1_1Stack(i,:)= g1_temp1(lowTriagDataIndx)';  
 end

 % Group 2 condition 1
 for i = 1:subjs_group2
     g2_temp1 = Group2_1Z(:,:,i);
     Group2_1Stack(i,:)= g2_temp1(lowTriagDataIndx)';  
 end

% If you have more than 1 condition per group, that's fine. The section 
% below provides an example of how you can format the script. 
 
 % Group 1 condition 2
%  for i = 1:subjs_group1
%      g1_temp2 = Group1_2Z(:,:,i);
%      Group1_2Stack(i,:)= g1_temp2(lowTriagDataIndx)';  
%  end
%   % Group 2 condition 2
%  for i = 1:subjs_group2
%      g2_temp2 = Group2_2Z(:,:,i);
%      Group2_2Stack(i,:)= g2_temp2(lowTriagDataIndx)';  
%  end

% making input data for command line pls
datamat_list{1} = cat(1, Group1_1Stack);% for more conditions do cat(Group1_1Stack,Group1_2Stack,Group1_3Stack);
datamat_list{2} = cat(1, Group2_1Stack);% for more conditions do cat(Group2_1Stack,Group2_2Stack,Group2_3Stack);
%% Behavior PLS analysis setup
% If you are interested in running a behavioral PLS, use this section to
% load your behavioral data. Please note that the data has to be organized 
% in the same subject order as the correlation matrices!!

% Add path to behavioral data. In the example below we are using MMSE and
% Age. Change these variables to your behaviors of interest. 
% load([dataDir,'data/MMSE.mat'])
% load([dataDir,'data/Age.mat'])
% behav = [MMSE, Age];
% behav = zscore(behav); 
%% Specify other parameter for pls analysis
%PLS analysis can be set upi in various ways. Here I provide code for two
%analysis types (non-rotated and rotated pls). For more help on what these
%means reference the plscmd script "pls_analysis.m". Use this section to 
%set the  PLS method that the program will use:
%			1. Mean-Centering Task PLS
%			2. Non-Rotated Task PLS
%			3. Regular Behavior PLS
%			4. Regular Multiblock PLS
%			5. Non-Rotated Behavior PLS
%			6. Non-Rotated Multiblock PLS


% switch AnalysisType
%     case 1
%         %%%%% non-rotated PLS (predefined contrasts) %%%%%
%         disp('Running PLS')
%         num_subj_lst = [subjs_group1 subjs_group2];
%         num_cond = 1;
%         option.method = 2; %non-rotated PLS
%         option.num_perm = num_perm;
%         option.num_boot = num_boot;
%         option.clim = 95;
%         option.stacked_designdata = contrast;
%         resultNonRotated = pls_analysis(datamat_list, num_subj_lst, num_cond,option);
%        
                    %
%    case 2
        %%%%% rotated PLS (no predefined contrast) %%%%%
        disp('Running PLS')
        num_subj_lst = [subjs_group1 subjs_group2];
        num_cond = numel(conditions);
        option.method = 1; %rotated PLS
        option.num_perm = num_perm; 
        option.num_boot = num_boot; 
        %option.stacked_behavdata = behav; %uncomment for behavioral PLS analysis
        option.clim = 95;
        %option.meancentering_type = 1; %removes grand mean (over all subjects and conditions)
        resultRotated = pls_analysis(datamat_list, num_subj_lst, num_cond,option);
%% Looking at results

% Plot p-values
pval = resultRotated.perm_result.sprob
nLV=numel(pval);
figure;
bar(pval,'r');
hold on;
h = zeros(nLV,1);
for i=1:nLV
    h(i)=plot(NaN,NaN, '.r');
end
legend(h,strcat('LV', num2str([1:nLV]'), {' - '} ,num2str(pval)));
title(['Permuted values greater than observed, ', num2str(option.num_perm), ' permutation tests']);
hold off;

% Plot effect sizes (% crossblock covariance)
pcov = resultRotated.s.^2 / sum(resultRotated.s.^2)
figure;
bar(pcov);
hold on;
h = zeros(nLV,1);
for i=1:nLV
    h(i)=plot(NaN,NaN, '.');
end
legend(h,strcat('LV', num2str([1:nLV]'), {' - '} ,num2str(pcov*100), '%'));
title('Percent covariance explained');
hold off;

%% General Notes for interpreting PLS output       
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%Rotated PLS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*Note*: result directory is listed as resultRotated. This is not standard,
%but simply a change we've made to this specific script. If the name of the
%results directory has not been changed, replacing 'resultRotated' with
%'result' should rectify any problems. 

%1. mean centered brainscores (condition): resultRotated --> boot_result --> orig_usc
	%- columns: latent variables(LVs)
	%- rows: condition
		%- cell: average brainscore for a specific condition in a specific LV
    %- NOTE: for behavioral pls, it is not orig_usc but orig_corr instead!
 
 
%2. mean centered brainscores (subject x condition): resultRotated --> boot_result --> orig_usc
	%- columns: LVs
	%- rows: condition
		%- cell: average brainscore for a specific condition in a specific LV


%3. confidence intervals: resultRotated --> boot_result --> ulusc/llusc
%(ulcorr/llcorr for behavioral pls)
	%- ulusc refers to the upper limit of the confidence interval
	%- llusc refers to the lower limit of the confidence interval
	%- column: LVs
	%- row: condition
		%- cell: confidence interval (distance either above or below usc value) for a specific condition in a specific LV
		%- NOTE: ulusc_adj and llusc_adj are the percentile adjusted versions

%4. proportion of variance explained: result --> s
	%- resultRotated.s(1)^2/sum(resultRotated.s.^2) % number in parentheses reflects the LV of interest
		%- resultRotated.s(1)^2/sum(resultRotated.s.^2) % LV 1
		%- resultRotated.s(2)^2/sum(resultRotated.s.^2) % LV 2
		%- resultRotated.s(3)^2/sum(resultRotated.s.^2) % LV 3		

%5. significance of LV: resultRotated --> perm_result --> sprob
	%- columns: should only be one column
	%- row: LV (first cell = LV 1, and increases downwards)
		%- cells: p-value of specific LV

%6. number of occurences it was replicated using random permutations: result --> perm_result --> sp
	%- columns: should only be one column
	%- row: LV (first cell = LV 1, and increases downwards)
		%- cell: number of times this LV appeared when using random permutations of the data 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Non-rotated PLS%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. brainscores: result --> boot_result --> orig_usc
	%- columns: LVs
	%- rows: condition
		%- cell: average brainscore for a specific condition in a specific LV


%2. confidence intervals: result --> boot_result --> ulusc/llusc
	%- ulusc refers to the upper limit of the confidence interval
	%- llusc refers to the lower limit of the confidence interval
	%- column: LVs
	%- row: condition
		%- cell: confidence interval (distance either above or below usc value) for a specific condition in a specific LV
		%- note: ulusc_adj and llusc_adj are the percentile adjusted versions

%3. proportion of variance explained: result --> s
	%- result.s(1)^2/sum(result.s.^2) % number in parentheses reflects the LV of interest
		%- result.s(1)^2/sum(result.s.^2) % LV 1
		%- result.s(2)^2/sum(result.s.^2) % LV 2
		%- result.s(3)^2/sum(result.s.^2) % LV 3		

%4. significance of LV: result --> perm_result --> sprob
	%- columns: should only be one column
	%- row: LV (first cell = LV 1, and increases downwards)
		%- cells: p-value of specific LV
		%- note: this probability likely does not work for the first LV (if it is one trying to measure the commonality -- e.g., 3 conditions: 1/3 1/3 1/3). It will, however, be fine for other contrasts. If the commonality is important, try to explain using number of occurences during permutations to emphasize it's reliability.

%5. number of occurences it was replicated using random permutations: result --> perm_result --> sp
	%- columns: should only be one column
	%- row: LV (first cell = LV 1, and increases downwards)
		%- cell: number of times this LV appeared when using random permutations of the data 
		%- note: you typically want low values (e.g., 0), to indicate your difference was not obtained spuriously; however, if again you are interested in commonality, you should want high values. This would indicate that the commonality is reliable across multiple permutations (i.e., the commonality is so large, in terms of explained variance, no matter how you permutate the data, it is still apparent). 
