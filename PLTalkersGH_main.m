
%% Set environment



close all
clearvars
clc
restoredefaultpath;
addpath(genpath('/Users/jpohaku/Documents/MATLAB/mTRF-Toolbox-master/'));
addpath  '/Users/jpohaku/Documents/MATLAB/speechenv_toolbox/';
addpath  '/Users/jpohaku/Documents/MATLAB/NoiseTools/'; % http://audition.ens.fr/adc/NoiseTools/
addpath  '/Users/jpohaku/Documents/MATLAB/permutationTest/'
addpath '/Users/jpohaku/Documents/MATLAB/eeglab/eeglab2022.1/';
addpath '/Users/jpohaku/Documents/MATLAB/fieldtrip-master/fieldtrip-master/';
addpath '/Users/jpohaku/Documents/MATLAB/fdr_bh/';
addpath '/Users/jpohaku/Documents/MATLAB/sigstar/';

homefolder='/Volumes/PKUHD/PLTalkers_finalJOCN/Github_PLTalkers/';

dsh=homefolder(end);
cd(homefolder)
addpath(genpath(homefolder))
format long g

setupfolder=[homefolder 'setupfiles' dsh];
outputfolder=[homefolder 'TRF_output' dsh];

% Parallel processing setup - open a parallel pool if not already open
if isempty(gcp('nocreate'))
   % parpool(7);  % Opens a parallel pool with default settings
    parpool('Threads', 14)
end
myrb = redblue(256);







%% Identify TRF output

fileidx=dir(outputfolder);
fileidx = fileidx(~[fileidx.isdir]);


loadname='Decoderoutput_wApV_';
loadname='Decoderoutput_byfreq_wApV_';

loadname='Yaletry_128hz_decfwdmkgpu_condLamAlltrlszp1_ApV_';
loadname='Yaletry_128hz_decfwdmkgpu_condLamAlltrlsbyfreq110zp1_ApV_';

fileidx(~contains({fileidx(:).name},loadname))=[];

%% Load TRF output


    clearvars TRF_study
    % Loop through subjects to load daqta
    for s=1:length({fileidx(:).name})     
        subname=strrep(fileidx(s).name,loadname,'');
        subname=strrep(subname,'.mat','');
        load([fileidx(s).folder dsh fileidx(s).name]);
         if contains(loadname,'freq')
         TRF_study.(subname)=TRF_study_apv;
         else
        TRF_study.(subname)=TRF_all;
        end
        clearvars TRF_master_all subID subname TRF_study_apv
    end
    subnames=fieldnames(TRF_study);
    studynames=fieldnames(TRF_study.(subnames{1}));
    studynames(end)=[];
    stims=fieldnames(TRF_study.(subnames{1}).(studynames{1}));
    stimnames.audio=stims(contains(stims,'env'));
    stimnames.video=stims(~contains(stims,'env'));
    targetstims=stimnames.audio; %
    clearvars stims

    condnames={'AVC','A','V'};

    for s=2:length(subnames)
        TRF_study_apv= TRF_study.(subnames{s});
        for c=1:length(condnames)
                        TRF_study_apv.(subnames{s}).trinity.envmTRF.(condnames{c}).stest=[];
            % 
            % for bp=1:10
            % 
            % end
            % TRF_all.(subnames{s}).trinity.envmTRF.(condnames{c}).strain=[]; 
        end
       % save([ outputfolder  'Decoderoutput_wApV_'  subnames{s} ],'TRF_all');
         save([ outputfolder  'Decoderoutput_byfreq_wApV_'  subnames{s} ],'TRF_study_apv');

        clearvars TRF_study_apv
    end




%% pull decoder results (4 core conditions) 


load([setupfolder 'AVsenBehdata20subs'])
AVsen{:,1}=zscore(AVsen{:,1});
AVsen{:,2}=zscore(AVsen{:,2});

clearvars barColorMap
barColorMap{1}=[0.1 0.5 0.2];	% green
barColorMap{2}=[0.6350 0.0780 0.1840];	% Maroon 
barColorMap{3}=[.84 .68 .22]; % dark yellow
barColorMap{4}=[.25 .55 .79];	% Light blue

condnames={'AVC','AVI','A','V'};
condlabels={'AVc','AVi','Ao','Vo'};
lambdas=TRF_study.(subnames{1}).readme.lambda;
nlambda = length(lambdas);

clearvars tbl tblkeep
current_row=0;
for s=1:length(subnames)
        for c=1:length(condnames)

trlnames=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).stimnames;
            for t=1:length(trlnames)
                current_row=current_row+1;
                tbl{current_row,1}=subnames{s};
                tbl{current_row,2}=(condnames{c});
                tbl{current_row,3}=targetstims{1};
                tbl{current_row,4}=trlnames{t}{1};
                    currtrl=strrep(trlnames{t}{1},'stim','clip');
                    currtrl=strrep(currtrl,'AVI','');
                    currtrl=strrep(currtrl,'AVC','');
                    currtrl=strrep(currtrl,'A','');
                    currtrl=strrep(currtrl,'V','');
                tbl{current_row,5}=currtrl;
                tbl{current_row,6}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).test{1,t}.r;
                tbl{current_row,7}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).test{1,t}.err;
                    currtrl=strrep(trlnames{t}{1},'stim','');
                findtrl=find(contains(behdtatrinityall.([ subnames{s}]).VideoPath,currtrl) & ismember(behdtatrinityall.([ subnames{s}]).conds,condnames{c}));
               if ~isempty(findtrl)
                   tbl{current_row,8}=behdtatrinityall.([ subnames{s}]).respnum{findtrl};
                   tbl{current_row,9}=behdtatrinityall.([ subnames{s}]).key_resp_7_rt(findtrl);
                   tbl{current_row,10}=['TrlCode_' num2str(behdtatrinityall.([subnames{s}]).Trialcode(findtrl))];
               else
                   tbl{current_row,8}=[];
                   tbl{current_row,9}=[];
                   tbl{current_row,10}=[];
               end
%                 tbl{current_row,8}=behdtatrinity.(subnames{s}).resp(ismember(behdtatrinity.(subnames{s}).VideoPath,currtrl));
%                 tbl{current_row,9}=behdtatrinity.(subnames{s}).key_resp_7_keys(ismember(behdtatrinity.(subnames{s}).VideoPath,currtrl));
                 tbl{current_row,11}=AVsen{s,1};
                 tbl{current_row,12}=AVsen{s,2};
                 tbl{current_row,13}=AVsen{s,4};
                 tbl{current_row,14}=std(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).rtest{1,t},[],'all');
            if isequal(condnames{c},'ApV')
                 tbl{current_row,15}=0;%TRF_study.(subnames{s}).(studynames{1}).(targetstims{y}).(condnames{c}).idx;
            else
                 tbl{current_row,15}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).idx;
            end
tbl{current_row,16}=length(trlnames);
tbl{current_row,17}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).lambdause;
%  [~,tbl{current_row,18}]=max(mean(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cv.All.r,1));
% tbl{current_row,19}=lambdas(tbl{current_row,18});
% tbl{current_row,20}=find(lambdas==tbl{current_row,17});


            end
        end
    

end

% once we have vectors, combine into table
tbl = cell2table(tbl,'VariableNames',{'subID','cond','stimtype','clipinfo','clipID','test_r','test_err','resp','resprt'...
    'TrlCode','AVsenC','AVsenI','AVsenD','eegnoise','lambda','numtrls','lambdause'});

% tbl = cell2table(tbl,'VariableNames',{'subID','cond','stimtype','clipinfo','clipID','test_r','test_err','resp','resprt'...
%     'TrlCode','AVsenC','AVsenI','AVsenD','eegnoise','lambda','numtrls','lambdause','Alllambdaidx','Alllambda','lambdauseidx'});


% Write the table to a CSV file
% writetable(tbl, 'PLTalkers_decoder_dtatable.csv');

% tbl.test_r(:)=abs(tbl.test_r(:));

tblkeep = PLTalkers_makedectbl(tbl,condnames,subnames,targetstims);


clearvars submeans suberrs suballlambs sublambdaused subcvs sublambdaidx subtrlcnt
for s=1: length(subnames)
    for c=1:length(condnames)
        submeans(s,c)=mean(tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,condnames{c})));
        suberrs(s,c)=std(tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,condnames{c})));
        sublambdaidx(s,c)=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).idx.(condnames{c});
        subcvs(s,:,:) =   TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cv.(condnames{c}).r;
        subtrlcnt(s,c)= length(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).test);
        sublambdaused(s,c)=find(lambdas==TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).lambdause);

    end
    %suballlambs(s,1)=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).idx.All;
end






% Plot mean decoder score by subject 
close all
for s=1: length(subnames)
    subplot(5,4,s)
    for c=1:length(condnames)
        h=bar(c,submeans(s,c));
        h.FaceColor=barColorMap{c};
        hold on
    end
    ylim([-0.1 0.3])
    xlim([0 6])
    title(subnames{s})
    grid on
end
set(gcf, 'Position', [494,60,510,801]);


% Plot lambda idx per condition - barchart
close all
for s=1: length(subnames)
    subplot(5,4,s)
    for c=1:length(condnames)
        h=bar(c,sublambdaidx(s,c));
        h.FaceColor=barColorMap{c};
        hold on
    end
     % h=bar(length(condnames)+1,suballlambs(s,1));
     % h.FaceColor=[0 0 0];
     %   h=bar(length(condnames)+2,find(lambdas==sublambdaused(s,1)));
     % h.FaceColor=[0 0 0];
    ylim([0 nlambda])
    xlim([0 length(condnames)+1])
    title(subnames{s})
    set(gca,'xtick',1:length(condnames),'xticklabel',{sublambdaidx(s,:)},...
        'ytick',1:nlambda,'yticklabel',{1:nlambda})
    grid on
end
set(gcf, 'Position', [494,60,510,801]);

% Plot lambda used per condition - barchart
close all

for s=1: length(subnames)
    subplot(5,4,s)
    for c=1:length(condnames)
        h=bar(c,sublambdaused(s,c));
        h.FaceColor=barColorMap{c};
        hold on
    end
     % h=bar(length(condnames)+1,suballlambs(s,1));
     % h.FaceColor=[0 0 0];
     %   h=bar(length(condnames)+2,find(lambdas==sublambdaused(s,1)));
     % h.FaceColor=[0 0 0];
    ylim([0 nlambda])
    xlim([0 length(condnames)+1])
    title(subnames{s})
    set(gca,'xtick',1:length(condnames),'xticklabel',{sublambdaused(s,:)},...
        'ytick',1:nlambda,'yticklabel',{1:nlambda})
    grid on
end
set(gcf, 'Position', [494,60,510,801]);


% Plot trials per condition - barchart
close all
for s=1: length(subnames)
    subplot(5,4,s)
    for c=1:length(condnames)
        h=bar(c,subtrlcnt(s,c));
        h.FaceColor=barColorMap{c};
        hold on
    end
     % h=bar(length(condnames)+1,suballlambs(s,1));
     % h.FaceColor=[0 0 0];
     %   h=bar(length(condnames)+2,find(lambdas==sublambdaused(s,1)));
     % h.FaceColor=[0 0 0];
    ylim([0 35])
    xlim([0 length(condnames)+1])
    title(subnames{s})
    set(gca,'xtick',1:length(condnames),'xticklabel',{subtrlcnt(s,:)},...
        'ytick',1:35,'yticklabel',{1:35})
    grid on
end
set(gcf, 'Position', [494,60,510,801]);




% Plot lambdas used scatter
% suballlambs
% sublambdaused
close all
for s=1: length(subnames)
    h=scatter(s,log(sublambdaused(s,1)),'filled','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.8);
    hold on
end
ylim(log([lambdas(13-5) lambdas(13+5)]))
xlim([0 length(subnames)+1])
title('Regularization parameter (λ) identified')
set(gca,'xtick',1:length(subnames),'xticklabel',strrep(subnames,'TCD_sub0','S '),...
    'ytick',[linspace(log(lambdas(13-5)),log(lambdas(13+5)),5)],'yticklabel',{'woo'})
grid on
set(gcf, 'Position', [220,383,624,116]);

%'$\lambda$', 'Interpreter', 'latex'


% Plot mean decoder score across subjects
close all
    for c=1:length(condnames)
        h=bar(c,mean(submeans(:,c)));
        h.FaceColor=barColorMap{c};
        hold on
    end
    ylim([-0.1 0.25])
    xlim([0 length(condnames)+1])
    title('Group Means')
    grid on
set(gcf, 'Position', [476,360,245,329]);



mybars = PLTalkers_barplotwdots(stimnames,condnames,condlabels,barColorMap,tblkeep,'yes');


%% Stats on naturalness ratings - ordinal regression

% Multinomial regression
% predictor = Conditions AVC AVI Ao Vo ; response = Likert 1 2 or 3
% mnrmdl = fitmnr(predictor, response); % Convert categorical Y to numeric

% Fit an ordinal multinomial regression model using Acceleration, 
% Displacement, Horsepower, and Weight as predictor variables 
% and Mileage as the response variable. 
% fitmnr ignores the rows of X containing the undefined entries in the Mileage column.
% This type of model is sometimes called a proportional odds model
% mileage = discretize(MPG,[9,19,29,39,48],"categorical");
% X = table(Acceleration,Displacement,Horsepower,Weight,mileage,VariableNames=["Acceleration","Displacement","Horsepower","Weight","Mileage"])
% MnrModel = fitmnr(X,"Mileage",ModelType="ordinal")

% clearvars mnr_dv mnr_iv
% for ii=1:size(tblkeep.core4conds,1)
% 
%       % mnr_iv{ii,1}=['R_' num2str(tblkeep.core4conds.resp(ii))];
%        mnr_iv(ii,1)=tblkeep.core4conds.resp(ii);
%        mnr_dv{ii,1}=tblkeep.core4conds.cond{ii};
% end

clearvars mnrtbl mnrtblnum respcounts MnrModelO MnrModelN
mnrtbl=vertcat(tblkeep.AVC,tblkeep.AVI); % mnrtbl=tblkeep.core4conds;
mnrtbl = removevars(mnrtbl,{'subID','stimtype','clipinfo','clipID','test_r','test_err',...
    'resprt','TrlCode','AVsenC','AVsenI','AVsenD','eegnoise','lambda','numtrls','lambdause'});

% {'subID','cond','stimtype','clipinfo','clipID','test_r','test_err','resp','resprt'...
%     'TrlCode','AVsenC','AVsenI','AVsenD','eegnoise','lambda','numtrls'}

% mnrtbl = removevars(mnrtbl,{'subID','cond','stimtype','clipinfo','clipID',...
%     'TrlCode','AVsenC','AVsenI'});
mnrtbl.cond = categorical(mnrtbl.cond);
mnrtbl.cond=reordercats(mnrtbl.cond,{'AVC','AVI'});
mnrtblnum=mnrtbl;
disp(mean(mnrtbl.resp(1:583)));
disp(mean(mnrtbl.resp(584:end)));
mnrtbl.resp = categorical(mnrtbl.resp);  
unique(mnrtbl.cond)
unique(mnrtbl.resp)

respcounts(1,1)=sum(mnrtblnum.resp(ismember(mnrtblnum.cond,'AVC'))==1);
respcounts(1,2)=sum(mnrtblnum.resp(ismember(mnrtblnum.cond,'AVC'))==2);
respcounts(1,3)=sum(mnrtblnum.resp(ismember(mnrtblnum.cond,'AVC'))==3);
respcounts(2,1)=sum(mnrtblnum.resp(ismember(mnrtblnum.cond,'AVI'))==1);
respcounts(2,2)=sum(mnrtblnum.resp(ismember(mnrtblnum.cond,'AVI'))==2);
respcounts(2,3)=sum(mnrtblnum.resp(ismember(mnrtblnum.cond,'AVI'))==3);

clc;disp(respcounts);

MnrModelO = fitmnr(mnrtbl,"resp",ModelType="ordinal"); % MnrModel = fitmnr(tbl,"Species ~ Petal_Length*Sepal_Length + Petal_Width + Sepal_Width")
clc;disp(MnrModelO);

% Example estimated parameters from fitmnr
alpha1 = -2.97015883419932;  % First intercept
alpha2 = -1.14954800107329;   % Second intercept
beta = 3.42907434417434;     % Condition effect (AVI vs. AVC)

% Define predictor values
X_AVC = 0;  % Reference category (AVC)
X_AVI = 1;  % Condition AVI

% Compute cumulative probabilities for AVC (reference)
P_leq1_AVC = exp(alpha1 + beta*X_AVC) / (1 + exp(alpha1 + beta*X_AVC));
P_leq2_AVC = exp(alpha2 + beta*X_AVC) / (1 + exp(alpha2 + beta*X_AVC));

% Compute cumulative probabilities for AVI
P_leq1_AVI = exp(alpha1 + beta*X_AVI) / (1 + exp(alpha1 + beta*X_AVI));
P_leq2_AVI = exp(alpha2 + beta*X_AVI) / (1 + exp(alpha2 + beta*X_AVI));

% Convert to category probabilities
P_Y1_AVC = P_leq1_AVC;  % P(Y = 1) = P(Y ≤ 1)
P_Y2_AVC = P_leq2_AVC - P_leq1_AVC; % P(Y = 2) = P(Y ≤ 2) - P(Y ≤ 1)
P_Y3_AVC = 1 - P_leq2_AVC; % P(Y = 3) = 1 - P(Y ≤ 2)

P_Y1_AVI = P_leq1_AVI;
P_Y2_AVI = P_leq2_AVI - P_leq1_AVI;
P_Y3_AVI = 1 - P_leq2_AVI;

% Display results
fprintf('Probabilities for AVC:\n');
fprintf('P(Y=1) = %.3f\n', P_Y1_AVC);
fprintf('P(Y=2) = %.3f\n', P_Y2_AVC);
fprintf('P(Y=3) = %.3f\n\n', P_Y3_AVC);

fprintf('Probabilities for AVI:\n');
fprintf('P(Y=1) = %.3f\n', P_Y1_AVI);
fprintf('P(Y=2) = %.3f\n', P_Y2_AVI);
fprintf('P(Y=3) = %.3f\n', P_Y3_AVI);

MnrModel_Null = fitmnr(mnrtbl(:,2), "resp", ModelType="ordinal");
clc;disp(MnrModel_Null);

% Extract log-likelihoods
LL_full = MnrModelO.LogLikelihood;
LL_null = MnrModel_Null.LogLikelihood;

LRT_stat = -2 * (LL_null - LL_full);
p_value = 1 - chi2cdf(LRT_stat, 1);
fprintf('Likelihood Ratio Test:\n');
fprintf('Log-Likelihood (Full Model) = %.3f\n', LL_full);
fprintf('Log-Likelihood (Intercept-Only Model) = %.3f\n', LL_null);
fprintf('Chi-square statistic = %.3f\n', LRT_stat);
% fprintf('Degrees of freedom = %d\n', df);
fprintf('p-value = %.5f\n', p_value);
clc;disp(MnrModelO);


% Define table title
tableTitle = 'Table 1. Predicted Naturalness Ratings';
% Column headers
columns = {'', 'AVc', 'AVi'}; % First column is empty for row labels
% Updated row labels with full response descriptions
row_labels = {
    'P(Response = "Unnatural")',...
    'P(Response = "Moderately Natural")',... 
    'P(Response = "Very Natural")'};
% Example predicted probabilities (modify as needed)
data = [
    0.049, 0.613;  % P(Response = "Unnatural")
    0.192, 0.294;  % P(Response = "Moderately Natural")
    0.759, 0.093  % P(Response = "Very Natural")
];
% Print table title
fprintf('%s\n', tableTitle);
fprintf(repmat('-', 1, 60)); % Horizontal line
fprintf('\n');
% Print column headers
fprintf('%-35s %-10s %-10s\n', columns{:});
fprintf(repmat('-', 1, 60)); % Horizontal line
fprintf('\n');
% Print each row of the table
for i = 1:size(data, 1)
    fprintf('%-35s %-10.2f %-10.2f\n', row_labels{i}, data(i, :));
end
% Print notes (optional)
fprintf('\nNote: Predicted probabilities are based on an ordinal logistic regression model.\n');








%% Make additive model weights (decoders)
tic
% % Specify additve DECODER model setup
pkucfg=[];
pkucfg.modeltype='backward'; 
pkucfg.runpredtest='yes';
pkucfg.Dir = -1; % direction of causality
pkucfg.tmin = -100; % minimum time lag (ms)
pkucfg.tmax = 500; % maximum time lag (ms)
pkucfg.fs=128;
pkucfg.lambda = 10.^(-24:2:24); % regularization parameters (ridge values?)
pkucfg.zeropad=0;
pkucfg.unimodels={'A','V'};
pkucfg.newmdlname={'ApV'};
pkucfg.newmdlcompare={'AVC'};
pkucfg.conds2pred={'ApV','AVC','AVI','A','V'};
pkucfg.cond2pullstim={'AVC','AVC','AVI','A','V'};
TRF_study = PLTalkers_makeadditive_sublambda(pkucfg,TRF_study);clc;
toc

    for s=1:length(subnames)     
        %load([ outputfolder 'TRF_study_apv_' subnames{s} ]);
        TRF_all= TRF_study.(subnames{s});
         % load([ outputfolder 'TRF_apv_' loadname subnames{s} ]);
         save([ outputfolder loadname 'ApV_'  subnames{s} ],'TRF_all');

        clearvars TRF_all 
    end



load([setupfolder 'AVsenBehdata20subs'])
AVsen{:,1}=zscore(AVsen{:,1});
AVsen{:,2}=zscore(AVsen{:,2});

clearvars barColorMap
barColorMap{1}=[0.1 0.5 0.2];	% green
barColorMap{2}=[0.6350 0.0780 0.1840];	% Maroon 
barColorMap{3}=[0.8 0.8 0.8];	% grey
barColorMap{4}=[.84 .68 .22]; % dark yellow
barColorMap{5}=[0 0.451 0.7412];	% Light blue

condnames={'AVC','AVI','ApV','A','V'};
condlabels={'AVc','AVi','[Ao + Vo]','Ao','Vo'};
lambdas=TRF_study.(subnames{1}).readme.lambda;
nlambda = length(lambdas);

clearvars tbl tblkeep
current_row=0;
for s=1:length(subnames)
        for c=1:length(condnames)
             if isequal(condnames{c},'ApV')
            trlnames=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).AVC.stimnames;
            else
            trlnames=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).stimnames;
            end
            for t=1:length(trlnames)
                clearvars lambdause
                current_row=current_row+1;
                tbl{current_row,1}=subnames{s};
                tbl{current_row,2}=(condnames{c});
                tbl{current_row,3}=targetstims{1};
                tbl{current_row,4}=trlnames{t}{1};
                    currtrl=strrep(trlnames{t}{1},'stim','clip');
                    currtrl=strrep(currtrl,'AVI','');
                    currtrl=strrep(currtrl,'AVC','');
                    currtrl=strrep(currtrl,'A','');
                    currtrl=strrep(currtrl,'V','');
                tbl{current_row,5}=currtrl;
                tbl{current_row,6}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).test{1,t}.r;
                tbl{current_row,7}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).test{1,t}.err;
                    currtrl=strrep(trlnames{t}{1},'stim','');
                findtrl=find(contains(behdtatrinityall.([ subnames{s}]).VideoPath,currtrl) & ismember(behdtatrinityall.([ subnames{s}]).conds,condnames{c}));
               if ~isempty(findtrl)
                   tbl{current_row,8}=behdtatrinityall.([ subnames{s}]).respnum{findtrl};
                   tbl{current_row,9}=behdtatrinityall.([ subnames{s}]).key_resp_7_rt(findtrl);
                   tbl{current_row,10}=['TrlCode_' num2str(behdtatrinityall.([subnames{s}]).Trialcode(findtrl))];
               else
                   tbl{current_row,8}=[];
                   tbl{current_row,9}=[];
                   tbl{current_row,10}=[];
               end
%                 tbl{current_row,8}=behdtatrinity.(subnames{s}).resp(ismember(behdtatrinity.(subnames{s}).VideoPath,currtrl));
%                 tbl{current_row,9}=behdtatrinity.(subnames{s}).key_resp_7_keys(ismember(behdtatrinity.(subnames{s}).VideoPath,currtrl));
                 tbl{current_row,11}=AVsen{s,1};
                 tbl{current_row,12}=AVsen{s,2};
                 tbl{current_row,13}=AVsen{s,4};
                 tbl{current_row,14}=std(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).rtest{1,t},[],'all');
            if isequal(condnames{c},'ApV')
                 tbl{current_row,15}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).AVC.idx.AVC;
                 lambdause=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).AVC.lambdause;
            else
                 tbl{current_row,15}=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).idx.(condnames{c});
                                  lambdause=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).lambdause;

            end
tbl{current_row,16}=length(trlnames);
tbl{current_row,17}=lambdause;
%             if isequal(condnames{c},'ApV')
%                  [~,tbl{current_row,18}]=max(mean(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).AVC.cv.All.r,1));
% 
%             else
%              [~,tbl{current_row,18}]=max(mean(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cv.All.r,1));
% 
%             end
% tbl{current_row,19}=lambdas(tbl{current_row,18});
% tbl{current_row,20}=find(lambdas==tbl{current_row,17});


            end
        end
    

end

% once we have vectors, combine into table
tbl = cell2table(tbl,'VariableNames',{'subID','cond','stimtype','clipinfo','clipID','test_r','test_err','resp','resprt'...
    'TrlCode','AVsenC','AVsenI','AVsenD','eegnoise','lambda','numtrls','lambdause'});




tblkeep = PLTalkers_makedectbl(tbl,condnames,subnames,targetstims);



clearvars barColorMap
barColorMap{1}=[0.5 0.5 0.5];	% grey
barColorMap{2}=[0.1 0.5 0.2];	% green
barColorMap{3}=[0.6350 0.0780 0.1840];	% Maroon 
barColorMap{4}=[.84 .68 .22]; % dark yellow
barColorMap{5}=[0 0.451 0.7412];	% Light blue

condnames={'ApV','AVC','AVI','A','V'};
condlabels={'[Ao + Vo]','AVc','AVi','Ao','Vo'};

mybars = PLTalkers_barplotwdots(stimnames,condnames,condlabels,barColorMap,tblkeep,'yes');

lambdas=10.^(-24:2:24);
nlambda = length(lambdas);
clearvars submeans suberrs suballlambs sublambdaused subcvs sublambdaidx subtrlcnt sublambdausedApV
for s=1: length(subnames)
    for c=1:length(condnames)
        submeans(s,c)=mean(tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,condnames{c})));
        suberrs(s,c)=std(tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,condnames{c})));
        subtrlcnt(s,c)= length(TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).test);
        if isequal(condnames{c},'ApV')
          % sublambdausedApV(s,1)=find(lambdas==TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cvidxA);
           %sublambdausedApV(s,2)=find(lambdas==TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cvidxV);
            sublambdausedApV(s,1)=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cvidxA;
           sublambdausedApV(s,2)=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cvidxV;
sublambdaused(s,c)=0;
sublambdaidx(s,c)=0;
        else
         subcvs(s,c,:,:) =   TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).cv.(condnames{c}).r;
           sublambdaused(s,c)=find(lambdas==TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).lambdause);
        sublambdaidx(s,c)=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).idx.(condnames{c});

        end

    end
    %suballlambs(s,1)=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).idx.All;
end








%% parametric assumption testing


close all
% Fit the linear mixed-effects model
mdl = fitlme(tbl, 'test_r ~ 1 + cond + (1|subID) + (1|clipID)', 'FitMethod','REML');
% Extract residuals
myresiduals = residuals(mdl);
% Histogram
subplot(2,1,1)
histogram(myresiduals);
title('Histogram of Residuals');
% QQ-plot
subplot(2,1,2)
qqplot(myresiduals);
title('QQ-plot of Residuals');
% Shapiro-Wilk test (alternative: lillietest, Kolmogorov-Smirnov)
[h, p] = swtest(myresiduals);
disp(['Shapiro-Wilk p-value: ', num2str(p)]);
% If p < 0.05, residuals significantly deviate from normality.



% Extract predictor matrix (excluding intercept)
X = dummyvar(categorical(tbl.cond)); % Convert to dummy variables
% Run Durbin-Watson test
[p_dw, stat_dw] = dwtest(myresiduals, X);
disp(['Durbin-Watson p-value: ', num2str(p_dw)]);
disp(['Durbin-Watson statistic: ', num2str(stat_dw)]);
% If p < 0.05, residuals are autocorrelated, violating the assumption.
% If Durbin-Watson ≈ 2, no autocorrelation.
% < 2 → Positive autocorrelation.
% > 2 → Negative autocorrelation.





% Example of Shapiro-Wilk test for normality (for single continuous variable)
clearvars shap_wilk
clc;
for c=1:length(condnames)

    [H, pValue, SWstatistic]=swtest(tbl.test_r, 0.05);
disp(H);
disp(pValue);
% [shap_wilk(c,1), shap_wilk(c,2)] = swtest(tblkeep.(condnames{c}).test_r, 0.05);

end

[H, pValue, SWstatistic]=swtest(tblkeep.AVC.test_r, 0.05);

%% permutation-based power analysis (simulation) 

% Define parameters
numSubjects = 20;
numTrialsPerCond = 32; % 32 trials per condition
numPermutations = 10000;
permcompare={'ApV','AVC'};
% data = randn(numSubjects, numTrialsPerCond * 2); % 20 subjects, 64 trials (32 per condition)
labels = [ones(1, numTrialsPerCond), 2*ones(1, numTrialsPerCond)]; % Condition labels (1 and 2)


for pp=1:size(permcompare,1)
    % Simulated data: Replace this with your real dataset
    clearvars data group1 group2
    for s = 1:numSubjects
        data.AVC(s,:)=tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'AVC'));
        data.ApV(s,:)=tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'ApV'));
        data.A(s,:)=tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'A'));
        data.AVI(s,:)=tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'AVI'));
    end
    data=  horzcat(data.(permcompare{pp,1}),data.(permcompare{pp,2}));

N = 100;   % Number of values
mu1 = 0.1;  % Desired mean
mu2 = 0.2;  % Desired mean

sim1 = 2 * rand(20,32) - 1;  % Generate uniform values between -1 and 1
sim1 = sim1 - mean(sim1) + mu1;   % Adjust mean to the desired value
sim2 = 2 * rand(20,32) - 1;  % Generate uniform values between -1 and 1
sim2 = sim2 - mean(sim2) + mu2;   % Adjust mean to the desired value
    simdata=  horzcat(sim1,sim2);

        % simdata=randn(20,64);
        % simdata = 2*rand(20,64) - 1;
        % simdata(:,33:end)=simdata(:,33:end) + 0.05;

data=simdata;

    % Compute real mean difference across conditions for each subject
    group1 = mean(data(:, labels == 1), 2);group2= mean(data(:, labels == 2), 2);

    realDiffs = mean(data(:, labels == 1), 2) - mean(data(:, labels == 2), 2);
    realEffect = mean(realDiffs); % Group-level effect size

    % Initialize null distribution
    nullEffects = zeros(numPermutations, 1);

    % Null permutation loop
    tic
    parfor p = 1:numPermutations
        permutedData = data; % Copy data for shuffling
        permutedDiffs=zeros(20,1);
        for s = 1:numSubjects
            shuffledLabels = labels(randperm(length(labels))); % Shuffle condition labels
            permutedDiffs(s,1) = mean(permutedData(s, shuffledLabels == 1)) - mean(permutedData(s, shuffledLabels == 2));
        end
        nullEffects(p) = mean(permutedDiffs); % Store null effect
    end
    toc

    % Compute p-value (two-tailed test)
    p_value = mean(abs(nullEffects) >= abs(realEffect));
    
    % Estimate post-hoc power (proportion of permutations with a more extreme test statistic)
    postHocPower = mean(abs(nullEffects) >= abs(realEffect));
    
    
    % Plot null distribution and real effect
    histogram(nullEffects, 50, 'FaceColor', 'k', 'EdgeColor', 'k', 'Normalization', 'probability');
    hold on;
    yLimits = ylim;
    plot([realEffect realEffect], yLimits, 'r', 'LineWidth', 2);
    xlabel('Mean Difference (Null)');
    ylabel('Probability');
    title(['Null Distribution (p = ' num2str(p_value) ') ' permcompare{pp,1} ':' permcompare{pp,2}]);
    legend({'Null Distribution', 'Real Effect'});

end

% Power Analysis for Two-Group Comparison
m1 = mean(group1);m2 = mean(group2);

% Calculate sample sizes
n1 = length(group1);n2 = length(group2);
% Calculate sample standard deviations
s1 = std(group1);s2 = std(group2);
% Calculate pooled standard deviation
sp = sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2));
% Define input parameters
alpha = 0.05;    % Significance level (α)
power = 0.80;    % Desired power (1 - β)
sigma = sp;      % Estimated standard deviation (σ)
cohensD=realEffect/sigma;
noncentrality=cohensD * sqrt((n1*n2)/(n1+n2));
delta = realEffect;       % Desired minimum detectable difference (Δ)
% Calculate Z-scores for alpha/2 and beta
Z_alpha = norminv(1 - alpha / 2);  % Z-score for significance level (two-tailed)
Z_beta = norminv(power);           % Z-score for desired power
% Calculate sample size per group
n = (2 * (Z_alpha + Z_beta)^2 * sigma^2) / delta^2;
% Round up to the nearest whole number
n = ceil(n);
% Display the result
fprintf('The required sample size per group is: %d\n', n);

disp(m1);
disp(m2);
disp(s1);
disp(s2);








% Subject-level analysis
close all
for s=1: length(subnames)
    subplot(5,4,s)
    group1 = tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'AVC'));
    group2 = tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'ApV'));
    group3 = tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'A'));
    histogram(group1,32,'FaceColor',barColorMap{1})
    xline(mean(group1),'Color',barColorMap{1},'LineWidth',2,'LineStyle','-')
    xline(mean(group1)+(2*std(group1)),'Color',barColorMap{1},'LineWidth',1,'LineStyle','--')
    xline(mean(group1)-(2*std(group1)),'Color',barColorMap{1},'LineWidth',1,'LineStyle','--')
    hold on
    histogram(group2,32,'FaceColor',barColorMap{2})
    xline(mean(group2),'Color',barColorMap{2},'LineWidth',2,'LineStyle','-')
    xline(mean(group2)+(2*std(group2)),'Color',barColorMap{2},'LineWidth',1,'LineStyle','--')
    xline(mean(group2)-(2*std(group2)),'Color',barColorMap{2},'LineWidth',1,'LineStyle','--')
    histogram(group3,32,'FaceColor',barColorMap{3})
    xline(mean(group3),'Color',barColorMap{3},'LineWidth',2,'LineStyle','-')
    xline(mean(group3)+(2*std(group3)),'Color',barColorMap{3},'LineWidth',1,'LineStyle','--')
    xline(mean(group3)-(2*std(group3)),'Color',barColorMap{3},'LineWidth',1,'LineStyle','--')
    ylim([0 5])
    xlim([-0.5 1])
    grid on
    % Compute Cohen's D
    d = PLT_computeCohensD(group1, group2);
    disp(['Mean Cohen''s D across time points: ', num2str(mean(d))]);
    % Estimate Required Sample Size
    power = 0.8;  % Desired power (80%)
    alpha = 0.05; % Significance level (5%)
    sample_size = PLT_estimateSampleSize(d, power, alpha); % (d, power, alpha)
    disp(['Required sample size per group: ', num2str(ceil(sample_size))]);
    title([' n=' num2str(round(sample_size)) ])
end
sgtitle('AVc vs ApV')
set(gcf, 'Position', [494,60,510,801]);




% Group-level analysis 
close all 
for s=1: length(subnames)
    if s==1
        group1 = tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'AVC'));
        group2 = tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'ApV'));
        group3 = tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'A'));
    else
        group1 = vertcat(group1,tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'AVC')));
        group2 = vertcat(group2,tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'ApV')));
        group3 = vertcat(group3,tblkeep.(subnames{s}).test_r(ismember(tblkeep.(subnames{s}).cond,'A')));
    end
end

histogram(group1,32,'FaceColor',barColorMap{1})
xline(mean(group1),'Color',barColorMap{1},'LineWidth',2,'LineStyle','-')
xline(mean(group1)+(2*std(group1)),'Color',barColorMap{1},'LineWidth',1,'LineStyle','--')
xline(mean(group1)-(2*std(group1)),'Color',barColorMap{1},'LineWidth',1,'LineStyle','--')
hold on
histogram(group2,32,'FaceColor',barColorMap{2})
xline(mean(group2),'Color',barColorMap{2},'LineWidth',2,'LineStyle','-')
xline(mean(group2)+(2*std(group2)),'Color',barColorMap{2},'LineWidth',1,'LineStyle','--')
xline(mean(group2)-(2*std(group2)),'Color',barColorMap{2},'LineWidth',1,'LineStyle','--')
histogram(group3,32,'FaceColor',barColorMap{3})
xline(mean(group3),'Color',barColorMap{3},'LineWidth',2,'LineStyle','-')
xline(mean(group3)+(2*std(group3)),'Color',barColorMap{3},'LineWidth',1,'LineStyle','--')
xline(mean(group3)-(2*std(group3)),'Color',barColorMap{3},'LineWidth',1,'LineStyle','--')
%ylim([0 5])
xlim([-0.5 1])
grid on
% Compute Cohen's D
d = PLT_computeCohensD(group1, group2);
disp(['Mean Cohen''s D across time points: ', num2str(mean(d))]);
% Estimate Required Sample Size
power = 0.8;  % Desired power (80%)
alpha = 0.05; % Significance level (5%)
sample_size = PLT_estimateSampleSize(d, power, alpha); % (d, power, alpha)
disp(['Required sample size per group: ', num2str(ceil(sample_size))]);
title(['AVc vs ApV n=' num2str(round(sample_size)) ])
set(gcf, 'Position', [494,60,510,801]);




% Fieldtrip Permutation Test for Statistical Power Estimation
ft_defaults;

% Simulate Example Data for Two Conditions
nSubjects = 20; % Change to test power at different sample sizes
nPermutations = 1000; % Number of Monte Carlo iterations

% Simulate EEG power values (normally distributed)
mu1 = 5; sigma1 = 1; % Condition 1 (Mean, SD)
mu2 = 6; sigma2 = 1; % Condition 2 (Effect: 1 unit)
condition1 = mu1 + sigma1 * randn(nSubjects, 1);
condition2 = mu2 + sigma2 * randn(nSubjects, 1);

% Compute Observed Effect Size (Cohen’s d)
effect_size = (mean(condition2) - mean(condition1)) / std([condition1; condition2]);
fprintf('Effect Size (Cohen’s d): %.3f\n', effect_size);

real_t = abs(ttest2(condition1, condition2)); % Compute observed t-statistic
perm_t = zeros(nPermutations, 1);

for i = 1:nPermutations
    % Shuffle condition labels
    all_data = [condition1; condition2];
    shuffled_labels = all_data(randperm(length(all_data)));
    
    % Split shuffled data into two groups
    shuffled1 = shuffled_labels(1:nSubjects);
    shuffled2 = shuffled_labels(nSubjects+1:end);
    
    % Compute t-statistic for shuffled data
    perm_t(i) = abs(ttest2(shuffled1, shuffled2));
end

% Compute Statistical Power
alpha = 0.05; % Significance level
p_value = mean(perm_t >= real_t); % Compute p-value
power = mean(p_value < alpha); % Power estimation (1 - β)

fprintf('Estimated Power: %.3f (for N = %d)\n', power, nSubjects);


%% By-frequency analysis (model re-train)
 


pkucfg=[];
pkucfg.modeltype='backward'; 
pkucfg.runpredtest='yes';
pkucfg.Dir = -1; % direction of causality
pkucfg.tmin = -100; % minimum time lag (ms)
pkucfg.tmax = 500; % maximum time lag (ms)
pkucfg.fs=128;
pkucfg.lambda = 10.^(-24:2:24); % regularization parameters (ridge values?)
pkucfg.zeropad=0;
pkucfg.unimodels={'A','V'};
pkucfg.newmdlname={'ApV'};
pkucfg.newmdlcompare={'AVC'};
pkucfg.conds2pred={'ApV','AVC','AVI','A','V'};
pkucfg.cond2pullstim={'AVC','AVC','AVI','A','V'};
pkucfg.freqs=1:10;

TRF_study = PLTalkers_makeadditive_sublambda_byfreq(pkucfg,TRF_study);clc;


for s=1:length(subnames) 
    tic
    clearvars TRF_study_apv
    TRF_study_apv=TRF_study.(subnames{s});
    save([ outputfolder loadname 'ApV_' subnames{s} ],'TRF_study_apv','-v7.3');
toc
end
    clearvars TRF_study_apv




clearvars barColorMap
barColorMap{1}=[0.1 0.5 0.2];	% green
barColorMap{2}=[0.6350 0.0780 0.1840];	% Maroon 
barColorMap{3}=[0.5 0.5 0.5];	% grey
barColorMap{4}=[.84 .68 .22]; % dark yellow
barColorMap{5}=[.25 .55 .79];	% Light blue

condnames={'ApV','AVC'};%,'A','V'};
condlabels={'[Ao + Vo]','AVc'};%,'Ao','Vo'};
clearvars barColorMap
barColorMap{1}=[0.5 0.5 0.5];	% grey
barColorMap{2}=[0.1 0.5 0.2];	% green 

condnames={'A','AVC'};%,'A','V'};
condlabels={'Ao','AVc'};%,'Ao','Vo'};
myfreqs=1:10;
clearvars barColorMap
barColorMap{1}=[.84 .68 .22];	% y
barColorMap{2}=[0.1 0.5 0.2];	% green 

clearvars byfreqtbl byfreqmeans byfreqerr byfreqsubmeans
counter=0;
for bp=myfreqs
    for s=1:length(subnames)
        for c=1:length(condnames)
            if isequal(condnames{c},'ApV')
                trlnames=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).AVC.stimnames;
            else
                trlnames=TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).stimnames;
            end
            for t=1:length(trlnames)
                counter=counter+1;
                byfreqtbl{counter,1}=    TRF_study.(subnames{s}).(studynames{1}).(targetstims{1}).(condnames{c}).(['Fq_' num2str(bp)]).test{1,t}.r;
                byfreqtbl{counter,2}= subnames{s};
                byfreqtbl{counter,3}= condnames{c};
                byfreqtbl{counter,4}= ['Fq_' num2str(bp)];
                byfreqtbl{counter,5}= trlnames{t};
            end
        end
    end
end
for bp=myfreqs
    for c=1:length(condnames)
        clearvars freqmeans
        freqmeans=  [byfreqtbl{ismember(byfreqtbl(:,4),['Fq_' num2str(bp)]) & ismember(byfreqtbl(:,3),condnames{c}),1}];
        freqmeans(isnan(freqmeans))=[];
        byfreqmeans(bp,c)=mean(freqmeans);
        byfreqerr(bp,c)=std(freqmeans)/sqrt(length(freqmeans));
    end
end
for bp=myfreqs
    for s=1:20
        for c=1:length(condnames)
            clearvars freqmeans
            freqmeans=  [byfreqtbl{ismember(byfreqtbl(:,2),subnames{s}) & ismember(byfreqtbl(:,4),['Fq_' num2str(bp)]) & ismember(byfreqtbl(:,3),condnames{c}),1}];
            freqmeans(isnan(freqmeans))=[];

            byfreqsubmeans(bp,c,s,:)=freqmeans;
            % byfreqerr(bp,c)=std(freqmeans)/sqrt(length(freqmeans));
        end
    end
end

% Define parameters
numSubjects = 20;
numTrialsPerCond = 32; % 32 trials per condition
numPermutations = 1000;
labels = [ones(1, numTrialsPerCond), 2*ones(1, numTrialsPerCond)]; % Condition labels (1 and 2)
clearvars byfreqstats
for bp=myfreqs

    data=  horzcat(squeeze(byfreqsubmeans(bp,1,:,:)),squeeze(byfreqsubmeans(bp,2,:,:)));
        % Compute real mean difference across conditions for each subject
    group1 = mean(data(:, labels == 1), 2);group2= mean(data(:, labels == 2), 2);

    realDiffs = mean(data(:, labels == 1), 2) - mean(data(:, labels == 2), 2);
    realEffect = mean(realDiffs); % Group-level effect size

    % Initialize null distribution
    nullEffects = zeros(numPermutations, 1);

    % Null permutation loop
    tic
    parfor p = 1:numPermutations
        permutedData = data; % Copy data for shuffling
        permutedDiffs=zeros(20,1);
        for s = 1:numSubjects
            shuffledLabels = labels(randperm(length(labels))); % Shuffle condition labels
            permutedDiffs(s,1) = mean(permutedData(s, shuffledLabels == 1)) - mean(permutedData(s, shuffledLabels == 2));
        end
        nullEffects(p) = mean(permutedDiffs); % Store null effect
    end
    toc

    % Compute p-value (two-tailed test)
    p_value = mean(abs(nullEffects) >= abs(realEffect));

    byfreqstats(bp,1)=mean(group1);
    byfreqstats(bp,2)=mean(group2);
    byfreqstats(bp,3)=std(group1)/sqrt(20);
    byfreqstats(bp,4)=std(group2)/sqrt(20);
    byfreqstats(bp,5)=realEffect;
    byfreqstats(bp,6)=p_value;

end


% plot frequency by condition means
close all
errorbar(myfreqs,byfreqstats(:,1),byfreqstats(:,3),'linewidth',1,'Color',[0 0 0],'LineStyle',':','CapSize',3)
hold on
plot(myfreqs,byfreqstats(:,1),'Color',barColorMap{1},'Marker','.','MarkerFaceColor',barColorMap{1},'MarkerSize',28)%,'LineWidth',.05)
errorbar(myfreqs,byfreqstats(:,2),byfreqstats(:,4),'linewidth',1,'Color',[0 0 0],'LineStyle',':','CapSize',3)
plot(myfreqs,byfreqstats(:,2),'Color',barColorMap{2},'Marker','.','MarkerFaceColor',barColorMap{2},'MarkerSize',28)%,'LineWidth',.05)
yline(0,'LineWidth',1,'Color',[0 0 0]);ylim([-0.01 0.1]);xlim([0.5 10.5]);%xline(0);
set(gca,'xtick',myfreqs,'xticklabel',{myfreqs},...
    'ytick',0:0.02:0.08,'yticklabel',{0:0.02:0.08}),  grid on
xlabel('Frequency')
ylabel('Mean Decoder Score')
%legend(condlabels)
set(gcf, 'Position', [332,313,560,226]);





for c=1:length(condnames)
       %errorbar(myfreqs+(c/4),byfreqmeans(:,c),zeros(10,1),byfreqerr(:,c),'linewidth',1,'Color',[0 0 0],'LineStyle',':')
       errorbar(myfreqs,byfreqmeans(:,c),byfreqerr(:,c),'linewidth',1,'Color',[0 0 0],'LineStyle',':','CapSize',3)
    hold on
     plot(myfreqs,byfreqmeans(:,c),'Color',barColorMap{c},'Marker','.','MarkerFaceColor',barColorMap{c},'MarkerSize',28)%,'LineWidth',.05)
end
% sigstar(starbars,starsig)
yline(0,'LineWidth',1,'Color',[0 0 0]);ylim([-0.01 0.1]);xlim([0.5 10.5]);%xline(0);
set(gca,'xtick',myfreqs,'xticklabel',{myfreqs},...
    'ytick',0:0.02:0.08,'yticklabel',{0:0.02:0.08}),  grid on
xlabel('Frequency')
ylabel('Mean Decoder Score')
%legend(condlabels)
set(gcf, 'Position', [332,313,560,226]);




