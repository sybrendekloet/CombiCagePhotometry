%% Analysis program for fiber photometry data in 5-CSRTT
% Contents:
% 
% - Data import and setup
% 
%   > Data has already been pre-processed, and needs to be loaded from
%   corresponding .mat file
%
%   > In this part, we set all the relevant parameters for the next step.
%
% - Data structuring
%
%   > Pre-processed data has been structured in a particular way. This part
%   of the script generates two outputs:
%
%       1)  data_in: This is a cell array of 2 cells. Each cell contains a
%       matrix of all individual traces. Traces correspond to single
%       trials.
%
%       2)  meta_in: This is a matrix with all the relevant task parameters
%       (e.g. session type ['variable delay' (vITI) or 'variable cue'
%       (vSD)] or outcome ['Correct', 'Incorrect', 'Omission',
%       'Premature']). This allows for easy access to trials needed for
%       statistical analysis, graphing or any other use.
% 
% - Behavioral analysis
% 
% - Neurophysiological analysis

%% Import and basic parameters

% Set directory for analysis output (graphs, tables)
exptname = 'group_data_projections'; % name folder to which analysis will belong
savefolder = ['D:\Publications\Projection paper\Data\Analysis\' exptname];%folder where this analysis will be saved
mkdir(savefolder); % creates directory where figures will end up
saveFig = 1; % 1 = save all figures, 0 = don't save figures    

%% Classification
%
norm = 1; % Normalization type. 0 = no transformation, 1 = z-transform
dev = 0; % Deviation criterion. 0 = no selection. Other inputs = only trials 
% where amount of frames that are >[dev] from baseline is larger than 10.
% This allows selection of only trials where population is active. 

[data_in, meta_in] = classify_FPdata(compiledData, norm, dev);

%% Comparison between ITIs
 
% take correct responses only
% get AUC for each trt-cue
% compare AUCs for each group w ANOVA
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9

[aucData, aucTest, aucTraces] = compareAUC(data_in, meta_in);


%% Basic metrics

% for each rat, get all 12.5s trials and determine 20-80 and 10-90 rise time and
% take mean
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9
% trial time

[test_out,test_stats_pop] = signalParameters(data_in, meta_in);



%% Save all plots
figHandles = findobj('Type', 'figure'); % Gather all figure handles
set(figHandles, 'renderer', 'painters'); % Set renderer to 'painters' for proper further processing in Illustrator
set(figHandles, 'Position', get(0, 'Screensize')); % Set figures to full screen size (easier to process images in Illustrator)

for fi = 1:numel(figHandles)
    saveas(figHandles(fi), [savefolder, '\FPbars-',num2str(figHandles(fi).Number)], 'epsc'); % Save each figure as .EPSC file 
end

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 3:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = ['Fig ' num2str(get(FigHandle, 'Number'))];
%   saveas(FigHandle, [savefolder, '/',FigName], 'epsc');
% end
% close all

%% Population activity DURING ITI
% Is population active during task and does it respond to changing task
% parameters?
%
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9
[win_store_kt, win_store_dt, win_store_pt] = compareTraces(data_in, meta_in)

% b) ITI duration & AUC

%% Population activity AROUND CUE
% Is population active during task and does it respond to changing task
% parameters?
%
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9
% a) Bootstrap to baseline
data_in_boot_cue = cell(4,4); %bootstrap input array
combs_resp = [1 3; 1 4; 3 4]; % combinations for permutation test
combs_pop = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % combinations for permutation test
combs_cue = [1 3]; % combinations for permutation test
perm_out_resp = cell(1,4); %permutation test for trial outcome array
perm_pop_resp = cell(1,4); %permutation test for population array
f3 = figure;
f4 = figure;
set(f3, 'renderer', 'painters')
set(f4, 'renderer', 'painters')

% Ztransform the data (on single trial lvl)
% baseline window = -4s : -2s
% bwin = [ceil(1*15.89)^0, ceil(4*15.89)];
% data_cue = (data_in{2}-mean(data_in{1}(:,bwin(1):bwin(2)),2))./...
%     std(data_in{1}(:,bwin(1):bwin(2)),[],2);

% Only trials that show deviation from BL
% data_in_cue_z = data_cue(sum(data_in_z>2,2)>10,:);
% meta_in_act_cue_z = meta_in(sum(data_in_z>2,2)>10,:);
% Sync win == -3: 5 around response

data_cue = data_in_z_act;
meta_cue = meta_in_act;


for g = 1:max(unique(meta_cue(:,1)))
    
    for o = [1 3]
        
        % Categorize all data per rat
        for r =1:numel(unique(meta_cue(meta_cue(:,1)==g,2)))
            ratno = unique(meta_cue(meta_cue(:,1)==g,2));
            ratId = ratno(r);
            
            % distinguish between response types, iti, rat, etc.
            cfj = (meta_cue(:,1)==g & meta_cue(:,2)==ratId & meta_cue(:,3)==1 &...
                meta_cue(:,5)==o & meta_cue(:,6)==3);
            data_tests = data_cue(cfj,ceil(14.5*15.89):ceil(22.5*15.89));
            
            % get baseline per rat
            %             bm_mean = mean(data_in{1}(cfj,1:ceil(4*15.89)),'all');
            %             bm_std= std(data_in{1}(cfj,1:ceil(4*15.89)),0,'all');
            %             tmp_mean = mean((data_in{1}(cfj,:)-bm_mean)./bm_std);
            
            %             tmp_mean = mean((data_in{1}(cfj,:)-mean(data_in{1}(cfj,1:ceil(3*15.89)),2))./...
            %                 std(data_in{1}(cfj,1:ceil(3*15.89)),[],2));
            
            tmp_mean = mean(data_tests);
            
            % bin data
            bin_w = 3;
            [~,~,idx] = histcounts(1:numel(tmp_mean),1:bin_w:numel(tmp_mean));
            tmp_mean(idx==0)=[]; idx(idx==0)=[];
            data_in_boot_cue{g,o}(r,:) = accumarray(idx(:),tmp_mean,[],@mean);

        end
        
        % Plot means for each outcome
        figure(f3)
        subplot(4,4,(g-1)*4+o)
        imagesc(data_in_boot_cue{g,o})
        xlim([0 42])
        caxis([-2 2])
        colormap(flipud(brewermap([],'RdBu')))
        line([16 16], get(gca,'ylim'), 'color', 'k')
        if g==4 
            if o ==4
            colorbar
            end
        end
%         colormap(jet)
        
        figure(f4)
        subplot(2,4,g)
        xlim([0 42])
        plot(movmean(mean(data_in_boot_cue{g,o}),3), 'color','k','linewidth',2)
         hold on
        
        % Bootstrap test for each rat
        sig = 0.0005; %bootstrap significance level
        fit_boots = 1; %bootstrap fit parameter
        num_boots = 5000; %number of bootstrap iterations
        data_bootstrap_analysis = data_in_boot_cue{g,o};
        
        [bootsCI, ...
          ~, kbootsCI, ...
          ~, ~, ...
          ~, ~] = ...
            bootstrap_data_v201(data_bootstrap_analysis,num_boots,fit_boots,sig);
 
        % bootstrap window shaded
        plotmean = mean(data_in_boot_cue{g,o});
%         jbfill(1:size(kbootsCI,2),...
%             movmean(bootsCI(1,:),3),...
%             movmean(bootsCI(2,:),3), col_resp(o),col_resp(o),1,0.25);
%         hold on

        % sem shaded
        jbfill(1:size(movmean(plotmean,3),2),...
            movmean(plotmean,3)+movmean(std(data_in_boot_cue{g,o},[],1),3)./sqrt(size(data_in_boot_cue{g,o},1)),...
            movmean(plotmean,3)-movmean(std(data_in_boot_cue{g,o},[],1),3)./sqrt(size(data_in_boot_cue{g,o},1)),...
            col_resp(o),col_resp(o),1,0.25);
        hold on
            line([16 16], [-2 4.5], 'color', 'k')
        
        % Plot bootstrap windows > 0
        kt_frames = 1:size(bootsCI,2);
        kt_frames(kbootsCI(1,:)<=0)=nan;

        plot(kt_frames, 2.5+o*0.06*ones(1,length(kt_frames)), 'color', col_resp(o),'linewidth',2)
        
        subplot(2,4,o+4)
        xlim([0 42])
        % plot comparison between populations
        plot(movmean(mean(data_in_boot_cue{g,o}),3), 'color','k','linewidth',2)
        hold on
            line([16 16], [-2 4.5],'color',  'k')
        % boots interval shaded
%         jbfill(1:size(kbootsCI,2),...
%             movmean(bootsCI(1,:),3),...
%             movmean(bootsCI(2,:),3), col_resp(g),col_resp(g),1,0.25);
%         hold on
        % sem shaded
        jbfill(1:size(movmean(plotmean,3),2),...
            movmean(plotmean,3)+movmean(std(data_in_boot_cue{g,o},[],1),3)./sqrt(size(data_in_boot_cue{g,o},1)),...
            movmean(plotmean,3)-movmean(std(data_in_boot_cue{g,o},[],1),3)./sqrt(size(data_in_boot_cue{g,o},1)),...
            col_rep(g),col_rep(g),1,0.25);
         hold on
       
    end
    
    % Permutation and bootstrap test between outcomes
    % in: data_in_boot{g,o}
    for c=1:size(combs_cue,1)
        subplot(2,4,g)
        % Permutation tests
        [perm_out_resp{g}(c,:), ~]=permTest_array(data_in_boot_cue{g,combs_cue(c,1)},...
            data_in_boot_cue{g,combs_cue(c,2)}, 1000);
        
        % Bootstrap difference between outcomes
        boot_diff_in = data_in_boot_cue{g,combs_cue(c,1)}-...
            data_in_boot_cue{g,combs_cue(c,2)};
        
        [bootsCI, ~, kbootsCI, ~, ~, ~, ~] = ...
            bootstrap_data_v201(boot_diff_in,num_boots,fit_boots,sig);
        dt_frames = ceil(1:8*15.89/3);
        dt_frames(kbootsCI(1,:)<=0)=nan;
        
        plot(dt_frames, 4+c*0.12*ones(1,length(dt_frames)), 'color', col_resp(combs_cue(c,1)),'linewidth',2)
        plot(dt_frames, 4.03+c*0.12*ones(1,length(dt_frames)), 'color', col_resp(combs_cue(c,2)),'linewidth',2)
        
    end

    if g == 4
        for o = 1:4
            % Permutation test between populations
            for c=1:size(combs_pop,1)
                subplot(2,4,o+4)
                try
                    [perm_pop_resp{o}(c,:), ~]=permTest_array(data_in_boot_cue{combs_pop(c,1),o},...
                        data_in_boot_cue{combs_pop(c,2),o}, 3000);
                    
                    ptp_frames = ceil(1:8*15.89/3);
                    ptp_frames(perm_pop_resp{o}(c,:)>=0.05)=nan;
                    plot(ptp_frames, 3.5+c*0.12*ones(1,length(ptp_frames)), 'color', col_rep(combs_pop(c,1)),'linewidth',2)
                    plot(ptp_frames, 3.53+c*0.12*ones(1,length(ptp_frames)), 'color', col_rep(combs_pop(c,2)),'linewidth',2)
                end
            end
        end
    end

    
end

% b) ITI duration & AUC
%% Save all plots
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = ['Fig ' num2str(get(FigHandle, 'Number'))];
  saveas(FigHandle, [savefolder, '/',FigName,'_cue'], 'epsc');
end

%% Population activity AROUND RESPONSE
% Is population active during task and does it respond to changing task
% parameters?
%
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9
% a) Bootstrap to baseline
data_in_boot_rsp = cell(4,4); %bootstrap input array
combs_resp = [1 3; 1 4; 3 4]; % combinations for permutation test
combs_pop = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % combinations for permutation test
combs_rsp = [1 3]; % combinations for permutation test
perm_out_resp = cell(1,4); %permutation test for trial outcome array
perm_pop_resp = cell(1,4); %permutation test for population array
f5 = figure;
f6 = figure;
set(f5, 'renderer', 'painters')
set(f6, 'renderer', 'painters')

% Ztransform the data (on single trial lvl)
% baseline window = -4s : -2s
% bwin = [ceil(1*15.89)^0, ceil(4*15.89)];
% data_rsp = (data_in{1}-mean(data_in{1}(:,bwin(1):bwin(2)),2))./...
%     std(data_in{1}(:,bwin(1):bwin(2)),[],2);
% 
% % Only trials that show deviation from BL
% data_in_rsp_z = data_rsp(sum(data_in_z>2,2)>10,:);
% meta_in_act_rsp_z = meta_in(sum(data_in_z>2,2)>10,:);
% % Sync win == -3: 5 around response

data_rsp = data_in_z_resp_act;
meta_rsp = meta_in_resp_act;

for g = 1:max(unique(meta_rsp(:,1)))
    
    for o = [1 3 4]
        
        % Categorize all data per rat
        for r =1:numel(unique(meta_rsp(meta_rsp(:,1)==g,2)))
            ratno = unique(meta_rsp(meta_rsp(:,1)==g,2));
            ratId = ratno(r);
            
            % distinguish between response types, iti, rat, etc.
            cfj = (meta_rsp(:,1)==g & meta_rsp(:,2)==ratId & meta_rsp(:,3)==1 &...
                meta_rsp(:,5)==o & meta_rsp(:,6)==3);
            data_tests = data_rsp(cfj,:);
            
            % get baseline per rat
            %             bm_mean = mean(data_in{1}(cfj,1:ceil(4*15.89)),'all');
            %             bm_std= std(data_in{1}(cfj,1:ceil(4*15.89)),0,'all');
            %             tmp_mean = mean((data_in{1}(cfj,:)-bm_mean)./bm_std);
            
            %             tmp_mean = mean((data_in{1}(cfj,:)-mean(data_in{1}(cfj,1:ceil(3*15.89)),2))./...
            %                 std(data_in{1}(cfj,1:ceil(3*15.89)),[],2));
            
            tmp_mean = mean(data_tests);
            
            % bin data
            bin_w = 3;
            [~,~,idx] = histcounts(1:numel(tmp_mean),1:bin_w:numel(tmp_mean));
            tmp_mean(idx==0)=[]; idx(idx==0)=[];
            data_in_boot_rsp{g,o}(r,:) = accumarray(idx(:),tmp_mean,[],@mean);

        end
        
        % Plot means for each outcome
        figure(f5)
        subplot(4,4,(g-1)*4+o)
        imagesc(data_in_boot_rsp{g,o})
        xlim([0 42])
        caxis([-2 2])
        colormap(flipud(brewermap([],'RdBu')))
            line([16 16], get(gca,'ylim'), 'color', 'k')
        if g==4 
            if o ==4
            colorbar
            end
        end
%         colormap(jet)
        
        figure(f6)
        subplot(2,4,g)
        xlim([0 42])
        plot(movmean(mean(data_in_boot_rsp{g,o}),3), 'color','k','linewidth',2)
         hold on
        
        % Bootstrap test for each rat
        sig = 0.0005; %bootstrap significance level
        fit_boots = 1; %bootstrap fit parameter
        num_boots = 5000; %number of bootstrap iterations
        data_bootstrap_analysis = data_in_boot_rsp{g,o};
        
        [bootsCI, ...
          ~, kbootsCI, ...
          ~, ~, ...
          ~, ~] = ...
            bootstrap_data_v201(data_bootstrap_analysis,num_boots,fit_boots,sig);
 
        % bootstrap window shaded
        plotmean = mean(data_in_boot_rsp{g,o});
%         jbfill(1:size(kbootsCI,2),...
%             movmean(bootsCI(1,:),3),...
%             movmean(bootsCI(2,:),3), col_resp(o),col_resp(o),1,0.25);
%         hold on

        % sem shaded
        jbfill(1:size(movmean(plotmean,3),2),...
            movmean(plotmean,3)+movmean(std(data_in_boot_rsp{g,o},[],1),3)./sqrt(size(data_in_boot_rsp{g,o},1)),...
            movmean(plotmean,3)-movmean(std(data_in_boot_rsp{g,o},[],1),3)./sqrt(size(data_in_boot_rsp{g,o},1)),...
            col_resp(o),col_resp(o),1,0.25);
        hold on
            line([16 16], [-2 4.5], 'color', 'k')
        
        % Plot bootstrap windows > 0
        kt_frames = 1:size(bootsCI,2);
        kt_frames(kbootsCI(1,:)<=0)=nan;

        plot(kt_frames, 2.5+o*0.06*ones(1,length(kt_frames)), 'color', col_resp(o),'linewidth',2)
        
        subplot(2,4,o+4)
        xlim([0 42])
        % plot comparison between populations
        plot(movmean(mean(data_in_boot_rsp{g,o}),3), 'color','k','linewidth',2)
        hold on
            line([16 16], [-2 4.5],'color',  'k')
        % boots interval shaded
%         jbfill(1:size(kbootsCI,2),...
%             movmean(bootsCI(1,:),3),...
%             movmean(bootsCI(2,:),3), col_resp(g),col_resp(g),1,0.25);
%         hold on
        % sem shaded
        jbfill(1:size(movmean(plotmean,3),2),...
            movmean(plotmean,3)+movmean(std(data_in_boot_rsp{g,o},[],1),3)./sqrt(size(data_in_boot_rsp{g,o},1)),...
            movmean(plotmean,3)-movmean(std(data_in_boot_rsp{g,o},[],1),3)./sqrt(size(data_in_boot_rsp{g,o},1)),...
            col_rep(g),col_rep(g),1,0.25);
         hold on
       
    end
    
    % Permutation and bootstrap test between outcomes
    % in: data_in_boot{g,o}
    for c=1:size(combs_rsp,1)
        subplot(2,4,g)
        % Permutation tests
        [perm_out_resp{g}(c,:), ~]=permTest_array(data_in_boot_rsp{g,combs_rsp(c,1)},...
            data_in_boot_rsp{g,combs_rsp(c,2)}, 1000);
        
        % Bootstrap difference between outcomes
        boot_diff_in = data_in_boot_rsp{g,combs_rsp(c,1)}-...
            data_in_boot_rsp{g,combs_rsp(c,2)};
        
        [bootsCI, ~, kbootsCI, ~, ~, ~, ~] = ...
            bootstrap_data_v201(boot_diff_in,num_boots,fit_boots,sig);
        dt_frames = ceil(1:8*15.89/3);
        dt_frames(kbootsCI(1,:)<=0)=nan;
        
        plot(dt_frames, 4+c*0.12*ones(1,length(dt_frames)), 'color', col_resp(combs_rsp(c,1)),'linewidth',2)
        plot(dt_frames, 4.03+c*0.12*ones(1,length(dt_frames)), 'color', col_resp(combs_rsp(c,2)),'linewidth',2)
        
    end

    if g == 4
        for o = 1:4
            % Permutation test between populations
            for c=1:size(combs_pop,1)
                subplot(2,4,o+4)
                try
                    [perm_pop_resp{o}(c,:), ~]=permTest_array(data_in_boot_rsp{combs_pop(c,1),o},...
                        data_in_boot_rsp{combs_pop(c,2),o}, 3000);
                    
                    ptp_frames = ceil(1:8*15.89/3);
                    ptp_frames(perm_pop_resp{o}(c,:)>=0.05)=nan;
                    plot(ptp_frames, 3.5+c*0.12*ones(1,length(ptp_frames)), 'color', col_rep(combs_pop(c,1)),'linewidth',2)
                    plot(ptp_frames, 3.53+c*0.12*ones(1,length(ptp_frames)), 'color', col_rep(combs_pop(c,2)),'linewidth',2)
                end
            end
        end
    end

    
end

%% Save all plots
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = ['Fig ' num2str(get(FigHandle, 'Number'))];
  saveas(FigHandle, [savefolder, '/',FigName,'_rsp'], 'epsc');
end


%% 
%% Population activity DURING ITI FOR SINGLE RATS
% Is population active during task and does it respond to changing task
% parameters?
%
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9
% a) Bootstrap to baseline
data_in_boot = cell(4,4); %bootstrap input array
combs_resp = [1 3; 1 4; 3 4]; % combinations for permutation test
combs_pop = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % combinations for permutation test
perm_out = cell(1,4); %permutation test for trial outcome array
perm_pop = cell(1,4); %permutation test for population array
f1 = figure;
f2 = figure;
set(f1, 'renderer', 'painters')
set(f2, 'renderer', 'painters')


data_pop = data_in_z;
meta_pop = meta_in;


for g = 1:max(unique(meta_pop(:,1)))
    
    for o = [1 3 4]
        
        % Categorize all data per rat
        for r =1:numel(unique(meta_pop(meta_pop(:,1)==g,2)))
            ratno = unique(meta_pop(meta_pop(:,1)==g,2));
            ratId = ratno(r);            
            
            % distinguish between response types, iti, rat, etc.
            if o == 4 % Filter prematures
            cfj = (meta_pop(:,1)==g & meta_pop(:,2)==ratId & meta_pop(:,3)==1 &...
                meta_pop(:,5)==o & meta_pop(:,6)==3 & (meta_pop(:,end)-meta_pop(:,end-1))>7.5);
            else
            cfj = (meta_pop(:,1)==g & meta_pop(:,2)==ratId & meta_pop(:,3)==1 &...
                meta_pop(:,5)==o & meta_pop(:,6)==3);
            end
            data_tests = data_pop(cfj,:);
            
            % get baseline per rat
            %             bm_mean = mean(data_in{1}(cfj,1:ceil(4*15.89)),'all');
            %             bm_std= std(data_in{1}(cfj,1:ceil(4*15.89)),0,'all');
            %             tmp_mean = mean((data_in{1}(cfj,:)-bm_mean)./bm_std);
            
            %             tmp_mean = mean((data_in{1}(cfj,:)-mean(data_in{1}(cfj,1:ceil(3*15.89)),2))./...
            %                 std(data_in{1}(cfj,1:ceil(3*15.89)),[],2));
            
            
            
            
            % Bootstrap test for each rat
            sig = 0.001; %bootstrap significance level
            fit_boots = 1; %bootstrap fit parameter
            num_boots = 2000; %number of bootstrap iterations
            data_bootstrap_analysis = data_tests(:,ceil(5*15.89):ceil(12.5*15.89));
            
            [bootsCI, ...
                ~, kbootsCI, ...
                ~, ~, ...
                ~, ~] = ...
                bootstrap_data_v201(data_bootstrap_analysis,num_boots,fit_boots,sig);
         
            data_out_boot{g,o}(r,:) = bootsCI(1,:)>0;
        end
    end
end



%% Difference between response types
% Is population activity different leading up to different response types?
%
% a) Permutation tests or bootstrap to 0 for differences between populations
%
% Individual rat lvl: permutation tests or bootstrap for each response type
%
%% Predictive value of signal
% Does signal hold any predictive information about trial outcome?
%
% a) Use auROC to show classification of trial outcomes
% b) train model and use it to predict trial outcomes 
%