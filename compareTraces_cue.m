function [kt_frames, dt_frames, ptp_frames] = compareTraces_cue(data_in, meta_in)

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
% meta_in_cue_z = meta_in(sum(data_in_z>2,2)>10,:);
% Sync win == -3: 5 around response

data_cue = data_in{1};
meta_cue = meta_in;


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

end