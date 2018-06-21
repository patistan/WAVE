
expIDs = {'2087_1R_170328','2088_1R_170330','2099_1R1L_170327','2099_1R1L_170331',...
    '2105_1L_170410','2105_1L_170411','2110_1R_170406'};

home = pwd;
for s = 1:length(expIDs)
    %load data
    cd ..
    H = load(sprintf('%s_WAVE_classification_result_class_HYPERBOLICFREQUENCY_GREEDY_FINAL.mat',expIDs{s}));
    S = load(sprintf('%s_WAVE_classification_result_class_SPIRALFREQUENCY_GREEDY_FINAL.mat',expIDs{s}));
    cd(sprintf('H:/ProcessedDataArchive/Pati/73WAVE_matchedSF_fromSparrowhawk/%s_analysis/pati_analysis/Gauss_fit',expIDs{s}));
    s_fit = load('gaussFit_results.mat');
    cd ..
    W = load(sprintf('%s_WAVE_dataOut.mat',expIDs{s}));
    cd(home);
    %get cells needed
    s_cells_H = H.selected_cells{30};
    s_cells_S = S.selected_cells{30};
    s_intersect = length(intersect(s_cells_S,s_cells_H));
    %get rank of cells
    s_cells_H_rank = zeros(1,30);
    s_cells_S_rank = zeros(1,30);
    for i = 1:30
        s_cells_H_rank(ismember(s_cells_H,H.selected_cells{i})==1 & s_cells_H_rank==0)=i;
        s_cells_S_rank(ismember(s_cells_S,S.selected_cells{i})==1 & s_cells_S_rank==0)=i; 
    end
    %get bw
    intersect_H = intersect(W.dataOut.stats.global.responsive_cells_p001_fdr_average,s_cells_H)';
    s_cells_H_isSig = ismember(s_cells_H,intersect_H);
    s_cells_H_bw = s_fit.fit_oriBW(s_cells_H);    
    s_cells_H_bw_passed = s_cells_H_bw;
    s_cells_H_bw_passed(s_cells_H_isSig==0 | s_cells_H_bw<4 | s_cells_H_bw>90)=nan;
    
    intersect_S = intersect(W.dataOut.stats.global.responsive_cells_p001_fdr_average,s_cells_S)';
    s_cells_S_isSig = ismember(s_cells_S,intersect_S);
    s_cells_S_bw = s_fit.fit_oriBW(s_cells_S);    
    s_cells_S_bw_passed = s_cells_S_bw;
    s_cells_S_bw_passed(s_cells_S_isSig==0 | s_cells_S_bw<4 | s_cells_S_bw>90)=nan;
    
    %store values
    session.ID = expIDs{s};
    session.cellsH = s_cells_H;
    session.cellsH_rank = s_cells_H_rank;
    session.cellsS = s_cells_S;
    session.cellsS_rank = s_cells_S_rank;
    session.HSintersect = s_intersect;
    session.H_bw = s_cells_H_bw;
    session.H_sig = s_cells_H_isSig;
    session.H_bw_passed = s_cells_H_bw_passed;
    session.S_bw = s_cells_S_bw;    
    session.S_sig = s_cells_S_isSig;
    session.S_bw_passed = s_cells_S_bw_passed;
    
    session_data{s} = session;    
end

save('sessions_data.mat','session_data');

%% get BW for top 10 in each session
H_store_bw = [];
S_store_bw = [];
for s = 1:length(expIDs)
    H_ind = find(session_data{s}.cellsH_rank<11);
    H_bws = session_data{s}.H_bw(H_ind);
    H_sig = session_data{s}.H_sig(H_ind);
    H_bws(H_sig==0) = nan;
    
    S_ind = find(session_data{s}.cellsS_rank<11);
    S_bws = session_data{s}.S_bw(S_ind);
    S_sig = session_data{s}.S_sig(S_ind);
    S_bws(S_sig==0) = nan;
    
    H_store_bw = [H_store_bw H_bws];
    S_store_bw = [S_store_bw S_bws];
end

H_store_bw2 = H_store_bw;
H_store_bw2(H_store_bw<4 | H_store_bw>90) = nan;
sum(isnan(H_store_bw2))

S_store_bw2 = S_store_bw;
S_store_bw2(S_store_bw<4 |S_store_bw>90) = nan;
sum(isnan(S_store_bw))

%look at plots
figure
histogram(H_store_bw,[0:5:max(H_store_bw)+5])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Hyperbolic decoding cells (must have sig resp to grating)',sprintf('n=%i, median BW=%.2f',sum(~isnan(H_store_bw)),nanmedian(H_store_bw))})
saveas(gcf,'hist_decoding_hyp_group10_sigResp.fig')
saveas(gcf,'hist_decoding_hyp_group10_sigResp.png')

figure
histogram(S_store_bw,[0:5:max(S_store_bw)+5])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Spiral decoding cells (must have sig resp to grating)',sprintf('n=%i, median BW=%.2f',sum(~isnan(S_store_bw)),nanmedian(S_store_bw))})
saveas(gcf,'hist_decoding_spir_group10_sigResp.fig')
saveas(gcf,'hist_decoding_spir_group10_sigResp.png')

figure
histogram(H_store_bw2,[0:5:90])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Hyperbolic decoding cells (must have sig resp to grat and 4>BW<90)',sprintf('n=%i, median BW=%.2f',sum(~isnan(H_store_bw2)),nanmedian(H_store_bw2))})
saveas(gcf,'hist_decoding_hyp_group10_sigResp_passBW.fig')
saveas(gcf,'hist_decoding_hyp_group10_sigResp_passBW.png')

figure
histogram(S_store_bw2,[0:5:90])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Spiral decoding cells (must have sig resp to grat and 4>BW<90)',sprintf('n=%i, median BW=%.2f',sum(~isnan(S_store_bw2)),nanmedian(S_store_bw2))})
saveas(gcf,'hist_decoding_spir_group10_sigResp_passBW.fig')
saveas(gcf,'hist_decoding_spir_group10_sigResp_passBW.png')



%% get BW for top 30 in each session
H_store_bw = [];
S_store_bw = [];
H_store_rank = [];
S_store_rank = [];
for s = 1:length(expIDs)
    %H_ind = find(session_data{s}.cellsH_rank<11);
    H_bws = session_data{s}.H_bw;
    H_sig = session_data{s}.H_sig;
    H_rank = session_data{s}.cellsH_rank;
    H_bws(H_sig==0) = nan;
    H_rank(H_sig==0) = nan;
    
    %S_ind = find(session_data{s}.cellsS_rank<11);
    S_bws = session_data{s}.S_bw;
    S_sig = session_data{s}.S_sig;
    S_rank = session_data{s}.cellsS_rank;
    S_bws(S_sig==0) = nan;
    S_rank(S_sig==0) = nan;
    
    H_store_bw = [H_store_bw H_bws];
    H_store_rank = [H_store_rank H_rank];
    S_store_bw = [S_store_bw S_bws];
    S_store_rank = [S_store_rank S_rank];
end

%look at plots
figure
histogram(H_store_bw,[0:5:max(H_store_bw)+5])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Hyperbolic decoding cells (must have sig resp to grating)',sprintf('n=%i, median BW=%.2f',sum(~isnan(H_store_bw)),nanmedian(H_store_bw))})
saveas(gcf,'hist_decoding_hyp_group30_sigResp.fig')
saveas(gcf,'hist_decoding_hyp_group30_sigResp.png')

figure
histogram(S_store_bw,[0:5:max(S_store_bw)+5])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Spiral decoding cells (must have sig resp to grating)',sprintf('n=%i, median BW=%.2f',sum(~isnan(S_store_bw)),nanmedian(S_store_bw))})
saveas(gcf,'hist_decoding_spir_group30_sigResp.fig')
saveas(gcf,'hist_decoding_spir_group30_sigResp.png')

figure
hold on
for t = 1:30
    scatter(H_store_rank(H_store_rank==t),H_store_bw(H_store_rank==t))
end
xlabel('rank')
ylabel('Bandwidth')
title({'BW of Hyperbolic decoding cells by rank (must have sig resp to grating)',sprintf('n=%i, median BW=%.2f',sum(~isnan(H_store_bw)),nanmedian(H_store_bw))})
saveas(gcf,'scatter_decoding_hyp_group30_sigResp.fig')
saveas(gcf,'scatter_decoding_hyp_group30_sigResp.png')

figure
hold on
for t = 1:30
    scatter(S_store_rank(S_store_rank==t),S_store_bw(S_store_rank==t))
end
xlabel('rank')
ylabel('Bandwidth')
title({'BW of Spiral decoding cells by rank (must have sig resp to grating)',sprintf('n=%i, median BW=%.2f',sum(~isnan(S_store_bw)),nanmedian(S_store_bw))})
saveas(gcf,'scatter_decoding_spir_group30_sigResp.fig')
saveas(gcf,'scatter_decoding_spir_group30_sigResp.png')

%look at plots with removing bws outside of normal range
H_store_bw2 = H_store_bw;
H_store_bw2(H_store_bw<4 | H_store_bw>90) = nan;
H_store_rank2 = H_store_rank;
H_store_rank2(H_store_bw<4 | H_store_bw>90) = nan;
sum(isnan(H_store_bw2))

S_store_bw2 = S_store_bw;
S_store_bw2(S_store_bw<4 |S_store_bw>90) = nan;
S_store_rank2 = S_store_rank;
S_store_rank2(S_store_bw<4 |S_store_bw>90) = nan;
sum(isnan(S_store_bw))

save('sessions_data.mat','H_store_rank','H_store_bw','S_store_rank','S_store_bw',...
    'H_store_rank2','H_store_bw2','S_store_rank2','S_store_bw2','-append');

figure
histogram(H_store_bw2,[0:5:90])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Hyperbolic decoding cells (must have sig resp to grat and 4>BW<90)',sprintf('n=%i, median BW=%.2f',sum(~isnan(H_store_bw2)),nanmedian(H_store_bw2))})
saveas(gcf,'hist_decoding_hyp_group30_sigResp_passBW.fig')
saveas(gcf,'hist_decoding_hyp_group30_sigResp_passBW.png')

figure
histogram(S_store_bw2,[0:5:90])
xlabel('Bandwidth')
ylabel('# of cells')
title({'BW of Spiral decoding cells (must have sig resp to grat and 4>BW<90)',sprintf('n=%i, median BW=%.2f',sum(~isnan(S_store_bw2)),nanmedian(S_store_bw2))})
saveas(gcf,'hist_decoding_spir_group30_sigResp_passBW.fig')
saveas(gcf,'hist_decoding_spir_group30_sigResp_passBW.png')

figure
hold on
for t = 1:30
    scatter(H_store_rank2(H_store_rank2==t),H_store_bw2(H_store_rank2==t))
end
xlabel('rank')
ylabel('Bandwidth')
ylim([0 90])
title({'BW of Hyperbolic decoding cells by rank (sig resp to grat and 4>BW<90)',sprintf('n=%i, median BW=%.2f',sum(~isnan(H_store_bw2)),nanmedian(H_store_bw2))})
saveas(gcf,'scatter_decoding_hyp_group30_sigResp_passBW.fig')
saveas(gcf,'scatter_decoding_hyp_group30_sigResp_passBW.png')

figure
hold on
for t = 1:30
    scatter(S_store_rank2(S_store_rank2==t),S_store_bw2(S_store_rank2==t))
end
xlabel('rank')
ylabel('Bandwidth')
ylim([0 90])
title({'BW of Spiral decoding cells by rank (sig resp to grat and 4>BW<90)',sprintf('n=%i, median BW=%.2f',sum(~isnan(S_store_bw2)),nanmedian(S_store_bw2))})
saveas(gcf,'scatter_decoding_spir_group30_sigResp_passBW.fig')
saveas(gcf,'scatter_decoding_spir_group30_sigResp_passBW.png')

