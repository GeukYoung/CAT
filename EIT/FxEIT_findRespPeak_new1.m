function [tidal,data] = FxEIT_findRespPeak_new1(data,fs,lpf_fc)
% input
%   data : EIT data 
%   time : EIT time vector
% output
%   tidal : tidal volume & time info

% if nargin < 2
%     time = 1:length(data);
% end
% 
% fs = 1/mean(diff(time)); % EIT framerate

%% set values
global rules
rules.Peak_wsize = 60/12/2; % 60s, 12bpm, /2
rules.Peak_sdwsize = 3;
rules.Peak_margin = 0.2;

% rules.Peak_wsize = 60/12; % 60s, 12bpm, /2
% rules.Peak_sdwsize = 3;
% rules.Peak_margin = 0.2;

rules.th_iTV = 0.1;
rules.th_eTV = 0.25;
rules.n_thw = 3;
rules.th_ExpRatio = 0.5;
rules.th_InspRatio = 0.5;
rules.disp_nplot = 4;
rules.disp_cExpTV = 'g';
rules.disp_cInspTV = 'r';
rules.disp_cExpRatio = 'm';
rules.disp_cInspRatio = 'c';
rules.disp_marker = 'x';

%%
if size(data,2) == 1
    data = data';
end
time = ((0:(length(data)-1))/fs)';

%% detrend data 
if nargin == 4
    figure; plot(data); hold on;
    N = 1;
    Wn = lpf_fc;
    [b,a] = butter(N, 2*Wn/fs,'low');
    data1 = filtfilt(b,a,data);
    plot(data1,'r');
end

%% Peak detection
ws = round(rules.Peak_wsize*fs);
MAC = zeros(length(data),1);
MAC(1:ws) = mean(data(1:ws));
if MAC(ws) < data(ws)
    state_itc = 1;
else
    state_itc = -1;
end

cnt_itc = ws;
cnt_peak = 1;
cnt_v = 1;
i_itc = 0;
roc = 0;

% figure; plot(data); hold on; h_mac = plot(MAC,'--m');
for i = (ws+1):length(data)
    if i < rules.Peak_sdwsize*ws + 1
        MAC(i) = mean(data((i-ws):i)) - state_itc*std(data((i-ws):i))*rules.Peak_margin;
    else
        MAC(i) = mean(data((i-ws):i)) - state_itc*std(data((i-ws*rules.Peak_sdwsize):i))*rules.Peak_margin;
    end
%     h_mac.YData(i) = MAC(i); % plot

    if state_itc == 1
        if MAC(i) >= data(i)
            [~,temp_roc] = max(data((i-cnt_itc):i));
            roc(cnt_peak) = i - cnt_itc + temp_roc - 1;
            i_itc(cnt_peak) = i;
            cnt_peak = cnt_peak + 1;
            cnt_itc = 0;
            state_itc = -1;
            
%             plot(i,MAC(i),'mo','MarkerFaceColor','m');
%             plot(roc(end),data(roc(end)),'rv','MarkerFaceColor','r');
        end
    elseif state_itc == -1
        if MAC(i) <= data(i)
            [~,temp_roc] = min(data((i-cnt_itc):i));
            roc(cnt_peak) = i - cnt_itc + temp_roc - 1;
            i_itc(cnt_peak) = i;
            cnt_peak = cnt_peak + 1;
            cnt_itc = 0;
            state_itc = 1;
            
%             plot(i,MAC(i),'mo','MarkerFaceColor','m');
%             plot(roc(end),data(roc(end)),'k^','MarkerFaceColor','k');
        end
    end
    cnt_itc = cnt_itc + 1;
    
    if mod(i,5) == 0
%         xlim([i-30*fs i-1]);
%         drawnow;
    end
end

if data(roc(1)) > data(roc(2))
    roc(1) = [];
end

if mod(length(roc),2) == 1
    roc(end) = [];
end

locs_V = roc(1:2:end)'; locs_P = roc(2:2:end)';

figure;
plot(data,'linewidth',1.5); hold on;
plot(locs_V,data(locs_V),'k^','MarkerFaceColor','k');
plot(locs_P,data(locs_P),'gv','MarkerFaceColor','g');
plot(MAC,'--.m','linewidth',0.2);
plot(movmean(data,[ws 0]),'-m','linewidth',0.1);
plot(i_itc,MAC(i_itc),'mo','MarkerFaceColor','m');
% axis([126175.025201613,138706.527217742,335.772109340382,1570.01947863379]

figure;
if nargin > 2
    %% Volume remove
    % visual_1 peak point
    h(1) = subplot(rules.disp_nplot,1,1);
    plot(data,'linewidth',1.5); hold on;
    plot(locs_V,data(locs_V),'k^','MarkerFaceColor','k','MarkerSize',3);
    plot(locs_P,data(locs_P),'kv','MarkerFaceColor','k','MarkerSize',3);
    title('RVS'); box on
    % visual_2 debug view 
    h(2) = subplot(rules.disp_nplot,1,2); hold on;
    title('Removed Peak'); box on
    % visual_3 TV result
    h(3) = subplot(rules.disp_nplot,1,3);
    stem(locs_P,data(locs_P)-data(locs_V),'k.','linewidth',2);
    title('TV Result'); box on
    % visual_4 RR
    h(4) = subplot(rules.disp_nplot,1,4); hold on;
%     plot(locs_P(2:end),60./(diff(locs_P)/100),'.-k','MarkerSize',15,'linewidth',3);
    plot(locs_P(2:end),diff(locs_P)/100,'.-k','MarkerSize',15,'linewidth',3);

    title('RR Result'); box on
    linkaxes(h,'x');
    
    npre = length(locs_P);
    
%     % check insp ratio
%     [locs_P,locs_V] = FxRemove_InspRatio(locs_P,locs_V,data);
%     subplot(rules.disp_nplot,1,4); plot(locs_P(2:end),60./(diff(locs_P)/100),['.-' rules.disp_cInspRatio],'MarkerSize',14,'linewidth',2.5);
%     nremove = npre - length(locs_P)
    
    % check exp ratio
%     [locs_P,locs_V] = FxRemove_ExpRatio(locs_P,locs_V,data);
%     subplot(rules.disp_nplot,1,4); plot(locs_P(2:end),60./(diff(locs_P)/100),['.-' rules.disp_cExpRatio],'MarkerSize',14,'linewidth',2.5);
%     subplot(rules.disp_nplot,1,4); plot(diff(locs_P)/100,['.-' rules.disp_cExpTV],'MarkerSize',13,'linewidth',2);

    % check exp TV
    [locs_P,locs_V] = FxRemove_ExpTV(locs_P,locs_V,data);
%     subplot(rules.disp_nplot,1,4); plot(locs_P(2:end),60./(diff(locs_P)/100),['.-' rules.disp_cExpTV],'MarkerSize',13,'linewidth',2);
    subplot(rules.disp_nplot,1,4); plot(locs_P(2:end),diff(locs_P)/100,['.-' rules.disp_cExpTV],'MarkerSize',13,'linewidth',2);

    % check insp TV
    [locs_P,locs_V] = FxRemove_InspTV(locs_P,locs_V,data);
%     subplot(rules.disp_nplot,1,4); plot(locs_P(2:end),60./(diff(locs_P)/100),['.-' rules.disp_cInspTV],'MarkerSize',12,'linewidth',1.5);
    subplot(rules.disp_nplot,1,4); plot(locs_P(2:end),diff(locs_P)/100,['.-' rules.disp_cExpTV],'MarkerSize',13,'linewidth',2);

end

%% Params
tidal.idx_exp = locs_V;
tidal.idx_insp = locs_P;
tidal.t_exp = time(tidal.idx_exp);
tidal.t_insp = time(tidal.idx_insp);
tidal.volume_insp = abs(data(tidal.idx_insp) - data(tidal.idx_exp))'; % modified by JG
tidal.volume_exp = abs(data(tidal.idx_insp(1:end-1)) - data(tidal.idx_exp(2:end)))'; % modified by JG
tidal.volume_residual = tidal.volume_insp(1:end-1) - tidal.volume_exp;
tidal.interval = [0; (time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)))];
tidal.td_insp = time(tidal.idx_insp) - time(tidal.idx_exp);
tidal.td_exp = time(tidal.idx_exp(2:end)) - time(tidal.idx_insp(1:end-1));
tidal.td_tidal = tidal.td_insp(1:end-1) + tidal.td_exp;
tidal.IE = tidal.td_insp(1:end-1)./tidal.td_tidal;
tidal.flow_insp = tidal.volume_insp./tidal.td_insp;
tidal.flow_exp = tidal.volume_exp./tidal.td_exp;
tidal.RR = 60./(time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)));
tidal.RR1 = 60./tidal.interval;  
end

% temp_time_RR = tidal.t_insp;
% temp_RR = tidal.RR1;
% 
% temp_time_RR1 = zeros(length(temp_time_RR)*2,1);
% temp_time_RR1(1:2:end) = temp_time_RR;
% temp_time_RR1(2:2:end) = temp_time_RR+1e-10;
% temp_RR1 = zeros(length(temp_RR)*2,1);
% temp_RR1(1:2:end) = temp_RR;
% temp_RR1(2:2:end-1) = temp_RR(2:end);
% % 
% figure;
% h(1) = subplot(211); stem(tidal.t_insp, tidal.volume_insp,'k','linewidth',2,'Marker','none'); xlim([500 1200]);
% h(2) = subplot(212); plot(temp_time_RR1, temp_RR1,'k','linewidth',2); xlim([500 1200]);
% linkaxes(h,'x');


% 
% close all;
% figure; 
% h(1) = subplot(211); stem(tidal.volume_insp,'linewidth',2,'Marker','none'); xlim([0 255]);
% h(2) = subplot(212); plot(tidal.RR1); xlim([0 255]);
% linkaxes(h,'x');
% % 


% 
%   total_event_mv = event_Moderate+event_Severe;
%             temp_time_mv = ((mv_t-selected_time(t_start)+((hh*60*60)+(mm*60)+ss)-start_time)/60)-eit_delay;
%             temp_event_mv = total_event_mv;
%             
%             temp_time_mv2 = zeros(length(temp_time_mv)*2,1);
%             temp_time_mv2(1:2:end) = temp_time_mv;
%             temp_time_mv2(2:2:end) = temp_time_mv+1e-10;
%             temp_event_mv2 = zeros(length(temp_event_mv)*2,1);
%             temp_event_mv2(1:2:end) = temp_event_mv;
%             temp_event_mv2(2:2:end-1) = temp_event_mv(2:end);




function [locs_P,locs_V] = FxRemove_InspTV(locs_P,locs_V,data)
    global rules
    stop = 0;
    i = rules.n_thw+1;
    TV_insp = abs(data(locs_P) - data(locs_V));

    while stop == 0
        pre_tidal = median(TV_insp(i-rules.n_thw:i-1));
        cur_tidal = abs(data(locs_P(i)) - data(locs_V(i)));
        
        if (rules.th_iTV*pre_tidal) > cur_tidal
            % visual
            subplot(rules.disp_nplot,1,1);
            plot(locs_V(i),data(locs_V(i)),[rules.disp_cInspTV rules.disp_marker],'LineWidth',3);
            plot(locs_P(i),data(locs_P(i)),[rules.disp_cInspTV rules.disp_marker],'LineWidth',3);
            
            subplot(rules.disp_nplot,1,2);
            stem(locs_P(i),pre_tidal,'k.','linewidth',5);
            stem(locs_P(i),cur_tidal,[rules.disp_cInspTV rules.disp_marker],'linewidth',5);
          
            locs_V(i) = [];
            locs_P(i) = [];
            TV_insp(i) = [];
        else
            i = i + 1;
        end
        
        if i > length(TV_insp)
            stop = 1;
        end
    end
    subplot(rules.disp_nplot,1,3);
    stem(locs_P,data(locs_P)-data(locs_V),'k.','linewidth',2); title('TV Result');
end

function [locs_P,locs_V] = FxRemove_ExpTV(locs_P,locs_V,data)
    global rules
    stop = 0;
    i = rules.n_thw+1;
    TV_exp = abs(data(locs_P(1:end-1)) - data(locs_V(2:end)));
    
    while stop == 0
        pre_tidal = median(TV_exp(i-rules.n_thw:i-1));
        cur_tidal = abs(data(locs_P(i)) - data(locs_V(i+1)));
        
        if (rules.th_eTV*pre_tidal) > cur_tidal
            % visual
            subplot(rules.disp_nplot,1,1);
            plot(locs_V(i+1),data(locs_V(i+1)),[rules.disp_cExpTV rules.disp_marker],'LineWidth',3);
            plot(locs_P(i),data(locs_P(i)),[rules.disp_cExpTV rules.disp_marker],'LineWidth',3);
            subplot(rules.disp_nplot,1,2);
            stem(locs_V(i+1),pre_tidal,'k.','linewidth',5);
            stem(locs_V(i+1),data(locs_P(i))-data(locs_V(i+1)),[rules.disp_cExpTV rules.disp_marker],'linewidth',3);
            
            locs_V(i+1) = [];
            locs_P(i) = [];
            TV_exp(i) = [];
        else
            i = i + 1;
        end
        
        if i > length(TV_exp)-1
            stop = 1;
        end
    end
    subplot(rules.disp_nplot,1,3);
    stem(locs_P,data(locs_P)-data(locs_V),'k.','linewidth',2); title('TV Result');
end

function [locs_P,locs_V] = FxRemove_ExpRatio(locs_P,locs_V,data)
    global rules
    stop = 0;
    i = 1;
        
    while stop == 0
        TV_insp = abs(data(locs_P(i)) - data(locs_V(i)));
        TV_exp = abs(data(locs_P(i)) - data(locs_V(i+1)));
        
        if rules.th_ExpRatio > TV_exp/TV_insp
            % visual
            subplot(rules.disp_nplot,1,1);
            plot(locs_V(i+1),data(locs_V(i+1)),[rules.disp_cExpRatio rules.disp_marker],'LineWidth',3); hold on;
            plot(locs_P(i),data(locs_P(i)),[rules.disp_cExpRatio rules.disp_marker],'LineWidth',3);
            subplot(rules.disp_nplot,1,2);
            stem(locs_P(i),TV_insp,'k.','linewidth',5);
            stem(locs_P(i),TV_exp,[rules.disp_cExpRatio rules.disp_marker],'linewidth',3);
            
            locs_V(i+1) = [];
            locs_P(i) = [];
        else
            i = i + 1;
        end
        
        if i > length(locs_P)-1
            stop = 1;
        end
    end
    subplot(rules.disp_nplot,1,3);
    stem(locs_P,data(locs_P)-data(locs_V),'k.','linewidth',2); title('TV Result');
end

function [locs_P,locs_V] = FxRemove_InspRatio(locs_P,locs_V,data)
    global rules
    stop = 0;
    i = 2;
        
    while stop == 0
        TV_exp = abs(data(locs_P(i-1)) - data(locs_V(i)));
        TV_insp = abs(data(locs_P(i)) - data(locs_V(i)));
        
        if rules.th_InspRatio > TV_insp/TV_exp
            % visual
            subplot(rules.disp_nplot,1,1);
            plot(locs_V(i),data(locs_V(i)),[rules.disp_cInspRatio rules.disp_marker],'LineWidth',3); hold on;
            plot(locs_P(i),data(locs_P(i)),[rules.disp_cInspRatio rules.disp_marker],'LineWidth',3);
            subplot(rules.disp_nplot,1,2);
            stem(locs_P(i),TV_exp,'k.','linewidth',5);
            stem(locs_P(i),TV_insp,[rules.disp_cInspRatio rules.disp_marker],'linewidth',3);
            
            locs_V(i) = [];
            locs_P(i) = [];
        else
            i = i + 1;
        end
        
        if i > length(locs_P)-1
            stop = 1;
        end
    end
    subplot(rules.disp_nplot,1,3);
    stem(locs_P,data(locs_P)-data(locs_V),'k.','linewidth',2); title('TV Result');
end