function TV = FxBTL_Detection(path)
FOV.plot_xmin =482; % pixel
FOV.plot_xmax = 1127; % pixel
FOV.plot_ymin = 69; % pixel
FOV.plot_ymax = 620; % pixel

% FOV.volume_xmin = 294;
% FOV.volume_xmax = 304;
% FOV.volume_ymin = 118;
% FOV.volume_ymax = 672;
FOV.volume_xmin = 1177;
FOV.volume_xmax = 1179;
FOV.volume_ymin = 69;
FOV.volume_ymax = 665;

FOV.t0 = 479;
FOV.t60 = 1046;
FOV.fs = (FOV.t60 - FOV.t0)/60;
FOV.t_offset = (FOV.plot_xmin-FOV.t0)/FOV.fs;

% folder_name = 'ResultSpiro';
volume.grictick_value = 1000; % mL 8L : 670 /  6L : 500  / 16: 1250 / 18 : 1330
th_absplotcolor = 140; % RGB
filtformat = 'jpg';
filename = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'}; 
for cnt = 1:length(filename)
    %% Volume grid detection
    rawData = importdata(strcat(path,'\BTL\',filename{cnt},'.',filtformat));
    image_grid = rawData;
    image_grid = rgb2gray(image_grid);
    image_grid = double(image_grid(FOV.volume_ymin:FOV.volume_ymax,FOV.volume_xmin:FOV.volume_xmax));
    image_grid = mean(image_grid');
    
    image_grid(image_grid<200) = 0;
    image_grid(image_grid>200) = 1;
    
    image_grid = ~image_grid.*(1:length(image_grid));
    volume.min_grid = find(image_grid>0,1,'first');
    volume.max_grid = find(image_grid>0,1,'last');
    
    image_grid(image_grid == 0) = [];
    image_grid = abs(diff(image_grid));
    image_grid(image_grid == 1) = [];  %figure; imagesc(rawData)
    volume.n_grid = length(image_grid);
    volume.tick_grid = (volume.max_grid-volume.min_grid)/volume.n_grid;
    volume.scale = volume.grictick_value/volume.tick_grid;
    
    %% Plot detection
    temp = rawData;
    for i = 72:628
        for j = 386:1200
            if (temp(i,j,1) == temp(i,j,2)) && (temp(i,j,2) == temp(i,j,3))
                temp(i,j,:) = 255;
            end
        end
    end
    temp = rgb2gray(temp);
    temp(temp>th_absplotcolor) = 255;
    
    % ext waveform
    temp = double(temp(FOV.plot_ymin:FOV.plot_ymax,FOV.plot_xmin:FOV.plot_xmax));
    temp(temp<255) = 1;
    temp(temp==255) = NaN;
    temp = temp .* repmat(1:size(temp,1),size(temp,2),1)';
    SpiroSignal = nanmean(temp);
    tp_noise = 1; %%
    tidal_for_filt = FxfindRespPeak(SpiroSignal(tp_noise+1:end));
    tidal_for_filt.idx_insp = tidal_for_filt.idx_insp + tp_noise;
    SpiroSignal_filt = cumsum(diff(SpiroSignal)-nanmean(diff(SpiroSignal(tidal_for_filt.idx_insp(2):tidal_for_filt.idx_insp(end-1)))))+SpiroSignal(1);
    SpiroSignal_filt = [SpiroSignal_filt nan];
    
    % cal TV
    tidal = FxfindRespPeak(-SpiroSignal_filt,(1:length(SpiroSignal_filt))/FOV.fs);
    t_insp = tidal.t_insp + FOV.t_offset;
    TV = tidal.volume_insp * volume.scale;
    RR = [NaN; tidal.RR];
    MV = [TV.*RR]*0.001;
    
    SpiroSignal_filt_save(:,cnt) = SpiroSignal_filt;

    %% Result
    figure('Units','normalized','Position',[0 0 1 1]);
    imagesc(rawData); axis image off; hold on;
    plot(FOV.plot_xmin:FOV.plot_xmax,SpiroSignal_filt+FOV.plot_ymin,'--.r');
    plot(tidal.idx_exp + FOV.plot_xmin,SpiroSignal_filt(tidal.idx_exp) + FOV.plot_ymin + 10,'k^','MarkerFaceColor','k');
    plot(tidal.idx_insp + FOV.plot_xmin,SpiroSignal_filt(tidal.idx_insp) + FOV.plot_ymin - 10,'rv','MarkerFaceColor','r');
    for i = 1:length(tidal.volume_insp)
        text(tidal.idx_insp(i) + FOV.plot_xmin,SpiroSignal_filt(tidal.idx_insp(i)) + FOV.plot_ymin - 25,num2str(i),'FontWeight','Bold','FontSize',13);
    end
    % plot TV value
    for i = 1:length(tidal.volume_insp)
        text(1270,80+(30*i),strcat('Breath',num2str(i),': ',num2str(round(TV(i))),' mL'),'FontWeight','Bold','FontSize',15);
    end
    % plot grid point
    for i = 1:volume.n_grid+1
        plot(FOV.plot_xmin,(volume.min_grid + volume.tick_grid.*(i-1)) + FOV.volume_ymin,'b.','markersize',12);
    end
    text(1600,80+(30*1),strcat('Average: ',num2str(round(mean(TV))),' mL'),'FontWeight','Bold','FontSize',15);
    text(1600,80+(30*2),strcat('Sum: ',num2str(round(sum(TV))),' mL'),'FontWeight','Bold','FontSize',15);
        
    % save figure file
    saveas(gcf,strcat(path,'\BTL\',filename{cnt},'_p.png'));

    % save TV RR MV data
    resultfilename = strcat(path, '\BTL.xlsx');
    clear temp;
    temp = {'Time','TV(mL)','RR(bpm)','MV(L/min)'};
    for cnt_table = 1:length(TV)
        temp(cnt_table+1,1) = {round(t_insp(cnt_table),2)};
        temp(cnt_table+1,2) = {round(TV(cnt_table),2)};
        temp(cnt_table+1,3) = {round(RR(cnt_table),2)};
        temp(cnt_table+1,4) = {round(MV(cnt_table),2)};
    end    
    xlswrite(resultfilename,temp,char(filename(cnt)));
    
    disp(strcat(num2str(cnt),'/',num2str(length(filename))));
end
% save waveform data
xlswrite(resultfilename,SpiroSignal_filt_save,'Waveform');
close all;

% Sheet1 delete
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(resultfilename);
% objExcel.Workbooks.Open('D:\OneDrive\2. Project\01. Research\Bilab work\사내임상데이터처리\20210126_RM\2021-01-26-신경훈1\BTL.xlsx');
try
      objExcel.ActiveWorkbook.Worksheets.Item('Sheet1').Delete;
catch
      % Do nothing.
end
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;


%% Sub functions
function [tidal] = FxfindRespPeak(data,time)
% input
%   data : EIT data 
%   time : EIT time vector
% output
%   tidal : tidal volume & time info

if nargin < 2
    time = 1:length(data);
end
% 
% fs = 1/mean(diff(time)); % EIT framerate
if size(time,1) == 1
    time = time';
end
if size(data,2) == 1
    data = data';
end

fs = 1/mean(diff(time));

%% detrend data  
% N = 1; 
% % Wn = 0.8;
% Wn = 5;
% [b,a] = butter(N, 2*Wn/fs,'low');
% Detrend_data = filtfilt(b,a,data);
% time = (1:length(data))/length(data)*60;

%%
% fs = 10;
ws = round(2*fs);
MAC = zeros(length(data),1);
MAC(1:ws) = mean(data(1:ws));
if MAC(ws) < data(ws)
    state_itc = 1;
else
    state_itc = -1;
end
margin = 0.1;
cnt_itc = ws;
cnt_peak = 1;
cnt_v = 1;
i_itc = 0;
roc = 0;

% figure; plot(data); hold on; h_mac = plot(MAC,'--m');
for i = (ws+1):length(data)
    if i < 5*ws + 1
        MAC(i) = mean(data((i-ws):i)) - state_itc*std(data((i-ws):i))*margin;
    else
        MAC(i) = mean(data((i-ws):i)) - state_itc*std(data((i-5*ws):i))*margin;
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

%%

locs_Peak2 = roc;
if data(locs_Peak2(1)) > data(locs_Peak2(2))
    locs_Peak2(1) = [];
end

if mod(length(locs_Peak2),2) == 1
    locs_Peak2(end) = [];
end

tidal_th = 0.3;
idx_rm = [];
% for i = 2:(length(locs_Peak2)/2)
%     pre_tidal = abs(Detrend_data(locs_Peak2(2*(i-1)-1)) - Detrend_data(locs_Peak2(2*(i-1))));
%     cur_tidal = abs(Detrend_data(locs_Peak2(2*i-1)) - Detrend_data(locs_Peak2(2*i)));
%     if (tidal_th*pre_tidal) > cur_tidal
%         idx_rm = [idx_rm i];
%     end
% end

if length(idx_rm) > 1
    locs_Peak2([(2*idx_rm-1) (2*idx_rm)]) = [];
end

% figure;
% h(1) = subplot(211);
% plot(data,'linewidth',1.5); hold on;
% plot(locs_Peak2(1:2:end),data(locs_Peak2(1:2:end)),'k^','MarkerFaceColor','k');
% plot(locs_Peak2(2:2:end),data(locs_Peak2(2:2:end)),'rv','MarkerFaceColor','r');
% plot(MAC,'--.m','linewidth',0.2);
% plot(i_itc,MAC(i_itc),'mo','MarkerFaceColor','m');

tidal.idx_exp = zeros(length(locs_Peak2),1);
tidal.idx_insp = zeros(length(locs_Peak2),1);
tidal.t_exp = zeros(length(locs_Peak2),1);
tidal.t_insp = zeros(length(locs_Peak2),1);
tidal.volume_insp = zeros(length(locs_Peak2),1);
tidal.volume_exp = zeros(length(locs_Peak2),1);
tidal.volume_residual = zeros(length(locs_Peak2),1);
tidal.interval = zeros(length(locs_Peak2),1);
tidal.td_insp = zeros(length(locs_Peak2),1);
tidal.td_exp = zeros(length(locs_Peak2),1);
tidal.td_tidal = zeros(length(locs_Peak2),1);
tidal.IE = zeros(length(locs_Peak2),1);
tidal.flow_insp = zeros(length(locs_Peak2),1);
tidal.flow_exp = zeros(length(locs_Peak2),1);
tidal.RR = zeros(length(locs_Peak2),1);

tidal.idx_exp = locs_Peak2(1:2:end)';
tidal.idx_insp = locs_Peak2(2:2:end)';
tidal.t_exp = time(tidal.idx_exp);
tidal.t_insp = time(tidal.idx_insp);
tidal.volume_insp = abs(data(tidal.idx_insp) - data(tidal.idx_exp))'; % modified by JG
tidal.volume_exp = abs(data(tidal.idx_insp(1:end-1)) - data(tidal.idx_exp(2:end)))'; % modified by JG
tidal.volume_residual = tidal.volume_insp(1:end-1) - tidal.volume_exp';
% tidal.interval = [0; (time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)))];
% tidal.td_insp = time(tidal.idx_insp) - time(tidal.idx_exp);
% tidal.td_exp = time(tidal.idx_exp(2:end)) - time(tidal.idx_insp(1:end-1));
% tidal.td_tidal = tidal.td_insp(1:end-1) + tidal.td_exp;
% tidal.IE = tidal.td_insp(1:end-1)./tidal.td_tidal;
% tidal.flow_insp = tidal.volume_insp./tidal.td_insp;
% tidal.flow_exp = tidal.volume_exp./tidal.td_exp;
tidal.RR = 60./(time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)));

% h(2) = subplot(212);
% stem(locs_Peak2(2:2:end), tidal.volume_insp,'LineWidth',3,'Marker','none');
% linkaxes(h,'x');
end
end