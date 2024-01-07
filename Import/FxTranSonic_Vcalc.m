function [TranSonic] = FxTranSonic_Vcalc(TranSonic)
fs = 100;
op.fs = 100; % 100 f/s
op.th_ws = 1; % 1 s
op.th_lv = 0.6; % 60%
op.th_fws = 0.2; % 0.2
op.th_flow = 0.15; % 20% th

temp_flow_filt = movmean(TranSonic.flow_raw,[3 3]);
% th_amp = movmax(temp_flow_filt,round([op.fs op.fs]*op.th_ws))*op.th_lv;
th_amp_max = movmean(movmax(temp_flow_filt,round([op.fs op.fs]*.5*op.th_ws)),[op.fs op.fs]*.5);
th_amp_min = movmin(temp_flow_filt,round([op.fs op.fs]*.5*op.th_ws));
th_amp = (th_amp_max-th_amp_min)*op.th_lv + th_amp_min;
[local_flag] = find(diff(temp_flow_filt > th_amp)==1)+1;
local_flag(1) = [];
local_flag_amp = temp_flow_filt(local_flag);
local_flag(local_flag_amp<0.01) = [];
dflow = diff(temp_flow_filt);
ddflow = diff(dflow);

% figure; 
% h(1) = subplot(311); plot(TranSonic.flow_raw); hold on; plot([0 length(TranSonic.flow_raw)], [0 0],'r');
% h(2) = subplot(312); plot(dflow); hold on; plot([0 length(dflow)], [0 0],'r');
% h(3) = subplot(313); plot(ddflow); hold on; plot([0 length(ddflow)], [0 0],'r');
% linkaxes(h,'x');

figure;
plot(TranSonic.flow_raw); hold on;
% plot(TranSonic.ECG_i,TranSonic.flow_raw(TranSonic.ECG_i),'rv');
plot(th_amp);
plot(local_flag,TranSonic.flow_raw(local_flag),'rv','markerfacecolor','r');
% plot(diff(local_flag))
clear mflow idx_mflow mflow tsys;
for cnt = 1:length(local_flag)-1
    % find peak point
    [mflow(cnt), idx_mflow(cnt)] = max(TranSonic.flow_raw(local_flag(cnt):local_flag(cnt+1)-20));
%     figure; plot(TranSonic.flow_raw(local_flag(cnt-1):local_flag(cnt+1)))
    idx_mflow(cnt) = idx_mflow(cnt) + local_flag(cnt)-1;
    mflow(cnt) = mflow(cnt) * TranSonic.Scale;
    
    % find start sys
%     temp_dflow = dflow((idx_mflow(cnt)-round(diff(idx_mflow(cnt-1:cnt))*op.th_fws)):(idx_mflow(cnt)-5));
%     temp_flow = TranSonic.flow_raw(local_flag(cnt)-round(diff(local_flag(cnt:cnt+1))*op.th_fws):local_flag(cnt));
    temp_dflow = dflow((idx_mflow(cnt)-50):(idx_mflow(cnt)-5));
    th_flow(cnt) = (max(temp_dflow) - min(abs(temp_dflow)))*op.th_flow + min(abs(temp_dflow));
    tsys(cnt) = find(temp_dflow<th_flow(cnt),1,'last');
    tsys(cnt) = tsys(cnt) + idx_mflow(cnt) - 50;
end
% figure;
% plot(TranSonic.flow_raw); hold on; plot(idx_mflow(2:end), TranSonic.flow_raw(idx_mflow(2:end)),'rv','markerfacecolor','r');
% plot(local_flag,TranSonic.flow_raw(local_flag),'kv','markerfacecolor','k');
% subplot(211)
% plot(temp_flow,'k.'); hold on; plot(temp_flow,'k');
% subplot(212)
% plot(temp_dflow,'k.'); hold on; plot(temp_dflow,'k');
% figure; plot(diff(tsys)); figure; plot(TranSonic.flow_raw)
for cnt = 1:length(tsys)-1
    temp = TranSonic.flow_raw(tsys(cnt):tsys(cnt+1));
    
    TranSonic.volume.idx_sys(cnt) = tsys(cnt);
    TranSonic.volume.idx_maxflow(cnt) = idx_mflow(cnt);
    TranSonic.volume.max_flow(cnt) = mflow(cnt);
    
    if min(temp) > 0
        TranSonic.volume.steady(cnt) = (min(temp) * length(temp));
    else
        TranSonic.volume.steady(cnt) = 0;
    end
    TranSonic.volume.forward(cnt) = sum(temp.*(temp>0)) - TranSonic.volume.steady(cnt);
    TranSonic.volume.backward(cnt) = -sum(temp.*(temp<0));
    
    [~, tp_minima] = min(TranSonic.flow_raw(idx_mflow(cnt):tsys(cnt+1)));
    tp_minima = tp_minima + idx_mflow(cnt) - tsys(cnt);
    temp(tp_minima:end) = [];
    TranSonic.volume.forward_2(cnt) = sum(temp(temp>=temp(1))-temp(1));
    TranSonic.volume.forward_3(cnt) = sum(temp(temp>=temp(1)));
    TranSonic.volume.forward_4(cnt) = sum(temp>=temp(1)) * TranSonic.volume.max_flow(cnt);
    TranSonic.volume.forward_5(cnt) = sum(temp>=temp(1)) * (TranSonic.volume.max_flow(cnt)-temp(1));
end
TranSonic.volume.forward = TranSonic.volume.forward * TranSonic.Scale * 1000;
TranSonic.volume.backward = TranSonic.volume.backward * TranSonic.Scale * 1000;
TranSonic.volume.steady = TranSonic.volume.steady * TranSonic.Scale * 1000;
TranSonic.volume.unit = 'mL';

% TranSonic.volume.t_hms = HV.EIT.t_hms(TranSonic.volume.idx_sys);
TranSonic.volume.disp_flow = "figure; plot(TranSonic.flow_raw); hold on; plot(TranSonic.volume.idx_maxflow,TranSonic.flow_raw(TranSonic.volume.idx_maxflow),'gv','markerfacecolor','g'); plot(TranSonic.volume.idx_sys,TranSonic.flow_raw(TranSonic.volume.idx_sys),'r^','markerfacecolor','r'); plot([0 length(TranSonic.flow_raw)],[0 0],'k');";
TranSonic.volume.disp_volume = "figure; h(1) = subplot(311); plot(TranSonic.volume.idx_sys,TranSonic.volume.forward,'.'); temp_ylim=ylim; ylim([0 temp_ylim(2)]); title('forward volume'); h(2) = subplot(312); plot(TranSonic.volume.idx_sys,TranSonic.volume.backward,'.'); ylim([0 temp_ylim(2)]); title('backward volume'); h(3) = subplot(313); plot(TranSonic.volume.idx_sys,TranSonic.volume.steady,'.'); ylim([0 temp_ylim(2)]); title('steady volume'); linkaxes(h,'x');";
TranSonic.volume.disp_volume_hms = "figure; h(1) = subplot(311); plot(HV.EIT.t_hms(TranSonic.volume.idx_sys),TranSonic.volume.forward,'.'); temp_ylim=ylim; ylim([0 temp_ylim(2)]); title('forward volume'); h(2) = subplot(312); plot(HV.EIT.t_hms(TranSonic.volume.idx_sys),TranSonic.volume.backward,'.'); ylim([0 temp_ylim(2)]); title('backward volume'); h(3) = subplot(313); plot(HV.EIT.t_hms(TranSonic.volume.idx_sys),TranSonic.volume.steady,'.'); ylim([0 temp_ylim(2)]); title('steady volume'); linkaxes(h,'x');";
TranSonic.volume.disp = 'eval(TranSonic.volume.disp_flow)';

% eval(TranSonic.volume.disp_volume)

figure; 
h(1) = subplot(411); plot(TranSonic.flow_raw); hold on; plot([0 length(TranSonic.flow_raw)],[0 0],'k'); plot(TranSonic.volume.idx_sys,TranSonic.flow_raw(TranSonic.volume.idx_sys),'r^','markerfacecolor','r');
h(2) = subplot(412); plot(TranSonic.volume.idx_sys,TranSonic.volume.forward,'.'); temp_ylim=ylim; ylim([0 temp_ylim(2)]); title('forward volume');
h(3) = subplot(413); plot(TranSonic.volume.idx_sys,TranSonic.volume.backward,'.'); ylim([0 temp_ylim(2)]); title('backward volume');
h(4) = subplot(414); plot(TranSonic.volume.idx_sys,TranSonic.volume.steady,'.'); ylim([0 temp_ylim(2)]); title('steady volume');
linkaxes(h,'x');

figure;
plot(TranSonic.flow_raw); hold on;
plot(th_amp);
% plot(local_flag,TranSonic.flow_raw(local_flag),'k.');
plot(TranSonic.volume.idx_maxflow,TranSonic.flow_raw(TranSonic.volume.idx_maxflow),'gv','markerfacecolor','g');
plot(TranSonic.volume.idx_sys,TranSonic.flow_raw(TranSonic.volume.idx_sys),'r^','markerfacecolor','r');
plot(TranSonic.ECG_i,TranSonic.flow_raw(TranSonic.ECG_i),'bo','markerfacecolor','b');
plot([0 length(TranSonic.flow_raw)],[0 0],'k');
legend({'flow' 'lv th' 'flow max' 'SYS' 'Rpeak' 'zero'})

end

