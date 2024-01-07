function [sqi] = sqi_ESU(Data,temp,CVS)
if nargin < 2
    temp = 1;
end

switch temp
    case 1
        wd.CV = 100 * 0.2; % 200 ms
        wd.th_t = 100 * 60; % 60 s
        wd.th_sqi = 5;
        margin = 100*5; % 5 s additional remove from start of sqi bad

        temp4 = movstd(Data',[wd.CV,0])'./movmean(Data',[wd.CV,0])';
        cvmax = movmean(max(temp4),[wd.CV 0]);

        th_sqi = movmin(cvmax,[wd.th_t 0])*wd.th_sqi;
        sqi = cvmax<th_sqi; % 1: good, 0: bad
        sqi = movmin(sqi,[margin margin]);
    case 2
%         wd.CV = 100 * 0.2; % 200 ms
%         wd.th_t = 100 * 60; % 10 s
%         wd.th_sqi = 5;
%         margin = 100*5; % 5 s additional remove from start of sqi bad
% 
%         temp4 = movstd(Data',[wd.CV,0])'./movmean(Data',[wd.CV,0])';
%         cvmax = movmean(max(temp4),[wd.CV 0]);
% 
%         th_sqi = movmin(cvmax,[wd.th_t 0])*wd.th_sqi;
%         sqi_cv = cvmax<th_sqi; % 1: good, 0: bad
% %         sqi_cv = movmin(sqi_cv,[margin margin]);
%         
%         
%         N = 2; fc = 10; butter = 1; bpf = 2;
%         Data2 = FxEIT_Filter(Data,100,N,fc,butter,bpf);
%         figure; plot(Data2(1,:));
%         Data2_abs = abs(Data2);
%         Data2_abs_avg = movmean(Data2_abs',[20 0])';
%         Data2_abs_avg_max = mean(Data2_abs_avg);
% 
%         th_lv = 0.005;
%         figure;
%         h(1) = subplot(611); plot(Data2(1,:)); title('Z_1(HF)');
%         h(2) = subplot(612); plot(Data2_abs_avg(1,:)); title('abs(Z_1(HF))');
%         h(3) = subplot(613); plot(Data2_abs_avg_max); yline(th_lv,'r'); title('max(abs(Z_1(HF)))');
%         h(4) = subplot(614); plot(Data2_abs_avg_max<th_lv); ylim([-0.25 1.25]); title('SQI');
%         h(5) = subplot(615); plot(cvmax); hold on; plot(th_sqi,'r'); title('CV_{max}');
%         h(6) = subplot(616); plot(sqi_cv); ylim([-0.25 1.25]); title('SQI_{CV}');
%         linkaxes(h,'x');
%         xlim([6.7 7.8]*1e4);
%         
%         
%         figure; plot(Data2_abs(1,:)); hold on; plot(Data2_abs_avg(1,:));
%         figure; plot(Data2_abs_avg_max);
%         xlim([6.7 7.8]*1e4);
%         
%         %%
%         CVS2 = FxEIT_Filter(CVS,100,N,fc,butter,bpf);
%         
%         th_lv = 0.0001;
%         figure;
%         h(1) = subplot(611); plot(CVS); title('CVS');
%         h(2) = subplot(612); plot(CVS2); title('CVS(HF)');
%         h(3) = subplot(613); plot(movmean(CVS2.^2,[20 0])*50); yline(th_lv,'r'); title('moveavg(CVS^2)'); ylim([0 1]);
%         h(4) = subplot(614); plot(movmean(CVS2.^2,[20 0])*50<th_lv); ylim([-0.25 1.25]); title('SQI');
%         h(5) = subplot(615); plot(cvmax); hold on; plot(th_sqi,'r'); title('CV_{max}');
%         h(6) = subplot(616); plot(sqi); ylim([-0.25 1.25]); title('SQI_{CV}');
%         linkaxes(h,'x');
%         xlim([6.7 7.8]*1e4);
%         
%         th_lv = 0.1;
%         figure;
%         h(1) = subplot(611); plot(CVS); title('CVS');
%         h(2) = subplot(612); plot(CVS2); title('CVS(HF)');
%         h(3) = subplot(613); plot(movmean(abs(CVS2),[20 0])*50); yline(th_lv,'r'); title('moveavg(abs(CVS))');
%         h(4) = subplot(614); plot(movmean(abs(CVS2),[20 0])*50<th_lv); ylim([-0.25 1.25]); title('SQI');
%         h(5) = subplot(615); plot(cvmax); hold on; plot(th_sqi,'r'); title('CV_{max}');
%         h(6) = subplot(616); plot(sqi); ylim([-0.25 1.25]); title('SQI_{CV}');
%         linkaxes(h,'x');
%         xlim([6.7 7.8]*1e4);
%         
%         
% 
%         
%         xlim([62291.6802159209,365698.775689386]);
end
% sqi1 = sqi & circshift(sqi,-margin);
% sqi2 = sqi & circshift(sqi,-margin) & circshift(sqi,margin);

% figure; plot(sqi*1.05); hold on;
% plot(sqi1);  plot(sqi2*0.95); ylim([-0.25 1.25]);
% figure; plot(temp7); hold on; plot(th_sqi);
% figure; plot(temp7); hold on;
% plot(th_sqi); plot(movmin(temp7,[wd.th_t 0]) + movstd(temp7,[wd.th_t 0]));
% figure; plot(temp5');
% figure; plot(temp7');
% figure; boxplot(temp5');
% figure;
% h(1) = subplot(211); plot(temp_0');
% h(2) = subplot(212); plot(temp7'); hold on; plot()
% linkaxes(h,'x');

% d_sqi = diff(sqi);
% sqi_x = []; sqi_2 = [];
% for cnt = 1:length(sqi)-1
%     if d_sqi(cnt) == 1
%         sqi_x = [sqi_x cnt cnt];
%         sqi_2 = [sqi_2 0 1];
%     elseif d_sqi(cnt) == -1
%         sqi_x = [sqi_x cnt cnt];
%         sqi_2 = [sqi_2 1 0];
%     end
%     if cnt == length(sqi)-1
%         sqi_x = [sqi_x cnt];
%         sqi_2 = [sqi_2 sqi_2(end)];
%     end
% end

% figure; h(1) = subplot(211); hold on;
% area(sqi_x,sqi_2*5000,'FaceAlpha',0.6,'FaceColor','k','EdgeAlpha',0); ylim([0 20]);
% plot(Data(1:1:end,:)');
% h(2) = subplot(212); plot(temp7,'.'); hold on;
% plot(unique(sqi.*[1:length(sqi)]),[0 temp7(sqi)],'r.');
% plot(th_sqi,'g');
% % h(3) = subplot(313); plot(temp7)
% linkaxes(h,'x'); ylim([0 1]); xlim([1.64 1.7]*1e5)
end

