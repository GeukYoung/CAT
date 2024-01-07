% x축 y축 데이터 선택
n_EIT_point = length(Result2.EIT.time_hms);
n_ECG_point = length(Result2.ECG.time_hms);
n_point = min([n_EIT_point n_ECG_point]); % check EIT ECG data length missmatch
Result2.graph.xdata = Result2.EIT.time_hms(1:n_point);
Result2.graph.ydata(:,1) = Result2.EIT.RVS_inv(1:n_point);
Result2.graph.ydata(:,2) = Result2.EIT.TV(1:n_point);
Result2.graph.ydata(:,3) = Result2.EIT.RR(1:n_point);
Result2.graph.ydata(:,4) = Result2.EIT.MV(1:n_point);
Result2.graph.ydata(:,5) = Result2.EIT.SV(1:n_point);
Result2.graph.ydata(:,6) = Result2.EIT.HR(1:n_point);
Result2.graph.ydata(:,7) = Result2.EIT.CO(1:n_point);
Result2.graph.ydata(:,8) = Result2.EIT.dCVSdt_max(1:n_point);
Result2.graph.ydata(:,9) = Result2.EIT.MV(1:n_point).*Result2.EIT.CO(1:n_point);
Result2.graph.ydata(:,10) = Result2.EIT.MV(1:n_point)./Result2.EIT.CO(1:n_point);
Result2.graph.ydata(:,11) = Result2.EV1000_EIT.SYS(1:n_point);
Result2.graph.name = {'RVS0', 'TV', 'RR', 'MV', 'SV', 'HR', 'CO', 'dCVSdt' ,'MVxCO', 'MV/CO', 'SYS'};

%
sub_info.ytitle_pos = [-0.1 0.5];
sub_info.ytitle_margin = -0.005;
sub_info.DefautWidth = 0.85 + sub_info.ytitle_margin;
sub_info.label_fontsize = 13;
sub_info.axis_fontsize = 11;
sub_info.line_width = 2.5;
sub_info.subplot_posy_diff = -0.014;
sub_info.subplot_posy_offset = 0.05 + sub_info.subplot_posy_diff;
sub_info.n_ch = size(Result2.graph.ydata,2)+1;

figure('Units','normalized','position',[0.1,0.1,0.5,0.8]); clear H_ax;
for i = [1:3 (sub_info.n_ch-1):-1:4]
    H_ax{i} = subplot(sub_info.n_ch,1,i); hold on;
    temp_ypos = sub_info.subplot_posy_offset + sub_info.subplot_posy_diff*i;
    set(gca, 'Position', [ H_ax{i}.Position(1)-sub_info.ytitle_margin H_ax{i}.Position(2)+temp_ypos sub_info.DefautWidth H_ax{i}.Position(4)+0.02 ] );
    plot(Result2.graph.xdata,Result2.graph.ydata(:,i),'.k','markersize',12);
    
    ylabel(Result2.graph.name{i},'Fontsize',sub_info.label_fontsize,'FontWeight','Bold','Units','normalized','Position',sub_info.ytitle_pos,'Rotation',0,'VerticalAlignment','middle');
    box off; grid off;
    set(gca, 'Color', [ 0.99 0.99 0.99 ]);
    set(gca, 'GridColor', [ 0.3 0.3 0.3 ], 'GridAlpha', 0.1);
    set(gca, 'TickDir', 'out', 'YTick',ylim)
    set(gca, 'FontWeight','Bold','Fontsize',sub_info.axis_fontsize);
    set(gca, 'TickLength', [0.005, 1]);
    if i ~= sub_info.n_ch-1
        set(gca, 'XTickLabel','','XColor', 'w');
    end
end
linkaxes([H_ax{:}],'x');
xlim([Result2.graph.xdata(1)+minutes(1) Result2.graph.xdata(end)]);
