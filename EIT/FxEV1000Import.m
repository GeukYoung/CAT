function [EV1000] = FxEV1000Import(path)
    if nargin < 1
        [path1, path2] = uigetfile('*.xls');
        path = [path2,path1];
    end
    [xlsx_txt] = readtable(path,'ReadVariableNames',0);

    op.nblank_num = find(strcmp(xlsx_txt.Var2,'Time'),1,'last')+1;
    cnt_EV1000 = 1;
    for i = 1+op.nblank_num:size(xlsx_txt,1)
        temp_date = strsplit(xlsx_txt.Var1{i},'.');
        temp_time = strsplit(xlsx_txt.Var2{i},{':',' '});
        if strcmp(temp_time{4},'pm') && str2double(temp_time{1})<12
            EV1000.t_hms(cnt_EV1000,1) = datetime(str2double(temp_date{3}),str2double(temp_date{2}),str2double(temp_date{1}), ...
                                    str2double(temp_time{1})+12,str2double(temp_time{2}),str2double(temp_time{3}));
        else
            EV1000.t_hms(cnt_EV1000,1) = datetime(str2double(temp_date{3}),str2double(temp_date{2}),str2double(temp_date{1}), ...
                                    str2double(temp_time{1}),str2double(temp_time{2}),str2double(temp_time{3}));
        end
        EV1000.CO(cnt_EV1000,1) = str2double(xlsx_txt.Var3{i});
        EV1000.CI(cnt_EV1000,1) = str2double(xlsx_txt.Var4{i});
        EV1000.SV(cnt_EV1000,1) = str2double(xlsx_txt.Var5{i});
        EV1000.SVI(cnt_EV1000,1) = str2double(xlsx_txt.Var6{i});
        EV1000.SVV(cnt_EV1000,1) = str2double(xlsx_txt.Var7{i});
        EV1000.PRV(cnt_EV1000,1) = str2double(xlsx_txt.Var8{i});
        EV1000.PR(cnt_EV1000,1) = str2double(xlsx_txt.Var13{i});
        EV1000.SYS(cnt_EV1000,1) = str2double(xlsx_txt.Var14{i});
        EV1000.DIA(cnt_EV1000,1) = str2double(xlsx_txt.Var15{i});
        EV1000.MAP(cnt_EV1000,1) = str2double(xlsx_txt.Var16{i});
        cnt_EV1000 = cnt_EV1000 + 1;
    end
end

