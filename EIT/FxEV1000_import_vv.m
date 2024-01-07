function [EV1000] = FxEV1000_import_vv(path_EV1000, delay)
% ex) delay 
% EIT: -103s
% EV1000: -150s
% delay = 253s (MP50 - EIT)

[xlsx_txt] = readtable(path_EV1000,'ReadVariableNames',0,'Sheet','TPTD');

cnt_EV1000 = 1;
cnt = 3;
while(1)
    if isempty(xlsx_txt.Var1{cnt})
        break;
    end
    temp_date = strsplit(xlsx_txt.Var1{cnt},'/');
    temp_time = strsplit(xlsx_txt.Var2{cnt},{':',' '});
    if strcmp(temp_time{4},'pm') && str2double(temp_time{1})<12
        EV1000.TPTD.t_hms(cnt_EV1000,1) = datetime(str2double(temp_date{3}),str2double(temp_date{1}),str2double(temp_date{2}), ...
                                str2double(temp_time{1})+12,str2double(temp_time{2}),str2double(temp_time{3}));
    else
        EV1000.TPTD.t_hms(cnt_EV1000,1) = datetime(str2double(temp_date{3}),str2double(temp_date{1}),str2double(temp_date{2}), ...
                                str2double(temp_time{1}),str2double(temp_time{2}),str2double(temp_time{3}));
    end
    EV1000.TPTD.iCO(cnt_EV1000,1) = str2double(xlsx_txt.Var6{cnt});
    EV1000.TPTD.iCI(cnt_EV1000,1) = str2double(xlsx_txt.Var7{cnt});
    EV1000.TPTD.iSV(cnt_EV1000,1) = str2double(xlsx_txt.Var8{cnt});
    EV1000.TPTD.iSVI(cnt_EV1000,1) = str2double(xlsx_txt.Var9{cnt});
    EV1000.TPTD.iSVR(cnt_EV1000,1) = str2double(xlsx_txt.Var10{cnt});
    EV1000.TPTD.iSVRI(cnt_EV1000,1) = str2double(xlsx_txt.Var11{cnt});
    EV1000.TPTD.GEDV(cnt_EV1000,1) = str2double(xlsx_txt.Var12{cnt});
    EV1000.TPTD.GEDI(cnt_EV1000,1) = str2double(xlsx_txt.Var13{cnt});
    EV1000.TPTD.ITBV(cnt_EV1000,1) = str2double(xlsx_txt.Var14{cnt});
    EV1000.TPTD.ITBI(cnt_EV1000,1) = str2double(xlsx_txt.Var15{cnt});
    EV1000.TPTD.EVLW(cnt_EV1000,1) = str2double(xlsx_txt.Var16{cnt});
    EV1000.TPTD.ELWI(cnt_EV1000,1) = str2double(xlsx_txt.Var17{cnt});
    EV1000.TPTD.PVPI(cnt_EV1000,1) = str2double(xlsx_txt.Var18{cnt});
    EV1000.TPTD.GEF(cnt_EV1000,1) = str2double(xlsx_txt.Var19{cnt});
    EV1000.TPTD.CFI(cnt_EV1000,1) = str2double(xlsx_txt.Var20{cnt});
    EV1000.TPTD.PR(cnt_EV1000,1) = str2double(xlsx_txt.Var21{cnt});
    EV1000.TPTD.MAP(cnt_EV1000,1) = str2double(xlsx_txt.Var22{cnt});
    EV1000.TPTD.CVP(cnt_EV1000,1) = str2double(xlsx_txt.Var22{cnt});
    EV1000.TPTD.BT(cnt_EV1000,1) = str2double(xlsx_txt.Var22{cnt});
    EV1000.TPTD.ivolume(cnt_EV1000,1) = str2double(xlsx_txt.Var4{cnt});
    cnt_EV1000 = cnt_EV1000 + 1;
    cnt = cnt + 1;
end

[xlsx_txt] = readtable(path_EV1000,'ReadVariableNames',0,'Sheet','EV1000');
op.nblank_num = find(strcmp(xlsx_txt.Var2,'Time'),1,'last')+1;
cnt_EV1000 = 1;
for i = 1+op.nblank_num:size(xlsx_txt,1)
    temp_date = strsplit(xlsx_txt.Var1{i},'/');
    temp_time = strsplit(xlsx_txt.Var2{i},{':',' '});
    if strcmp(temp_time{4},'pm') && str2double(temp_time{1})<12
        EV1000.APCO.t_hms(cnt_EV1000,1) = datetime(str2double(temp_date{3}),str2double(temp_date{1}),str2double(temp_date{2}), ...
                                str2double(temp_time{1})+12,str2double(temp_time{2}),str2double(temp_time{3}));
    else
        EV1000.APCO.t_hms(cnt_EV1000,1) = datetime(str2double(temp_date{3}),str2double(temp_date{1}),str2double(temp_date{2}), ...
                                str2double(temp_time{1}),str2double(temp_time{2}),str2double(temp_time{3}));
    end
    EV1000.APCO.CO(cnt_EV1000,1) = str2double(xlsx_txt.Var3{i});
    EV1000.APCO.CI(cnt_EV1000,1) = str2double(xlsx_txt.Var4{i});
    EV1000.APCO.SV(cnt_EV1000,1) = str2double(xlsx_txt.Var5{i});
    EV1000.APCO.SVI(cnt_EV1000,1) = str2double(xlsx_txt.Var6{i});
    EV1000.APCO.SVV(cnt_EV1000,1) = str2double(xlsx_txt.Var7{i});
    EV1000.APCO.PRV(cnt_EV1000,1) = str2double(xlsx_txt.Var8{i});
    EV1000.APCO.PR(cnt_EV1000,1) = str2double(xlsx_txt.Var13{i});
    EV1000.APCO.SYS(cnt_EV1000,1) = str2double(xlsx_txt.Var14{i});
    EV1000.APCO.DIA(cnt_EV1000,1) = str2double(xlsx_txt.Var15{i});
    EV1000.APCO.MAP(cnt_EV1000,1) = str2double(xlsx_txt.Var16{i});
    cnt_EV1000 = cnt_EV1000 + 1;
end
end

