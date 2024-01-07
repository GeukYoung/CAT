function [MP50] = FxImport_MP50(path_MP50, delay)
% ex) delay 
% EIT: -103s
% MP50: -1s
% delay = 102s (MP50 - EIT)
if nargin < 2
    delay = seconds(0);
end

fid = fopen(path_MP50);
k = fscanf(fid,'%c',34);
fclose(fid);
k = strsplit(k,',');
k = strsplit(k{2},'.');
k = strsplit(k{1},' ');
YMD = strsplit(k{1},'-');
HMS = strsplit(k{2},':');
D = datetime(str2double(YMD{1}),str2double(YMD{2}),str2double(YMD{3}), ...
             str2double(HMS{1}),str2double(HMS{2}),str2double(HMS{3})) - delay;
X = convertTo(D,'epochtime','TicksPerSecond',1e3);

opts = detectImportOptions(path_MP50);
for cnt = 1:length(opts.VariableTypes)
    if strcmp(opts.VariableNames{cnt}, 'Event')
        opts.VariableTypes{cnt} = 'char';
        opts.VariableOptions(cnt).Type = 'char';
    else
        opts.VariableTypes{cnt} = 'double';
        opts.VariableOptions(cnt).Type = 'double';
    end
end

T = readtable(path_MP50,opts);
for cnt = 1:length(opts.VariableTypes)
    if strcmp(opts.VariableNames{cnt}, 'Time')
        X2 = X + int64(round(T.Time*1000));
        MP50.t_hms = datetime(X2','ConvertFrom','epochtime','TicksPerSecond',1e3)';
    elseif strcmp(opts.VariableNames{cnt}, 'Event')
        TF = ~ismissing(T.Event);
        MP50.Event.idx = unique([1:length(TF)]'.*TF);
        MP50.Event.idx(MP50.Event.idx==0) = [];
        MP50.Event.t_hms = MP50.t_hms(MP50.Event.idx);
        MP50.Event.tag = T.Event(MP50.Event.idx);
        
%         MP50.Event.t_hms = MP50.t_hms(TF);
%         MP50.Event.tag = T.Event(TF);
%         MP50.Event.idx = unique([1:length(TF)]'.*TF);
    else
        eval(['MP50.' opts.VariableNames{cnt} '= T.' opts.VariableNames{cnt} ';']);
    end
end
end

