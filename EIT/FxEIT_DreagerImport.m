function [EIT_data] = FxEIT_DreagerImport(EIT_data_path)
cd(EIT_data_path);
dirlist = dir('.');
for i = 1:length(dirlist)
    data_name{i}=dirlist(i).name;
end
data_name(1:2) = [];

EIT_data = [];
for i = 1:length(data_name);
    file_path = strcat(EIT_data_path,'\',data_name{i});
    [vv, auxdata, stim] = eidors_readdata(file_path);
    EIT_data = [EIT_data vv];
    disp(['data import : ' num2str(i) ' / ' num2str(length(data_name))]);
end   
end