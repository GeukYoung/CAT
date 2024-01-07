function [EIT_data] = FxEIT_txtImport(EIT_data_path,filename)
if nargin < 2
    filename = 'Scan';
end

cd(EIT_data_path);
dirlist = dir('.');
fid = fopen(strcat([EIT_data_path,'\',dirlist(3,1).name]));
temp = textscan(fid, '%f%f%f%f');
fclose(fid);

if size(temp{1,1},1) < 16^2
    ch = 8;
elseif size(temp{1,1},1)< 32^2
    ch = 16;
else
    ch = 32;
end
disp(['ch : ', int2str(ch)]);

cnt = 1;
stop = 0;
miss = 0;
cnt_stop = 0;
while(stop~=1)
    try
        try
            path = strcat(EIT_data_path,'\',int2str(cnt),filename,'.txt' );
            fid = fopen(path);
            temp = textscan(fid, '%f%f%f%f');
            fclose(fid);
        catch
            path = strcat(EIT_data_path,'\',filename,int2str(cnt),'.txt' );
            fid = fopen(path);
            temp = textscan(fid, '%f%f%f%f');
            fclose(fid);
        end
        temp2 = temp{1,3} + 1i.*temp{1,4};
        if size(temp2,1) > ch^2
            temp2(end) = [];  %%% for java data header line
            temp{1,3}(ch^2+1:end) = [];
        end
        EIT_data(:,cnt) = temp2;
        
        cnt = cnt + 1;
        cnt_stop = 0;
    catch
        miss = miss + 1;
        cnt_stop = cnt_stop + 1;
        if cnt_stop > 4
            stop = 1;
        end
    end
    disp(['data import : ' num2str(cnt) ' / ' num2str(length(dirlist))]);
end

disp(['miss data : ', int2str(miss-5)]);