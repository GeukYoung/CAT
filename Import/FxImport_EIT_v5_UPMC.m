function [Data] = FxImport_EIT_v5_UPMC(varargin)
% - inputs
% filepath: '...\RawData'
% block: 'first'/'last'/[2 3 4]
% mask: 'inable'/'disable'/[1 2 16 17...]
% type: 'mag'/'comp'

if nargin < 1
    varargin{1} = uigetdir; 
end

if ischar(varargin{1})
    Path_main = varargin{1};
    if ~contains(Path_main, '\RawData')
        temp_dir = dir(Path_main);
        if contains([temp_dir.name], 'RawData')
            Path_main = [Path_main '\RawData'];
        end
        clear temp_dir;
    end
    varargin = varargin(2:end);
%     argOffset = 1;
else
    error('ERROR: Data path');
end

opt.block = 0;
opt.mask = 1;
opt.type = 'magnitude';

numOrigInputArgs = numel(varargin);
if numOrigInputArgs ~= 0
    while numOrigInputArgs
        switch varargin{1}
            case 'block'
                if isnumeric(varargin{2})
                    opt.block = varargin{2}+1;
                else
                    if contains(varargin{2},'first')
                        opt.block = 1;
                    elseif contains(varargin{2},'last')
                        opt.block = 99;
                    end
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'mask'
                if contains(varargin{2},'disable')
                    opt.mask = 0;
                elseif contains(varargin{2},'inable')
                    opt.mask = 1;
                else
                    opt.mask = 1;
                    Data.ch_mask = varargin{2};
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'type'
                if contains(varargin{2},'mag') || contains(varargin{2},'abs')
                    opt.type = 'magnitude';
                elseif contains(varargin{2},'com')
                    opt.type = 'complex';
                elseif contains(varargin{2},'phase') || contains(varargin{2},'angle')
                    opt.type = 'phase';
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            otherwise
                msgbox(['Unvalid input :' varargin{1}]);
                numOrigInputArgs = 0;
        end
    end
end

% isCreatedUIfigure = false;
% if(nargin<2)
%     fig = uifigure;
%     isCreatedUIfigure = true;
% end

f = waitbar(0,'Please Wait...','Name','Loading...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

% dlg = uiprogressdlg(fig,'Title','Please Wait',...
%                     'Message','Loading EIT data','Indeterminate','on', 'Cancelable','on');
% cd(Path_main);
% dirlist = dir('.');

dirlist = dir(Path_main);

cnt.eit = 0;
cnt.voltage = 0;
cnt.phase = 0;
cnt.ecg = 0;
cnt.cm = 0;
cnt.setting = 0;
cnt.pressure = 0;

for cnt_file = 1:length(dirlist)
    temp_string = strsplit(dirlist(cnt_file).name,{'.', '['});
    if contains(temp_string{1},'EIT_','IgnoreCase',true)
        idx.eit(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.eit = cnt.eit + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).eit = fullfile(Path_main, dirlist(cnt_file).name);
        opt.flag_Z = 0;
    end
    
    if contains(temp_string{1},'EITImpedance_','IgnoreCase',true)
        idx.eit(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.eit = cnt.eit + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).eit = fullfile(Path_main, dirlist(cnt_file).name);
        opt.flag_Z = 1;
    end

    if contains(temp_string{1},'Voltage','IgnoreCase',true)
        idx.voltage(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.voltage = cnt.voltage + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).voltage = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Phase','IgnoreCase',true)
        idx.phase(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.phase = cnt.phase + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).phase = fullfile(Path_main,dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'ECG','IgnoreCase',true)
        idx.ecg(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.ecg = cnt.ecg + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).ecg = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'EITCM','IgnoreCase',true)
        idx.cm(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.cm = cnt.cm + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).cm = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Setting','IgnoreCase',true)
        idx.setting(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.setting = cnt.setting + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).setting = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Pressure','IgnoreCase',true)
        idx.pressure(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.pressure = cnt.pressure + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).pressure = fullfile(Path_main, dirlist(cnt_file).name);
    end
end

total_EIT = cnt.eit; clear cnt_file;

if cnt.eit > 0
    if sum(idx.eit) ~= length(idx.eit)
        msgbox(['EIT_' num2str(find(idx.eit==0)) 'bin data missing']);
    end
end
if cnt.voltage > 0
    if sum(idx.voltage) ~= length(idx.voltage)
        msgbox(['Voltage_' num2str(find(idx.voltage==0)) 'data missing']);
    end
end
if cnt.phase > 0
    if sum(idx.phase) ~= length(idx.phase)  
        msgbox(['Phase_' num2str(find(idx.phase==0)) 'data missing']);
    end
end
if cnt.ecg > 0
    if sum(idx.ecg) ~= length(idx.ecg)
        msgbox(['ECG_' num2str(find(idx.ecg==0)) 'data missing']);
    end
end
if cnt.cm > 0
    if sum(idx.cm) ~= length(idx.cm)
        msgbox(['EITCM_ ' num2str(find(idx.cm==0)) 'data missing']);
    end
end

if cnt.pressure > 0
    if sum(idx.pressure) ~= length(idx.pressure)
        msgbox(['Pressure_ ' num2str(find(idx.pressure==0)) 'data missing']);
    end
end

if opt.block ~= 0
    if size(opt.block) == 1
        switch opt.block
            case 1
                Path = Path(1);
                total_EIT = 1;
            case 99
                Path = Path(end);
                total_EIT = 1;
            otherwise
                Path = Path(opt.block);
                total_EIT = 1;
        end
    elseif length(opt.block) > 1         
        try
            Path_setting = Path(1).setting;
            Path = Path(opt.block);
            Path(1).setting = Path_setting;
            total_EIT = length(opt.block);
        catch
            msgbox('block data is not match');
        end
    end
end

% check data version
fid = fopen(Path(1).eit);
version = fread(fid,4,'uint8');
version = version(4)*2^3 + version(3)*2^2 + version(2)*2^1 + version(1)*2^0;
fclose(fid);

Data.version = version;
Data.path = Path_main;
Data.ch_mask = FxEIT_mask(16);
% Data.mask_2skip = FxEIT_mask(16,[16 1 3 4]);

% set(dlg,'Indeterminate','off');
switch version
%% v0000 for Old version
    case 0 % 0000
        header.meas_num = 256; % 256 measurement
        header.header_num = 72; % 72 byte header of 1 data
        header.data_num = 1536;
        
        header.amp = 5;
        header.gain = 13;
        header.ci_flag = 53;
        header.as_flag = 54;
        header.time = 55;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            if getappdata(f,'canceling')
                break;
            end
            
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            for cnt_scan2 = 1:scan_num
                if(mod(cnt_scan2, 1000) == 0)
                    waitbar(cnt_scan2/scan_num, ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)]);
                end
                
                if 0
                    temp = raw_data(:,cnt_scan2);
                    temp = reshape(temp,6,256);
                    temp(1,:)=temp(1,:)-128;
                    temp(1:2,:) = temp([2 1],:);
                    temp(3,:) = temp(4,:).*256 + temp(3,:);
                    temp(4,:) = temp(6,:).*256 + temp(5,:);
                    temp(5:6,:) = [];
                    temp(temp>32767) = temp(temp>32767) - 65536;
                    temp = temp';
                    Data.EIT_raw(:,cnt_scan) = sqrt(temp(:,3).^2 + temp(:,4).^2);
                end
                % header
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:4
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                hms.Y(cnt_scan) = double(typecast(uint8(temp2(header.time:header.time+1)),'uint16'));
                hms.M(cnt_scan) = temp2(header.time+2);
                hms.D(cnt_scan) = temp2(header.time+3);
                hms.H(cnt_scan) = temp2(header.time+4);
                hms.MI(cnt_scan) = temp2(header.time+5);
                hms.S(cnt_scan) = temp2(header.time+6);
                hms.MS(cnt_scan) = double(typecast(uint8(temp2(header.time+7:header.time+8)),'uint16'));
                
                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if 1
                if cnt.voltage >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).voltage);
                    raw_data = fread(fid,'double');
                    fclose(fid);
                    scan_num = size(raw_data,1)/header.meas_num;
                    if cnt_EIT == 1
                        Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                    else
                        Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                    end
                end
            end
            
            % phase
            if 0
                if cnt.voltage >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).phase);
                    raw_data = fread(fid,'double');
                    fclose(fid);
                    scan_num = size(raw_data,1)/header.meas_num;
                    if cnt_EIT == 1
                        Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
                    else
                        Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                    end
                end
            end
        end
        Data.t_hms = datetime(hms.Y, hms.M, hms.D, ...
            hms.H, hms.MI, hms.S, hms.MS);

%% v0010 for AirTom-R2
    case 2 % 0010
        header.meas_num = 256; % 256 measurement
        header.header_num = 168; % 168 byte header of 1 data
        header.data_num = 1536;
        
        header.amp = 5;
        header.gain = 13;
%         header.vsf = 45;
        header.ci_flag = 149;
        header.as_flag = 150;
        header.time = 151;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            if getappdata(f,'canceling')
                break;
            end
            
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            for cnt_scan2 = 1:scan_num
                
                if(mod(cnt_scan2, 1000) == 0)
                    if getappdata(f,'canceling')
                        break;
                    end
                    
                    if(mod(cnt_scan2, 1000) == 0)
                        waitbar(cnt_scan2/scan_num, ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)]);
                    end
                end
                
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM   L H     
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM   L H     
                temp(5:6,:) = [];
                temp(temp>32767) = temp(temp>32767) - 65536;
                temp = temp';
                Data.EIT_raw(:,cnt_scan) = sqrt(temp(:,3).^2 + temp(:,4).^2);
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:16
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                hms.Y(cnt_scan) = double(typecast(uint8(temp2(header.time:header.time+1)),'uint16'));
                hms.M(cnt_scan) = temp2(header.time+2);
                hms.D(cnt_scan) = temp2(header.time+3);
                hms.H(cnt_scan) = temp2(header.time+4);
                hms.MI(cnt_scan) = temp2(header.time+5);
                hms.S(cnt_scan) = temp2(header.time+6);
                hms.MS(cnt_scan) = double(typecast(uint8(temp2(header.time+7:header.time+8)),'uint16'));
                
                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).voltage);
                raw_data = fread(fid,'double');
                fclose(fid);
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
            % phase
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).phase);
                raw_data = fread(fid,'double');
                fclose(fid);
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                end
            end
        end
        Data.t_hms = datetime(hms.Y, hms.M, hms.D, ...
            hms.H, hms.MI, hms.S, hms.MS);
        
%% v0011 for AirTom/HemoVista
    case 3 % 0011
        header.meas_num = 256; % 256 measurement
        header.volt_num = 512; % 256 raw + 256 filt
        header.header_num = 168; % 168 byte header of 1 data
        header.data_num = 1536*2; % add filtered data
        header.fs_EIT = 100;
        header.scan_num = 360000;
        header.CI_ch = [2:17:256 241];
        
        header.amp = 5;
        header.gain = 13;
        header.idx_scan = 141;
        header.ci_flag = 149;
        header.as_flag = 150;
        header.time = 151+1;
        
        header.cm_size = 5;
        header.cm_scan = header.cm_size*16;
        header.cm_d2v = 0.0000221;
        header.cm_gain = 21201.06/4980;
        header.cm_Rs = 10;
        header.cm_d2i = header.cm_d2v/header.cm_gain/header.cm_Rs;
        
        header.ecg_size = 1;
        header.ecg_idx = 2;
        header.ecg_data = 6;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            if getappdata(f,'canceling')
                break;
            end
            
            % EIT raw
            if 1
                if(mod(cnt_scan2, 1000) == 0)
                    waitbar(cnt_EIT/total_EIT, ['Loading EIT_bin data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)]);
                end
                
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            
            if rem(size(raw_data,1)/(header.header_num + header.data_num),10000) == 0
                header.scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            end
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            if mod(size(raw_data,1),(header.header_num + header.data_num)) ~= 0
                raw_data(floor(scan_num)*(header.header_num + header.data_num)+1:end) = [];
                scan_num = floor(scan_num);
            end
            if cnt_EIT < total_EIT && scan_num ~= header.scan_num
                msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num) ' scan)']);
            end
            raw_data = reshape(raw_data,(header.header_num + header.data_num),scan_num);

            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            % header
            Data.EIT.Amp(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.amp:header.amp+7,:)),1,[]),'double');
            for cnt_gain=1:16
                Data.EIT.Gain(cnt_gain,cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7,:)),1,[]),'double');
            end
            Data.EIT.idx_scan(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.idx_scan:header.idx_scan+7,:)),1,[]),'uint64');
            Data.EIT.CI_flag(cnt_scan:cnt_scan+scan_num-1) = header_data(header.ci_flag,:);
            Data.EIT.AS_flag(cnt_scan:cnt_scan+scan_num-1) = header_data(header.as_flag,:);
            epochtime(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.time:header.time+7,:)),1,[]),'uint64');
            clear header_data
%             Data.ECG. HR
            end
            
            % raw_data
            if 1
                raw_data = raw_data(:);
                V_R = raw_data(3:6:end) + raw_data(4:6:end).*256;
                V_R(V_R>32767) = V_R(V_R>32767) - 65536;
                V_I = raw_data(5:6:end) + raw_data(6:6:end).*256;
                V_I(V_I>32767) = V_I(V_I>32767) - 65536;
                raw_data = V_R + V_I*1i;
                raw_data = reshape(raw_data,header.volt_num,length(raw_data)/header.volt_num);

                if 1 % mag
                    Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data(1:header.volt_num/2,:));
                    Data.EIT.raw_filt(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data(header.volt_num/2+1:header.volt_num,:));
                end

%               if 1 % complex
%                   Data.EIT_raw(:,cnt_scan:cnt_scan+scan_num-1) = raw_data(1:header.volt_num/2,:);
%                   Data.EIT_raw_filt(:,cnt_scan:cnt_scan+scan_num-1) = raw_data(header.volt_num/2+1:header.volt_num,:);
%               end
                clear raw_data;
            end
            
            % cm
            if 1
%                 scan_num = 360000;
                if cnt.cm >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).cm);
                    raw_data = fread(fid,'int16');
                    fclose(fid);
                    
                    scan_num_cm = size(raw_data,1)/header.cm_scan;
                    if mod(size(raw_data,1),header.meas_num) ~= 0 || scan_num_cm ~= scan_num
                        if cnt_EIT < total_EIT
                            disp('cm data miss');
                        end
                        if scan_num_cm > scan_num
                            disp(['EITCM data larger than raw(' num2str(scan_num_cm-scan_num) 'scans)']);
                            raw_data(scan_num*header.cm_scan+1:end) = [];
                            scan_num_cm = scan_num;
                        else
                            raw_data(scan_num_cm*header.cm_scan+1:scan_num*header.cm_scan) = NaN;
                            scan_num_cm = scan_num;
                        end
%                         raw_data(floor(scan_num_vol)*header.meas_num+1:end) = [];
%                         scan_num_vol = floor(scan_num_vol);
                     end
                    if cnt_EIT < total_EIT && scan_num_cm ~= header.scan_num
                        msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num_cm) ' scan)']);
                    end
                    
                    Data.EIT.cm_src(:,cnt_scan:cnt_scan+scan_num-1) = ...
                        reshape((raw_data(2:header.cm_size:end).^2 + raw_data(3:header.cm_size:end).^2).^0.5,16,scan_num) *  header.cm_d2i;
                    Data.EIT.cm_snk(:,cnt_scan:cnt_scan+scan_num-1) = ...
                        reshape((raw_data(4:header.cm_size:end).^2 + raw_data(5:header.cm_size:end).^2).^0.5,16,scan_num) *  header.cm_d2i;
                    Data.EIT.cm_mean(:,cnt_scan:cnt_scan+scan_num-1) = (Data.EIT.cm_src(:,cnt_scan:cnt_scan+scan_num-1) + Data.EIT.cm_snk(:,cnt_scan:cnt_scan+scan_num-1)) * 0.5;
                    clear raw_data;
                end
            end
            
            % voltage
            if 1
                if cnt.voltage >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).voltage);
                    raw_data = fread(fid,'double');
                    fclose(fid);
                    scan_num_vol = size(raw_data,1)/header.meas_num;
                    if mod(size(raw_data,1),header.meas_num) ~= 0 || scan_num_vol ~= scan_num
                        if cnt_EIT < total_EIT
                            disp('voltage data miss');
                        end
                        if scan_num_vol > scan_num
                            disp(['voltage data larger than raw(' num2str(scan_num_vol-scan_num) 'scans)']);
                            scan_num_vol = scan_num;
                     raw_data(header.meas_num*scan_num+1:end) = [];
                        end
%                         raw_data(floor(scan_num_vol)*header.meas_num+1:end) = [];
%                         scan_num_vol = floor(scan_num_vol);
                        raw_data(scan_num_vol*header.meas_num+1:scan_num*header.meas_num) = NaN;
                     end
                    if cnt_EIT < total_EIT && scan_num_vol ~= header.scan_num
                        msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num_vol) ' scan)']);
                    end
                    Data.EIT.volt(:,cnt_scan:cnt_scan+scan_num-1) = reshape(raw_data,header.meas_num,scan_num);
                    clear raw_data;
                end
            end
            
            % phase
            if 0
                if cnt.voltage >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).phase);
                    raw_data = fread(fid,'double');
                    fclose(fid);
                    scan_num_phase = size(raw_data,1)/header.meas_num;
                    if mod(size(raw_data,1),header.meas_num) ~= 0 || scan_num_phase ~= scan_num
                        if cnt_EIT < total_EIT
                            disp('phase data miss');
                        end
                        if scan_num_phase > scan_num
                            scan_num_phase = scan_num;
                            disp(['voltage data larger than raw(' num2str(scan_num_phase-scan_num) 'scans']);
                        end
%                         raw_data(floor(scan_num_phase)*header.meas_num+1:end) = [];
%                         scan_num_phase = floor(scan_num_phase);
                        raw_data(scan_num_phase*header.meas_num+1:scan_num*header.meas_num) = NaN;
                    end
                    if cnt_EIT < total_EIT && scan_num_phase ~= header.scan_num
                            msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num_phase) ' scan)']);
                    end
                    Data.EIT.phase(:,cnt_scan:cnt_scan+scan_num_phase-1) = reshape(raw_data,header.meas_num,scan_num_phase);
                    clear raw_data;
                end
            end
            
            % ecg
            if 1
                if cnt.ecg >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).ecg);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    
                    if mod(length(raw_data),37) ~= 0
                        if cnt_EIT < total_EIT
                            disp('ECG data miss');
                        end
                        raw_data(end-mod(length(raw_data),37)+1:end) = [];
                    end
                    raw_data = reshape(raw_data,37,[]);
                    Data.ECG.idx_ecg(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                    Data.ECG.ECG_fs(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = raw_data(header.ecg_size,:) * header.fs_EIT;
                    Data.ECG.ECG_hv(:,(cnt_scan-1)*Data.ECG.ECG_fs/header.fs_EIT+1:(cnt_scan+size(raw_data,2)-1)*Data.ECG.ECG_fs/header.fs_EIT) = ...
                        double(typecast(reshape(uint8(raw_data(header.ecg_data:header.ecg_data+2*Data.ECG.ECG_fs(cnt_scan)/header.fs_EIT-1,:)),1,[]),'uint16'));
                    clear raw_data;
                end
            end
            cnt_scan = cnt_scan + scan_num;
        end
        Data.EIT.t_hms = datetime(epochtime,'ConvertFrom','epochtime','TicksPerSecond',1e3);
        
        
        % Calc CI
        if sum(Data.EIT.CI_flag == 1) > 0
            [~, CI_idx]=find(Data.EIT.CI_flag==1);
            Data.EIT.CI=Data.EIT.volt(header.CI_ch,CI_idx)./Data.EIT.Amp(CI_idx); % contact impedance
        end
        
        % volt to ohm
        if isfield(Data.EIT, 'volt')
            Data.EIT.ohm = Data.EIT.volt./repmat(Data.EIT.Amp,256,1);
            Data.EIT.ohm(isinf(Data.EIT.ohm)) = NaN;
        end

%% v0100 for AirTom/HemoVista
    case 4 % 0100
        header.meas_num = 256; % 256 measurement
        header.volt_num = 512; % 256 raw + 256 filt
        header.header_num = 2776; % 2778 byte header of 1 data
        header.data_num = 1536*2; % add filtered data
        header.fs_EIT = 100;
        header.scan_num = 360000;
        header.CI_ch = [2:17:256 241];
        
        header.amp = 5;
        header.gain = 13;
        header.idx_scan = 141;
        header.ci_flag = 149;
        header.as_flag = 150;
        header.Codkey = 167;
        header.time = 152; % bug fixed(x) 151 + 1

        header.cm_size = 5;
        header.cm_scan = header.cm_size*16;
        header.cm_d2v = 0.0000221;
        header.cm_gain = 21201.06/4980;
        header.cm_Rs = 10;
        header.cm_d2i = header.cm_d2v/header.cm_gain/header.cm_Rs;
        
        header.ecg_size = 1;
        header.ecg_idx = 2;
        header.ecg_data = 6;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            if getappdata(f,'canceling')
                break;
            end
            
            % EIT raw
            if 1
                waitbar(cnt_EIT/total_EIT, ['Loading EIT_bin data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)]);
            
                fid = fopen(Path(cnt_EIT).eit);
                raw_data = fread(fid,'uint8');
                fclose(fid);

                if rem(size(raw_data,1)/(header.header_num + header.data_num),10000) == 0
                    header.scan_num = size(raw_data,1)/(header.header_num + header.data_num);
                end
                scan_num = size(raw_data,1)/(header.header_num + header.data_num);
                if mod(size(raw_data,1),(header.header_num + header.data_num)) ~= 0
                    raw_data(floor(scan_num)*(header.header_num + header.data_num)+1:end) = [];
                    scan_num = floor(scan_num);
                end
                if cnt_EIT < total_EIT && scan_num ~= header.scan_num
                    msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num) ' scan)']);
                end
                raw_data = reshape(raw_data,(header.header_num + header.data_num),scan_num);

                header_data=raw_data(1:header.header_num,:);
                raw_data(1:header.header_num,:) = [];

                % header
                Data.EIT.Amp(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.amp:header.amp+7,:)),1,[]),'double');
                for cnt_gain=1:16
                    Data.EIT.Gain(cnt_gain,cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7,:)),1,[]),'double');
                end
                Data.EIT.idx_scan(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.idx_scan:header.idx_scan+7,:)),1,[]),'uint64');
                Data.EIT.CI_flag(cnt_scan:cnt_scan+scan_num-1) = header_data(header.ci_flag,:);
                Data.EIT.AS_flag(cnt_scan:cnt_scan+scan_num-1) = header_data(header.as_flag,:);
                epochtime(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(header_data(header.time:header.time+7,:)),1,[]),'uint64');
                clear header_data
    %             Data.ECG. HR
            end
            
            % raw_data
            if 1
                raw_data = raw_data(:);
                V_R = raw_data(3:6:end) + raw_data(4:6:end).*256;
                V_R(V_R>32767) = V_R(V_R>32767) - 65536;
                V_I = raw_data(5:6:end) + raw_data(6:6:end).*256;
                V_I(V_I>32767) = V_I(V_I>32767) - 65536;
                raw_data = V_R + V_I*1i;
                raw_data = reshape(raw_data,header.volt_num,length(raw_data)/header.volt_num);

                if 1 % mag
                    Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data(1:header.volt_num/2,:));
                    Data.EIT.raw_filt(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data(header.volt_num/2+1:header.volt_num,:));
                end

%               if 1 % complex
%                   Data.EIT_raw(:,cnt_scan:cnt_scan+scan_num-1) = raw_data(1:header.volt_num/2,:);
%                   Data.EIT_raw_filt(:,cnt_scan:cnt_scan+scan_num-1) = raw_data(header.volt_num/2+1:header.volt_num,:);
%               end
                clear raw_data;
            end
            
            % cm
            if 1
%                 scan_num = 360000;
                if cnt.cm >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).cm);
                    raw_data = fread(fid,'int16');
                    fclose(fid);
                    
                    scan_num_cm = size(raw_data,1)/header.cm_scan;
                    if mod(size(raw_data,1),header.meas_num) ~= 0 || scan_num_cm ~= scan_num
                        if cnt_EIT < total_EIT
                            disp('cm data miss');
                        end
                        if scan_num_cm > scan_num
                            disp(['EITCM data larger than raw(' num2str(scan_num_cm-scan_num) 'scans)']);
                            raw_data(scan_num*header.cm_scan+1:end) = [];
                            scan_num_cm = scan_num;
                        else
                            raw_data(scan_num_cm*header.cm_scan+1:scan_num*header.cm_scan) = NaN;
                            scan_num_cm = scan_num;
                        end
%                         raw_data(floor(scan_num_vol)*header.meas_num+1:end) = [];
%                         scan_num_vol = floor(scan_num_vol);
                     end
                    if cnt_EIT < total_EIT && scan_num_cm ~= header.scan_num
                        msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num_cm) ' scan)']);
                    end
                    
                    Data.EIT.cm_src(:,cnt_scan:cnt_scan+scan_num-1) = ...
                        reshape((raw_data(2:header.cm_size:end).^2 + raw_data(3:header.cm_size:end).^2).^0.5,16,scan_num) *  header.cm_d2i;
                    Data.EIT.cm_snk(:,cnt_scan:cnt_scan+scan_num-1) = ...
                        reshape((raw_data(4:header.cm_size:end).^2 + raw_data(5:header.cm_size:end).^2).^0.5,16,scan_num) *  header.cm_d2i;
                    Data.EIT.cm_mean(:,cnt_scan:cnt_scan+scan_num-1) = (Data.EIT.cm_src(:,cnt_scan:cnt_scan+scan_num-1) + Data.EIT.cm_snk(:,cnt_scan:cnt_scan+scan_num-1)) * 0.5;
                    clear raw_data;
                end
            end
            
            % voltage
            if 1
                if cnt.voltage >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).voltage);
                    raw_data = fread(fid,'double');
                    fclose(fid);
                    scan_num_vol = size(raw_data,1)/header.meas_num;
                    if mod(size(raw_data,1),header.meas_num) ~= 0 || scan_num_vol ~= scan_num
                        if cnt_EIT < total_EIT
                            disp('voltage data miss');
                        end
                        if scan_num_vol > scan_num
                            disp(['voltage data larger than raw(' num2str(scan_num_vol-scan_num) 'scans)']);
                            scan_num_vol = scan_num;
                     raw_data(header.meas_num*scan_num+1:end) = [];
                        end
%                         raw_data(floor(scan_num_vol)*header.meas_num+1:end) = [];
%                         scan_num_vol = floor(scan_num_vol);
                        raw_data(scan_num_vol*header.meas_num+1:scan_num*header.meas_num) = NaN;
                     end
                    if cnt_EIT < total_EIT && scan_num_vol ~= header.scan_num
                        msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num_vol) ' scan)']);
                    end
                    Data.EIT.volt(:,cnt_scan:cnt_scan+scan_num-1) = reshape(raw_data,header.meas_num,scan_num);
                    clear raw_data;
                end
            end
            
            % phase
            if 0
                if cnt.voltage >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).phase);
                    raw_data = fread(fid,'double');
                    fclose(fid);
                    scan_num_phase = size(raw_data,1)/header.meas_num;
                    if mod(size(raw_data,1),header.meas_num) ~= 0 || scan_num_phase ~= scan_num
                        if cnt_EIT < total_EIT
                            disp('phase data miss');
                        end
                        if scan_num_phase > scan_num
                            scan_num_phase = scan_num;
                            disp(['voltage data larger than raw(' num2str(scan_num_phase-scan_num) 'scans']);
                        end
%                         raw_data(floor(scan_num_phase)*header.meas_num+1:end) = [];
%                         scan_num_phase = floor(scan_num_phase);
                        raw_data(scan_num_phase*header.meas_num+1:scan_num*header.meas_num) = NaN;
                    end
                    if cnt_EIT < total_EIT && scan_num_phase ~= header.scan_num
                            msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num_phase) ' scan)']);
                    end
                    Data.EIT.phase(:,cnt_scan:cnt_scan+scan_num_phase-1) = reshape(raw_data,header.meas_num,scan_num_phase);
                    clear raw_data;
                end
            end
            
            % ecg
            if 1
                if cnt.ecg >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).ecg);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    
                    if mod(length(raw_data),37) ~= 0
                        if cnt_EIT < total_EIT
                            disp('ECG data miss');
                        end
                        raw_data(end-mod(length(raw_data),37)+1:end) = [];
                    end
                    raw_data = reshape(raw_data,37,[]);
                    Data.ECG.idx_ecg(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                    Data.ECG.ECG_fs(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = raw_data(header.ecg_size,:) * header.fs_EIT;
                    Data.ECG.ECG_hv(:,(cnt_scan-1)*Data.ECG.ECG_fs/header.fs_EIT+1:(cnt_scan+size(raw_data,2)-1)*Data.ECG.ECG_fs/header.fs_EIT) = ...
                    double(typecast(reshape(uint8(raw_data(header.ecg_data:header.ecg_data+2*Data.ECG.ECG_fs(cnt_scan)/header.fs_EIT-1,:)),1,[]),'uint16'));
                    clear raw_data;
                end
            end
            cnt_scan = cnt_scan + scan_num;
        end
        Data.EIT.t_hms = datetime(epochtime,'ConvertFrom','epochtime','TicksPerSecond',1e3);
        
        % Calc CI
        if sum(Data.EIT.CI_flag == 1) > 0
            [~, CI_idx]=find(Data.EIT.CI_flag==1);
            Data.EIT.CI=Data.EIT.volt(header.CI_ch,CI_idx)./Data.EIT.Amp(CI_idx); % contact impedance
        end
        
        % volt to ohm
        if isfield(Data.EIT, 'volt')
            Data.EIT.ohm = Data.EIT.volt./repmat(Data.EIT.Amp,256,1);
            Data.EIT.ohm(isinf(Data.EIT.ohm)) = NaN;
        end
        
%% v0101 for AirTom/HemoVista
    case 5 % 0101      
        % raw data
        header.meas_num = 256; % 256 measurement
        header.header_num = 24; % 24 byte header of 1 data
        header.fs_EIT = 100;
        header.scan_num = 360000;
        header.protocolkey1 = 5;
        header.protocolkey2 = 6;
        header.ci_flag = 7;
        header.as_flag = 8;
        header.idx_scan = 9:12;
        header.timestamp = 13+1; % bug fixed(x) 13 + 1
        header.bValid = 23;
        header.bSQI = 24;
        
        fid = fopen(Path(1).eit);
        data_version = fread(fid,10000,'uint8');
        data_version = typecast(uint8(data_version),'int32');
        fclose(fid);
        
        % check data format
        if sum(diff(data_version(1:2328/4:end))) == 0 % Impedance (4 byte for 1 data)
            opt.data_ver = 1;
            header.data_num = 2328 - header.header_num;
        elseif sum(diff(data_version(1:2584/4:end))) == 0 % w/ gain
            opt.data_ver = 2;
            header.data_num = 2584 - header.header_num;
        elseif sum(diff(data_version(1:1536/4:end))) == 0 % pre
            opt.data_ver = 0;
            header.data_num = 1536 - header.header_num;
        end
        
        % setting file
        header.set_scale_volt = 0.0000221;
        header.set_total_num = 2632;
        header.set_size_interch = 2608;
        header.set_Amp = 1:8;
        header.set_Gain_table_CodeKey = 9;
        header.set_ScanIdx = 2618:2621;
        header.set_TimeStamp = (2622:2629)+1; % idx error +1
        header.set_GainX = 1:128;
        header.set_GainDigi = 129:144;
        header.set_GainDigi2 = 145:160;
        header.set_Protocoltype = 161;
        header.set_ProjNum = 162;
        header.set_Update = 163;
        
        % ecg file
        header.ecg_size = 69;
        header.ecg_idx = 2;
        header.ecg_data = 6;

        % pressure file
        header.pres_size = 16;
%         header.pre_idx = 2;
%         header.pre_data = 6;

        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            if getappdata(f,'canceling')
                break;
            end
            
            % EIT raw
            if 1
            waitbar(cnt_EIT/total_EIT, f, ['Loading EIT bin data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)]);

            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            
            if rem(size(raw_data,1)/(header.header_num + header.data_num),10000) == 0
                header.scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            end
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            if mod(size(raw_data,1),(header.header_num + header.data_num)) ~= 0
                raw_data(floor(scan_num)*(header.header_num + header.data_num)+1:end) = [];
                scan_num = floor(scan_num);
            end
            if cnt_EIT < total_EIT && scan_num ~= header.scan_num
                msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num) ' scan)']);
            end
            raw_data = reshape(raw_data,(header.header_num + header.data_num),scan_num);

%             header_data=raw_data(1:header.header_num,:);
            
            Data.EIT.ProtocolKey1(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.protocolkey1,:);
            Data.EIT.ProtocolKey2(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.protocolkey2,:);
            Data.EIT.CI_flag(cnt_scan:cnt_scan+scan_num-1) = logical(raw_data(header.ci_flag,:));
            Data.EIT.AS_flag(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.as_flag,:);   
            Data.EIT.idx_scan(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.idx_scan(4),:)*256^3 + raw_data(header.idx_scan(3),:)*256^2 + raw_data(header.idx_scan(2),:)*256^1 + raw_data(header.idx_scan(1),:);
            Data.EIT.Valid(cnt_scan:cnt_scan+scan_num-1) = logical(raw_data(header.bValid,:));
            Data.EIT.SQI(cnt_scan:cnt_scan+scan_num-1) = logical(raw_data(header.bSQI,:));
                   
            % find time data point & block last byte
            temp = sum((diff(raw_data(header.timestamp-1:header.timestamp+7,:)'))>0);
            [~,tp] = max(temp);
            header.timestamp = header.timestamp-2+tp;
            raw_data(header.timestamp+7,:) = 0;
            
            epochtime(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(raw_data(header.timestamp:header.timestamp+7,:)),1,[]),'uint64');
%             header.timestamp = 16;
%             datetime(typecast(reshape(uint8(header_data(header.timestamp:header.timestamp+7,:)),1,[]),'uint64'),'ConvertFrom','epochtime','TicksPerSecond',1e3);
            raw_data(1:header.header_num,:) = [];
%             datetime(epochtime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3)
            end
            
            % raw_data
            if 1
                raw_data = raw_data(:);
                switch opt.data_ver
                    case 1
                        raw_data = reshape(raw_data,9,[]);
                        Data.EIT.OV_flag(:,cnt_scan:cnt_scan+scan_num-1) = reshape(logical(raw_data(1,:)),header.meas_num,[]);
                        Data.EIT.phase = reshape(typecast(reshape(uint8(raw_data(2:5,:)),1,[]),'single'),header.meas_num,[]);
                        Data.EIT.ohm = reshape(typecast(reshape(uint8(raw_data(6:9,:)),1,[]),'single'),header.meas_num,[]);
                    case 0
                        Data.EIT.OV_flag(:,cnt_scan:cnt_scan+scan_num-1) = reshape(raw_data(2:6:end),header.meas_num,length(raw_data)/header.meas_num/6);

                        V_R = raw_data(3:6:end) + raw_data(4:6:end).*256;
                        V_R(V_R>32767) = V_R(V_R>32767) - 65536;
                        V_I = raw_data(5:6:end) + raw_data(6:6:end).*256;
                        V_I(V_I>32767) = V_I(V_I>32767) - 65536;
                        raw_data = V_R + V_I*1i;
                        raw_data = reshape(raw_data,header.meas_num,length(raw_data)/header.meas_num);
                        switch opt.type
                            case 'magnitude'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data);
                            case 'complex'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = raw_data;
                            case 'phase'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data);
                                Data.EIT.phase(:,cnt_scan:cnt_scan+scan_num-1) = angle(raw_data);
                                if opt.mask == 1
                                    Data.EIT.phase(Data.ch_mask,:) = [];
                                end
                        end
                end
                clear raw_data;
            end
            
            % ecg
            if isfield(Path,'ecg')
                if cnt.ecg >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).ecg);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    
                    if cnt_EIT == 1
                        header.ecg_size = round(length(raw_data)/length(epochtime));
%                         header.ECG_fs = (header.ecg_size-1-4)/4; % 1:num, 4:idx
                    end
                
                    if mod(length(raw_data),header.ecg_size) ~= 0
                        if cnt_EIT < total_EIT
                            disp('ECG data miss');
                        end
                        raw_data(end-mod(length(raw_data),37)+1:end) = [];
                    end
                    raw_data = reshape(raw_data,header.ecg_size,[]);
                    Data.ECG.idx_ecg(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                    Data.ECG.ECG_fs = header.fs_EIT;
                    Data.ECG.ECG_hv(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = sum(raw_data(header.ecg_data:header.ecg_data+3,:).*[1 256 256^2 256^3]');
%                     Data.ECG.ECG_hv(:,(cnt_scan-1)*Data.ECG.ECG_fs/header.fs_EIT+1:(cnt_scan+size(raw_data,2)-1)*Data.ECG.ECG_fs/header.fs_EIT) = ...
%                         double(typecast(reshape(uint8(raw_data(header.ecg_data:end,:)),1,[]),'uint16'));
                    clear raw_data;
                end
            end
        
            % pressure
            if isfield(Path,'pressure')
                if cnt.pressure >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).pressure);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    
                    if mod(length(raw_data),header.pres_size) ~= 0
                        if cnt_EIT < total_EIT
                            disp('Pressure data miss');
                        end
                        raw_data(end-mod(length(raw_data),header.pres_size)+1:end) = [];
                    end
                    raw_data = reshape(raw_data,header.pres_size,[]);
                    
                    header.pres_idxscan = 13:16;
                    Data.Pressure.idx_press(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.pres_idxscan,:)),1,[]),'uint32'));
                    
                    temp_flow = raw_data(3,:) + raw_data(4,:).*256;
                    temp_flow(temp_flow>32767) = temp_flow(temp_flow>32767) - 65536;
                    Data.Pressure.Paw(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = temp_flow;
                    
                    temp_flow = raw_data(9,:) + raw_data(10,:).*256;
                    temp_flow(temp_flow>32767) = temp_flow(temp_flow>32767) - 65536;
                    Data.Pressure.Faw(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = temp_flow;
                    clear raw_data;
                end
            end
            
            cnt_scan = cnt_scan + scan_num;
        end
        
        % idx checking (start number is 1 for matlab index)
        if Data.EIT.idx_scan(1) == 0
            Data.EIT.idx_scan = Data.EIT.idx_scan + 1;
        end
        
%         if sum(diff(Data.EIT.idx_scan) ~= 1)
%             if  find(diff(Data.EIT.idx_scan) ~= 1) == 1 || Data.EIT.idx_scan(2) == 1
%                 Data.EIT.idx_scan(1) = Data.EIT.idx_scan(2)-1;
%             else
%                 msgbox('scan index missing');
%             end
%         end
        
        if isfield(Path,'ecg')
            if Data.ECG.idx_ecg(1) == 0
                Data.ECG.idx_ecg = Data.ECG.idx_ecg + 1;
            end
            
            if sum(diff(Data.ECG.idx_ecg) ~= 1)
%                 if  find(diff(Data.ECG.idx_ecg) ~= 1) == 1 || Data.ECG.idx_ecg(2) == 1
%                     Data.ECG.idx_ecg(1) = Data.ECG.idx_ecg(2)-1;
%                 else
%                     msgbox('scan index missing');
%                 end
            end
        end
        
        if isfield(Path,'pressure')
            if Data.Pressure.idx_press(1) == 0
                Data.Pressure.idx_press = Data.Pressure.idx_press + 1;
            end
            
%             if sum(diff(Data.Pressure.idx_press) ~= 1)
%                 if  find(diff(Data.Pressure.idx_press) ~= 1) == 1 || Data.Pressure.idx_press(2) == 1
%                     Data.Pressure.idx_press(1) = Data.Pressure.idx_press(2)-1;
%                 else
%                     msgbox('scan index missing');
%                 end
%             end
        end
        
        % setting file
        if isfield(Path,'setting')
            fid = fopen(Path(1).setting);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            
            if rem(size(raw_data,1),(header.set_total_num)) ~= 0
                if rem(size(raw_data,1),(header.set_total_num)+8) ~= 0
                    msgbox('data missing: Setting_bin');
                    header.set_loadfail = 1;
                else
                    header.set_total_num = header.set_total_num + 8;
                    header.set_nAS = size(raw_data,1)/header.set_total_num;
                    header.set_loadfail = 0;
                    header.set_eleconfig = 1;
                end
            else
                header.set_nAS = size(raw_data,1)/header.set_total_num;
                header.set_loadfail = 0;
                header.set_eleconfig = 0;
            end
            
            if header.set_loadfail ~= 1
                raw_data = reshape(raw_data,header.set_total_num,header.set_nAS);
                
                if header.set_eleconfig == 1
                    Data.EIT.Settings.Elec_config = raw_data(1,1);
                    raw_data([1:4 end-3:end],:) = [];
                end
                
                Data.EIT.Settings.Amp = typecast(reshape(uint8(raw_data(header.set_Amp,:)),1,[]),'double');
                Data.EIT.Settings.table_CodeKey = raw_data(header.set_Gain_table_CodeKey,:);
                Data.EIT.Settings.idx_scan = raw_data(header.set_ScanIdx(4),:)*256^3 + raw_data(header.set_ScanIdx(3),:)*256^2 + raw_data(header.set_ScanIdx(2),:)*256^1 + raw_data(header.set_ScanIdx(1),:);
                
                try
                    Data.EIT.Settings.t_hms = datetime(typecast(reshape(uint8(raw_data(header.set_TimeStamp+1,:)),1,[]),'uint64'),'ConvertFrom','epochtime','TicksPerSecond',1e3);
                catch
                    try % error 'TicksPerSecond' option (matlab version)
                        Data.EIT.Settings.t_hms = datetime(typecast(reshape(uint8(raw_data(header.set_TimeStamp+1,:)),1,[]),'uint64')/1e3,'ConvertFrom','epochtime');
                    catch
                        msgbox('Datetime data is corrupted');
                    end
                end
                
                raw_data = raw_data([10:10+header.set_size_interch-1],:);
                for cnt_AS = 1:header.set_nAS
                    gain_data = reshape(raw_data(:,cnt_AS),[],16);
                    Data.EIT.Settings.ProtocolType(:,cnt_AS) = gain_data(header.set_Protocoltype,:);
                    Data.EIT.Settings.ProjNum(:,cnt_AS) = gain_data(header.set_ProjNum,:);
                    Data.EIT.Settings.Update(:,cnt_AS) = gain_data(header.set_Update,:);
                    
                    Data.EIT.Settings.GainX(:,cnt_AS) = typecast(uint8(reshape(gain_data(header.set_GainX,:),[],1)),'double');
                    Data.EIT.Settings.GainDigi1(:,cnt_AS) = reshape(gain_data(header.set_GainDigi,:),[],1);
                    Data.EIT.Settings.GainDigi2(:,cnt_AS) = reshape(gain_data(header.set_GainDigi2,:),[],1);
                end
                
                % EIT raw -> ohm
                if opt.data_ver == 0
                    Data.EIT.raw(:,~logical(Data.EIT.Valid)) = NaN;
                    scaler_ohm = zeros(size(Data.EIT.raw));
                    for cnt_AS = 1:header.set_nAS
                        tp_s = Data.EIT.Settings.idx_scan(cnt_AS)+1;
                        if cnt_AS ~= header.set_nAS % check next AS exist!
                            tp_e = Data.EIT.Settings.idx_scan(cnt_AS+1);
                        else % apply whole data
                            tp_e = size(Data.EIT.raw,2);
                        end
                        scaler_ohm(:,tp_s:tp_e) = repmat(header.set_scale_volt./Data.EIT.Settings.GainX(:,cnt_AS)/Data.EIT.Settings.Amp(cnt_AS),1,tp_e-tp_s+1);
                    end
                    Data.EIT.ohm = Data.EIT.raw .* scaler_ohm;
                end
                
                %             Data.EIT = rmfield(Data.EIT,'raw');
                if opt.mask == 1
                    Data.EIT.ohm(Data.ch_mask,:) = [];
                end
            end
        end
        
        try
            Data.EIT.t_hms = datetime(epochtime,'ConvertFrom','epochtime','TicksPerSecond',1e3);
        catch
            try % error 'TicksPerSecond' option (matlab version)
                Data.EIT.t_hms = datetime(epochtime/1e3,'ConvertFrom','epochtime');
            catch
                msgbox('Datetime data is corrupted');
            end
        end
        
        if isfield(Path,'pressure')
            % pressure calibration
            Data.Pressure.Paw = Data.Pressure.Paw*0.26702779+0.13351485; % Gage->Pa
            Data.Pressure.Paw = Data.Pressure.Paw*0.0101972; % Pa->cmH2O
            
            idx_inv = Data.Pressure.Faw>0;
            Data.Pressure.Faw = (Data.Pressure.Faw*0.00009719); % Pa->cmH2O
            Data.Pressure.Faw = -0.941*(pow2(Data.Pressure.Faw))+(52.38*abs(Data.Pressure.Faw))+1; % cmH2O->L/m
            Data.Pressure.Faw(idx_inv) = -Data.Pressure.Faw(idx_inv);
            
            if Data.Pressure.idx_press(1) == 1 % first data error
                Data.Pressure.idx_press = [0 Data.Pressure.idx_press];
                Data.Pressure.Paw = [0 Data.Pressure.Paw];
                Data.Pressure.Faw = [0 Data.Pressure.Faw];
            end
            
            if length(Data.Pressure.idx_press) ~= length(epochtime) % data number missmatch
                if length(Data.Pressure.idx_press) < length(epochtime)
                    Data.Pressure.idx_press(length(Data.Pressure.idx_press)+1:length(epochtime)) = 0;
                    Data.Pressure.Paw(length(Data.Pressure.idx_press)+1:length(epochtime)) = 0;
                    Data.Pressure.Faw(length(Data.Pressure.idx_press)+1:length(epochtime)) = 0;
                else
                    Data.Pressure.idx_press(length(epochtime)+1:end) = [];
                    Data.Pressure.Paw(length(epochtime)+1:end) = [];
                    Data.Pressure.Faw(length(epochtime)+1:end) = [];
                end
            end
        end
        
        if sum(diff(epochtime)<=0) > 1
            rm_idx = [false diff(epochtime)<=0];
            FN = fieldnames(Data.EIT);
            for cnt_field = 1:length(FN)
                if ~isstruct(Data.EIT.(FN{cnt_field}))
                    Data.EIT.(FN{cnt_field})(:,rm_idx) = [];
                end
            end
            if isfield(Path,'pressure')
                FN = fieldnames(Data.Pressure);
                for cnt_field = 1:length(FN)
                    if ~isstruct(Data.Pressure.(FN{cnt_field}))
                        Data.Pressure.(FN{cnt_field})(rm_idx) = [];
                    end
                end
            end
        end

%% v1000 for AirTom-RM
    case 8 % 1000
        header.meas_num = 256; % 256 measurement
        header.header_num = 72; % 72 byte header of 1 data
        header.data_num = 1536;
        
        header.amp = 5;
        header.gain = 13;
        header.vsf = 45;
        header.ci_flag = 53;
        header.as_flag = 54;
        header.time = 55;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            fclose(fid);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            
            for cnt_scan2 = 1:scan_num
                
                if(mod(cnt_scan2, 1000) == 0)
                    if getappdata(f,'canceling')
                        break;
                    end
                    waitbar(cnt_scan2/scan_num, f, ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT) ' (' num2str(cnt_scan2) ' of ' num2str(scan_num) ')']);
                end
                
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM   L H     
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM   L H     
                temp(5:6,:) = [];
                temp(temp>32767) = temp(temp>32767) - 65536;
                temp = temp';
                Data.EIT_raw(:,cnt_scan) = sqrt(temp(:,3).^2 + temp(:,4).^2);
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:4
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                
                Data.VSF(cnt_scan)=typecast(uint8(temp2(header.vsf:header.vsf+7)),'double');
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                hms.Y(cnt_scan) = double(typecast(uint8(temp2(header.time:header.time+1)),'uint16'));
                hms.M(cnt_scan) = temp2(header.time+2);
                hms.D(cnt_scan) = temp2(header.time+3);
                hms.H(cnt_scan) = temp2(header.time+4);
                hms.MI(cnt_scan) = temp2(header.time+5);
                hms.S(cnt_scan) = temp2(header.time+6);
                hms.MS(cnt_scan) = double(typecast(uint8(temp2(header.time+7:header.time+8)),'uint16'));
                
                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).voltage);
                raw_data = fread(fid,'double');
                fclose(fid);
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
            % phase
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).phase);
                raw_data = fread(fid,'double');
                fclose(fid);
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                end
            end
        end
        Data.t_hms = datetime(hms.Y, hms.M, hms.D, ...
            hms.H, hms.MI, hms.S, hms.MS);       
end

    delete(f);
end

function [mask] = FxEIT_mask(ch,first_sat)
    temp = zeros(ch,ch);
    if nargin == 1
        first_sat = [ch,1,2];
    end

    if length(first_sat) == 4
        cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3); cnt4 = first_sat(4);
        for i = 1:ch
            temp([cnt1,cnt2,cnt3,cnt4],i) = 1;
            cnt1 = cnt1 + 1;
            cnt2 = cnt2 + 1;
            cnt3 = cnt3 + 1;
            cnt4 = cnt4 + 1;
            if cnt1 > ch
                cnt1 = 1;
            end
            if cnt2 > ch
                cnt2 = 1;
            end
            if cnt3 > ch
                cnt3 = 1;
            end
            if cnt4 > ch
                cnt4 = 1;
            end
        end
    else
        cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3);
        for i = 1:ch
            temp([cnt1,cnt2,cnt3],i) = 1;
            cnt1 = cnt1 + 1;
            cnt2 = cnt2 + 1;
            cnt3 = cnt3 + 1;
            if cnt1 > ch
                cnt1 = 1;
            end
            if cnt2 > ch
                cnt2 = 1;
            end
            if cnt3 > ch
                cnt3 = 1;
            end
        end
    end

    temp2 = reshape(temp,ch*ch,1);
    [mask, ~] = find(temp2 == 1);
end