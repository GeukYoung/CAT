function [Data] = FxImport_EIT_v5_0(varargin)

narginchk(1, inf);

if ischar(varargin{1})
    Path_main = varargin{1};
    varargin = varargin(2:end);
%     argOffset = 1;
else
    error('ERROR: Data path');
end

op.block = 0;
% op.type = 'mag';

numOrigInputArgs = numel(varargin);
if numOrigInputArgs ~= 0
    for cnt_Args = 1:2:numOrigInputArgs
        switch varargin{cnt_Args}
            case 'block'
                if strcmpi(varargin{cnt_Args+1},'first')
                    op.block = 1;
                elseif strcmpi(varargin{cnt_Args+1},'last')
                    op.block = 99;
                elseif isnumeric(varargin{cnt_Args+1})
                    op.block = varargin{cnt_Args+1}+1;
                end
        end
    end
end

% isCreatedUIfigure = false;

% if(nargin<2)
%     fig = uifigure;
%     isCreatedUIfigure = true;
% end

if(1) % if(isempty(fig))
    fig = uifigure;
    isCreatedUIfigure = true;
end

dlg = uiprogressdlg(fig,'Title','Please Wait',...
                    'Message','Loading EIT data','Indeterminate','on', 'Cancelable','on');
drawnow;

% cd(Path_main);
% dirlist = dir('.');

dirlist = dir(Path_main);

cnt.eit = 0;
cnt.voltage = 0;
cnt.phase = 0;
cnt.ecg = 0;
cnt.cm = 0;

for cnt_file = 1:length(dirlist)
    temp_string = strsplit(dirlist(cnt_file).name,{'.', '['});
    if contains(temp_string{1},'EIT_','IgnoreCase',true)
        idx.eit(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.eit = cnt.eit + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).eit = fullfile(Path_main, dirlist(cnt_file).name);
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

if op.block ~= 0
    if size(op.block) == 1
        switch op.block
            case 1
                Path = Path(1);
                total_EIT = 1;
            case 99
                Path = Path(end);
                total_EIT = 1;
        end
    elseif length(op.block) > 1         
        try
            Path = Path(op.block);
            total_EIT = length(op.block);
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
Data.mask_0skip = FxEIT_mask(16);
% Data.mask_2skip = FxEIT_mask(16,[16 1 3 4]);

set(dlg,'Indeterminate','off');
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
            if dlg.CancelRequested
                Data = [];
                close(dlg)
                if(isCreatedUIfigure)
                    close(fig)
                end
                return
            end
            
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            dlg.Message = ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
            
            for cnt_scan2 = 1:scan_num
                if(mod(cnt_scan2, 1000) == 0)
                    dlg.Value = cnt_scan2/scan_num;
                    dlg.Message = ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT) ' (' num2str(cnt_scan2) ' of ' num2str(scan_num) ')'];
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
            if dlg.CancelRequested
                Data = [];
                close(dlg)
                if(isCreatedUIfigure)
                    close(fig)
                end
                return
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
                    
                    if dlg.CancelRequested
                        Data = [];
                        close(dlg)
                        if(isCreatedUIfigure)
                            close(fig)
                        end
                        return
                    end
                    
                    dlg.Value = cnt_scan2/scan_num;
                    dlg.Message = ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT) ' (' num2str(cnt_scan2) ' of ' num2str(scan_num) ')'];
                end
                
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM�� L H ����
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM�� L H ����
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
            if dlg.CancelRequested
                Data = [];
                close(dlg)
                if(isCreatedUIfigure)
                    close(fig)
                end
                return
            end
            
            % EIT raw
            if 1
            dlg.Message = ['Loading EIT_bin data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
            dlg.Value = cnt_EIT/total_EIT;
            
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
                dlg.Message = ['Loading cm data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
                dlg.Message = ['Loading voltage data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
                dlg.Message = ['Loading phase data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
                dlg.Message = ['Loading ecg data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
%         Data.EIT.t_hms = datetime(epochtime,'ConvertFrom','epochtime','TicksPerSecond',1e3);
        
        
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
            if dlg.CancelRequested
                Data = [];
                close(dlg)
                if(isCreatedUIfigure)
                    close(fig)
                end
                return
            end
            
            % EIT raw
            if 1
            dlg.Message = ['Loading EIT_bin data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
            dlg.Value = cnt_EIT/total_EIT;
            
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
                dlg.Message = ['Loading cm data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
                dlg.Message = ['Loading voltage data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
                dlg.Message = ['Loading phase data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
                dlg.Message = ['Loading ecg data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
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
            
            dlg.Message = ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)];
            
            for cnt_scan2 = 1:scan_num
                
                if(mod(cnt_scan2, 1000) == 0)
                    if dlg.CancelRequested
                        Data = [];
                        close(dlg)
                        if(isCreatedUIfigure)
                            close(fig)
                        end
                        return
                    end
                    
                    dlg.Value = cnt_scan2/scan_num;
                    dlg.Message = ['Loading EIT data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT) ' (' num2str(cnt_scan2) ' of ' num2str(scan_num) ')'];
                end
                
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM�� L H ����
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM�� L H ����
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

    close(dlg)
    if(isCreatedUIfigure)
        close(fig)
    end

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