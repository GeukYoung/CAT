function Z = FxMbis_Dataimport(path)
    
freqnum = [100, 1, 5, 10, 50, 100];

for cnt0 = 1:length(freqnum)
    freq_str = int2str(freqnum(1,cnt0));
    filepath = strcat(path,'\',freq_str,'kHz');

    a1=strcat(filepath,'\Scan1.txt');
    a2=strcat(filepath,'\Source.txt');
    file1 = textread(a1);
    file2 = textread(a2);

    datamag(:,cnt0) = sqrt(file1(:,3).^2+file1(:,4).^2);
    sourcemag(:,cnt0) = sqrt(file2(:,3).^2+file2(:,4).^2);
    datareal(:,cnt0) = file1(:,3);
    dataquad(:,cnt0) = file1(:,4);
    sourcereal(:,cnt0) = file2(:,3);
    sourcequad(:,cnt0) = file2(:,4);
 
end

i = sqrt(-1);
Rs = 260;

datareal=datareal';
dataquad=dataquad';
sourcereal=sourcereal';
sourcequad=sourcequad';

for cnt1 = 1:28
    data{cnt1,1}(:,1) = datareal(:,cnt1)+dataquad(:,cnt1)*i;
    source{cnt1,1}(:,1) = sourcereal(:,cnt1)+sourcequad(:,cnt1)*i;
    RL(:,cnt1) = Rs*data{cnt1,1}(:,1)./source{cnt1,1}(:,1);
end

Z = RL'