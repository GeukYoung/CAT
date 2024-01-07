function output = FxAnisoRecon(Ref_path,Data_path)

% import RawSmat & Simdata
func_path = which('FxMbis_SigmaRecon');
func_path(end-18:end) = [];
load(strcat(func_path,'RawSmat_Sigma_Recon'));
load(strcat(func_path,'RawSmat_aniso'));

data_num = 6;
degree = '0';

Ref_path = strcat(Ref_path,'\');
%% find califactor

RL=FxMbis_Dataimport(Ref_path);
RL=RL';
for i = 1:8
    califactor(:,i) = sim0./abs(RL(:,i));
end

%% Recon
folder_path=strcat(Data_path,'\');
RL=FxMbis_Dataimport(folder_path);
RL=RL';
for i = 1:8
    calibrated(:,i) = RL(:,i).*califactor(:,i);
end

for cnt = 1:8
    InvSigmaComp = abs(RawSmat)\abs(calibrated(:,cnt));
    
    InvSigma22=[
        InvSigmaComp(1) InvSigmaComp(4)  ;
        InvSigmaComp(4) InvSigmaComp(2)  ];
    ReconSigma22=inv(InvSigma22);
    sigma = real(ReconSigma22);
    allsigma{cnt,1} = sigma;
    [eigvector,eigvalue]  = eig(sigma);
    alleigvector{cnt,1} = eigvector;
    alleigvalue{cnt,1} = eigvalue;
    ratio = abs(eigvalue(1,1))/abs(eigvalue(2,2));
    allratio(cnt,1) = ratio;
end
output = allratio;
