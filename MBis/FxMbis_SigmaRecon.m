function output = FxMbis_SigmaRecon(Ref_path,Data_path)

% import RawSmat & Simdata
func_path = which('FxMbis_SigmaRecon');
func_path(end-18:end) = [];
load(strcat(func_path,'RawSmat_Sigma_Recon'));
load(strcat(func_path,'Simdata_Sigma_Recon'));

Z = FxMbis_Dataimport(Ref_path);
for cnt1 = 1:size(Z,2)
    calibrationfactor(:,cnt1) = Simdata./Z(:,cnt1);
end
clear Z;

Z = FxMbis_Dataimport(Data_path);
for cnt1 = 1:size(Z,2)
    Calibrateddata(:,cnt1) = Z(:,cnt1).*calibrationfactor(:,cnt1);
end

for cnt = 1:size(Calibrateddata,2)
    TargetPosition=14;  % We only interested in 14th position in 18 voxels
    NewSmat=zeros(4,size(RawSmat,2));
    NewData=zeros(4,1);
    RawData = Calibrateddata(:,cnt);
    for ComputeBlock = 1:4
        % Load block
        SmatB=RawSmat(ComputeBlock:4:6*4+ComputeBlock,:);
        DataB=RawData(ComputeBlock:4:6*4+ComputeBlock,1);
        % Coefficient xi
        Cvec = SmatB(:,TargetPosition);
        normV= norm(Cvec,2);
        alpha= normV;  % Parameter of minimization of method
        Amat = SmatB; Amat(:,TargetPosition)=0;
        Mmat = (Amat*Amat'+alpha*eye(size(SmatB,1)));
        Coef = alpha*Mmat\Cvec;
        % LinCombSensMat
        NewSmat(ComputeBlock,1:size(SmatB,2))=Coef'*SmatB;
        NewData(ComputeBlock,1)=Coef'*DataB;
    end
    result(:,cnt) = abs(pinv(NewSmat)*NewData);
end
output = 1./result(14,:)';