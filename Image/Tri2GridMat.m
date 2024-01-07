clear all; clc
EIDORS_path = 'D:\Dropbox\#Lecture\7th\20170417_Dataprocessing_TidalVolume&RespirationRate\eidors-v3.8\startup.m';
EITtool_path = 'D:\Dropbox\#Lecture\7th\20170417_Dataprocessing_TidalVolume&RespirationRate\EIT tool\EITtools.m';
run(EIDORS_path);
run(EITtool_path);

% Mesh Info
load('H:\1.Data\2017.10.17_PSG_EJW\matlab\imdl.mat');

% Tri2Grid function
test_data = 1:length(imdl.fwd_model.elems);

elems = imdl.fwd_model.elems;
nodes = imdl.fwd_model.nodes;
nPixel = 64;
temp_image = FxEIT_Tri2Grid(elems, nodes, test_data, nPixel);
imagesc(temp_image);

%% Generate Triange to Grid Convert Matrix
nTelem = length(imdl.fwd_model.elems);
T2G_ConMat = zeros(nPixel^2,nTelem);
parfor i = 1:nTelem
    Impulse = zeros(1,nTelem);
    Impulse(i) = 1;
    Impulse_response = FxEIT_Tri2Grid(elems, nodes, Impulse, nPixel);
    T2G_ConMat(:,i) = reshape(Impulse_response,size(Impulse_response,1)*size(Impulse_response,2),1);
end

%% Test
% test data
test_data = 1:length(imdl.fwd_model.elems);

% Triangle 2 Grid
test_image = FxEIT_Tri2Grid(elems, nodes, test_data, nPixel); % Conv using func
test_image2 = reshape(T2G_ConMat * test_data', nPixel, nPixel); % Conv using ConvMat

% Display
subplot(131);
imagesc(test_image); axis image; caxis on;
subplot(132);
imagesc(test_image2); axis image; axis on;
subplot(133);
imagesc(test_image-test_image2); axis image; axis on;