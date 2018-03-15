% % --------------------------------------------------------------
% % 27.Mar.2017 (c) Arash Mirhashemi
% % An example code for:
% %     loading a SpecTex HyperSpectral Image Cube --------------------- using t_tif_read.m
% %     simulating the RGB image of the spectral cube ------------------ using t_spd2sth.m
% %     calculating some of its moment features and visualize them ----- using t_mom_ftr
% %     saving the SpecTex HyperSpectral Image Cube -------------------- using t_tif_write.m
% %     loading and locating SpecTex sample spectra -------------------- with  SpecTex_spectra.mat
% % --------------------------------------------------------------

%% --------------------------------------------------------------

% loading the SpecTex sample T01 image cube
cube = t_tif_read('T01.tif');

% loading a sub-number of bands, plus the tif-stack structure of the cube as well
[im, tumb, ifd, tif] = t_tif_read('T01.tif','bands',[10, 5, 6, 39]);

% accessing the wavelength information of the image cube from the tif structure
lmbd = tif(1+~isempty(tumb)).WaveLength;

% accessing any ASCII meta data, if it exist in the file 
% here, since no meta data is saved with the file an empty string '' is returned 
meta = char(ifd(1+~isempty(tumb)).value{ifd(1+~isempty(tumb)).tag==65111});

% checking the tumbnail image that was saved with the sample (simulated RGB under D65)
% and one of the bands as a greyscale image
figure('Units','Normalized','Position',[0 0 1 1]);
subplot(1,3,1); imshow(tumb); title('RGB-D65');
subplot(1,3,2); imshow(im(:,:,1)); title('band-#10');

% simulating the RGB tristimulus values under equal-energy illumination
RGB_EE = t_spd2sth(cube,lmbd,'EE',1964,'RGB');
subplot(1,3,3); imshow(RGB_EE); title('RGB-EE');

% calculating the moment feature for 
MoR1 = t_mom_ftr(cube, 'MoR1');

% visualizing the feature image and its histogram
figure('Units','Normalized','Position',[0 0 1 1]);
subplot(1,2,1); imshow( uint8((2^8-1)*mat2gray(MoR1)) ); title('MoR1 feature image');
subplot(1,2,2); hist(MoR1(:),1000); title('MoR1 histogram of values');

% calculating several moment features, with the original wavelength values (and not the default Natural Numbers)
mom_cube = t_mom_ftr(cube, {'Mir3','Sor5','MoR2'},'lmbd',lmbd);

% saving the spectral cube to a new tif-stack file with uint8 datatype, and with the RGB_EE image as the thumbnail image
t_tif_write(uint8(cube*(2^8-1)),'T01-uint8-EE.tif','tumb',uint8(RGB_EE*(2^8-1)),'datatype','uint8');
% saving the spectral cube to a new tif-stack file with uint16 datatype, and with an RGB made of 
% triangular sensitivity functions as the thumbnail image
t_tif_write(uint16(cube*(2^16-1)),'T01-uint16-triangular.tif','tumb','triangular','datatype','uint16');

% clearing up 
clear cube im tumb ifd tif lmbd meta RGB_EE MoR1 mom_cube 

%% --------------------------------------------------------------

% loading the .mat file, containing the sample spectra from SpecTex database
load('SpecTex_spectra.mat');
% The file contains a structure called SpecTex
% SpecTex.selection contains the sample spectra from each of the classes 
% to access one sample use this syntax
class_number = 1; % that is T01
example = SpecTex.selection{class_number}; % there are 1366 sample spectra in T01
plot(example); % this plots all the spectra 
% SpecTex.selec_ind contains the index of the pixel location where the spectra in SpecTex.selection are selected from
% the first number in SpecTex.selec_ind is the "linear index" of the pixel 
% the next two numbers in SpecTex.selec_ind are the (x,y) subscripts of the pixel
% for example SpecTex.selection{1}(:,700) is coming from T01-pixel-(579,325) 
% therefore SpecTex.selec_ind{1}(700,2:3) is equal to (579,325)

% % Random Grid Lab space sampling:
% The SpecTex.selection samples are collected by calculating the Lab values of all pixels in each SpecTex image cube,
% and sampling the Lab space with a random grid with uniform distribution.
% Furthermore, another collection is provided in SpecTex.spectra by sampling the SpecTex.selection spectra for a 2nd time
% This time, all the spectra in SpecTex.selection are collected (that is a total of 178684 spectra), 
% then their Lab space is sampled with another random grid with unifirm distribution.
% This sampling provided a set of 10k spectra that is saved as SpecTex.selection
% The Lab values of the SpecTex.selection spectra is also provided in SpecTex.gamut

%% --------------------------------------------------------------













