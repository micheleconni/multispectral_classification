function ifd = t_tif_write(im, filename, varargin)
% % function ifd = t_tif_write(im, filename, varagin)
% % --------------------------------------------------------------
% % 27.Mar.2017 (c) Arash Mirhashemi
% % A tiff writer for spectral image cubes
% % --------------------------------------------------------------
% % Inputs:
% % im         the input image cube structure (in floating point [0,1] or unsigned integer datatypes)
% % filename   the output file address/name
% % Name-------Value------------------------(optional)
% % tumb       the method to make the thumbnail image, default method is 'average' if not provided
% %            'none' or [] or '', no thumbnail image will be availale
% %            ['average'] = a Grey-scale is used by taking the mean of all bands. This is the default behaviour. 
% %            'triangular' = semi-triangular sensitivities are used to make a psuedo-color (for RGB images it equals to the RGB itself)
% %            'flux' = for FluxData camera, sum of the first two CCDs is used 
% %            [page numbers] = 1-or-3 numbers, the page numbers which are used to make the Grey-scale or the rgb thumbnail
% %            [wavelengths] = the CIE-XYZ vales is used as sensitivities for the given wavelength range (the #wavelengths should be equal to size(im,3))
% %            [sensitivities] = (a 3xsiz.bands matrix) of sensitivities, for makeing the thumbnail by a weighted sum of the pages
% %            [an RGB image] a given RGB image will be force-casted to uint8 [0,255] without checking the input range and will be used as tumb
% % datatype   (default = class(im(1)) ) the datatype that the image cube will be saved as
% % lmbd       (default = 1:size(im,3)) a vector containing the wavelength range of the im-cube (this extra data is saved in tag#65000)
% % minmax     (default [intmin(datatype), intmax(datatype)]) the [minimum, maximum] values that will be written in SMinSampleValue,SMaxSampleValue tags
% % resolution [XResolution(1:2),YResolution(1:2),ResolutionUnit(1)] where ResolutionUnit = 1(no absolute unit)/2(inch)/3(cm)
% %            (default = 72/1(dpi)-->[72,1,72,1,2])
% % num_strip  (default = 1) any integer number between 1 and size(im,1) specifies that in how many strips the data is saved 
% % nvxml      the Natural Vision XML data, or any (ASCII) data will be written in the first IFD after the thumbnail as meta-data
% % 
% % Outputs:
% % ifd        the whole image-file-directory(IFD) of the written tif file in a Matlab structure
% % 
% % --------------------------------------------------------------
% % usage example: tif = t_tif_write(imread('lena.png'), 'test.tif', 'datatype', 'uint8', 'tumb', 'triangular')
% % --------------------------------------------------------------
% % 
% % --------------------------------------------------------------
% % Note-1 Image cube order: This code assumes the image cube is given in a 3D Matlab matrix where 1st-dim is vertical
% %        spatial dimension, 2nd-dim is horizontal spatial dimension, and 3rd-dim is spectral dimension in INCREASING wavelengths
% %        EXCEPTION: for an RGB image we assume 1st page is R, 2nd is G, and 3rd is B because that's the norm to present them
% % Note-2 If lmbd is not provided numbers 1:size(im,3) will be saved as the wavelength range of the im-cube
% % Note-3 The main image cube (im) is not processed (scale,etc), but only force-type-casted and is saved directly as requested
% %        by 'datatype'. It is the users responsibility to take care that the datatype and values of input 'im' matches
% %        for example if one saves a double [0,1] 'im' cube as int8, data will be lost! one should instead inputs im = uint8(255*(im));
% %        also if one saves a uint16 im with uint8 the bigger than 255 values will be saturated, and so on
% %        'datatype' can be all these matlab numeric data-types: uint8-uint16-uint32-uint64-int8-int16-int32-int64-single-double
% % Note-4 When you are providing tumb input as [page numbers], 1st number is used for R, 2nd for G, and 3rd for B
% %        When you are providing tumb input as [wavelengths] it should be in INCREASING order
% %        When you are providing tumb input as [sensitivities], always 1st row is for R, 2nd row for G, and 3rd for B
% %        plus, in this case if the given [wavelengths] have no overlap with visual wavelength range, the image will be black obviously
% % Note-5 If the image cube's 3rd dimension size is 1-or-2, no matter what tumb option you give, no tumb will be saved and only the 1-or-2 pages are saved
% % Note-5 X,Y Resolutions are given as rationales, that is numerator/denominator ex. 72/1dpi-->[72,1,72,1,2]
% % Note-6 The image cube is saved in config = Band-Interleaved by-Pixel(BIP). BSQ is also possible by using (PlanarConfiguration = 2)
% % 
%% ---------------------------------------------------------------
% % Initializing -------------------------------------------------
% % parcing input
tumb = 'average';
datatype = class(im(1));
resolution = [];
num_strip = 1;
lmbd = [];
minmax = [];
nvxml = [];
if mod(nargin,2)
    error('Name/Value pairs are not matching!');
end
for c_var = 1:2:nargin-2
    switch varargin{c_var}
        case 'tumb'
            tumb = varargin{c_var+1};
        case 'datatype'
            datatype = varargin{c_var+1};
        case 'lmbd'
            if ~isnumeric(varargin{c_var+1})
                error('Input lmbd is not valid!');
            end
            lmbd = varargin{c_var+1};
        case 'minmax'
            if ~isnumeric(varargin{c_var+1})
                error('Input min/max is not valid!');
            end
            minmax = varargin{c_var+1};
        case 'resolution'
            if ~isnumeric(varargin{c_var+1}) || length(varargin{c_var+1})~=45
                error('Input resolution is not valid!');
            end
            resolution = varargin{c_var+1};
        case 'num_strip'
            if isempty(varargin{c_var+1}) || ~isnumeric(varargin{c_var+1}) || length(varargin{c_var+1})>1 || varargin{c_var+1}<1
                error('Input num_strip is not valid!');
            end
            num_strip = varargin{c_var+1};
        case 'nvxml'
            nvxml = varargin{c_var+1};
        otherwise
            error('Unknown given Name/Value pair!');
    end
end

config = 'BIP';
if isempty(tumb); tumb = 'none'; end
tif.res.x = [72,1]; tif.res.y = [72,1]; tif.res.unit = 2;
if ~isempty(resolution); tif.res.x = resolution(1:2); tif.res.y = resolution(3:4); tif.res.unit = resolution(5); end
siz.lines   = size(im,1);
siz.samples = size(im,2);
siz.bands   = size(im,3);
% % this is the formula in tiff manual to calculate num_strip from RowsPerStrip
% num_strip = floor((siz.lines + RowsPerStrip - 1) ./ RowsPerStrip); 
RowsPerStrip = ceil(siz.lines/num_strip);
last_strip = siz.lines - (num_strip-1)*RowsPerStrip; % number of rows in the last image strip (not necessarilly equal others)

ByteCount_tumb = length(typecast(uint8(0),'uint8')); % the Byte Count of a single unit of the tumb
ByteCount_im   = length(typecast(cast(im(1,1,1),datatype),'uint8')); % the Byte Count of a single unit of the im
if ~isempty(minmax)
    SMinSampleValue = minmax(1);
    SMaxSampleValue = minmax(2);
elseif isinteger(cast(0,datatype))
    SMinSampleValue = intmin(datatype);
    SMaxSampleValue = intmax(datatype);
else % single or double
    SMinSampleValue = 0;
    SMaxSampleValue = 1;
end
SMinSampleValue = cast(SMinSampleValue,datatype);
SMaxSampleValue = cast(SMaxSampleValue,datatype);
switch datatype
    case {'uint8','uint16','uint32','uint64'}
        SampleFormat = 1;
    case {'int8','int16','int32','int64'}
        SampleFormat = 2;
    case {'single','double'}
        SampleFormat = 3;
    otherwise
        SampleFormat = 4;
end
if isempty(lmbd)
    lmbd = 1:siz.bands;
end
if length(lmbd) ~= siz.bands
    warning('Input wavelenght range does not match the image-cube channels (3rd dimension size)!!');
end

%% ---------------------------------------------------------------
% % Making the tumb ---------------------------------------------
% prcntl = 99;
rgb = [];
if siz.bands < 3
    sens = [];
elseif ischar(tumb)
    switch lower(tumb)
        case 'none'
            sens = [];
        case 'flux'
            if siz.bands ~= 7; error('Given data size does not match FLUX camera! error in making the tumb image!'); end
            sens = [1,0,0,1,0,0,0; 0,1,0,0,1,0,0; 0,0,1,0,0,1,0]./2;
        case 'triangular' % triangular rgb sensitivities
            x = 0:(1/(1000-1)):1; b = -2*x+1; b(b<0) = 0; g1 = x; g2 = -x+1; g = min(g1,g2); r = 2*x -1; r(r<0) = 0;
            sens = [r;g;b];
            if siz.bands == 3
                sens = flip(sens,2);                
            end
            sens = interp1(x,sens',(0:(siz.bands-1))./(siz.bands-1)); sens = sens';
            sens = sens./repmat(sum(sens,2),[1, size(sens,2)]); % normalizing to sum = 1
            clear x r g g1 g2 b
        case 'average'
            sens = ones(3,size(im,3))./size(im,3);
        otherwise
            error('Unkown thumbnail is asked!')
    end
elseif isnumeric(tumb) % % NOTE! tumb's 1st dimension assumed to relate to RGB, and 2nd dim to band numbers
    if length(size(tumb)) == 3 % 3d tumb --> tumb is given as an RGB image
        if size(tumb,1) == size(im,1) && size(tumb,2) == size(im,2) && size(tumb,3) == 3
            sens = []; % rgb = tumb; ---> this is not applied here, but below, for the sake of code beauty!
        else
            error('Given thumbnail image size does not match the image cube spatial size!')
        end
    elseif length(tumb) == 1 || length(tumb(:)) == 3 % % the page-numbers to make the RGB, or the Grey image from
        if length(tumb) == 1; tumb = repmat(tumb,[3,1]); end
        sens = zeros(3,siz.bands);
        sens(1,tumb(1)) = 1; sens(2,tumb(2)) = 1; sens(3,tumb(3)) = 1;
    elseif length(tumb(:)) == siz.bands && siz.bands > 3 % % visible wavelength range is given --> use CIE-XYZ
        sens = 0;
    elseif size(tumb,1) == 3 && size(tumb,2) == siz.bands % % RGB weigths are given for each wavelength
        sens = tumb;
    else
        error('Unknown thumbnail option!');
    end
end
if isinteger(im)
    im_single = single(im)./single(intmax(class(im)));
elseif isfloat(im)
    im_single = single(im);
end
if ~sens
    rgb = t_spd2sth(im_single,tumb,'D65',1964,'RGB');
elseif ~isempty(sens)
    r = bsxfun(@times,im_single,permute(sens(1,:),[1,3,2])); r = sum(r,3); % r = r./max(r(:));
    g = bsxfun(@times,im_single,permute(sens(2,:),[1,3,2])); g = sum(g,3); % g = g./max(g(:));
    b = bsxfun(@times,im_single,permute(sens(3,:),[1,3,2])); b = sum(b,3); % b = b./max(b(:));
    rgb = cat(3,r,g,b); % rgb = rgb./max(rgb(:));
    clear r g b
end
% rgb = rgb./prctile(rgb(:),prcntl);
rgb = uint8((2^8-1)*rgb);
if length(size(tumb)) == 3 % 3d tumb --> tumb is given as an RGB image
    % rgb = uint8((2^8-1)*mat2gray(tumb));
    rgb = uint8(tumb);
end

flag_tumb = ~isempty(rgb);
if ~isempty(rgb)
    data.tumb = i_config(rgb,config,num_strip,RowsPerStrip);
end
data.im = struct([]);
for c_bands = 1:siz.bands
    data.im{c_bands,1} = i_config(im(:,:,c_bands),config,num_strip,RowsPerStrip);
end

%% --------------------------------------------------------------
% % Making the IFDs ---------------------------------------------
ifd = struct([]); tag = struct([]);
if flag_tumb
    c_tag = 1;
    ifd(1).num = []; ifd(1).name = {}; ifd(1).tag = []; ifd(1).type = []; ifd(1).len = []; ifd(1).value = {};
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'ImageWidth',       256,i_type('SHORT',   'nam2num'), 1,         siz.samples  }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'ImageLength',      257,i_type('SHORT',   'nam2num'), 1,         siz.lines    }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'BitsPerSample',    258,i_type('SHORT',   'nam2num'), 3,         ones(3,1).*8.*ByteCount_tumb }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'Compression',      259,i_type('SHORT',   'nam2num'), 1,         1            }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'PhotometricInterp',262,i_type('SHORT',   'nam2num'), 1,         2            }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'StripOffsets',     273,i_type('LONG',    'nam2num'), num_strip, zeros(num_strip,1)}); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1; % value is given later
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'SamplesPerPixel',  277,i_type('SHORT',   'nam2num'), 1,         3            }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'RowsPerStrip',     278,i_type('SHORT',   'nam2num'), 1,         RowsPerStrip }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'StripByteCounts',  279,i_type('LONG',    'nam2num'), num_strip, [ones(num_strip-1,1).*RowsPerStrip*siz.samples.*3.*ByteCount_tumb; last_strip*siz.samples.*3.*ByteCount_tumb] }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'XResolution',      282,i_type('RATIONAL','nam2num'), 1,         tif.res.x    }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'YResolution',      283,i_type('RATIONAL','nam2num'), 1,         tif.res.y    }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(1) = i_set_ifd(ifd(1),{c_tag, 'ResolutionUnit',   296,i_type('SHORT',   'nam2num'), 1,         tif.res.unit }); tag(1).(strcat('t',num2str(ifd(1).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
end

for c_ifd = flag_tumb + (1:siz.bands)
    c_tag = 1;
    ifd(c_ifd).num = []; ifd(c_ifd).name = {}; ifd(c_ifd).tag = []; ifd(c_ifd).type = []; ifd(c_ifd).len = []; ifd(c_ifd).value = {};
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'ImageWidth',       256  ,i_type('SHORT',   'nam2num'), 1,                      siz.samples  });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'ImageLength',      257  ,i_type('SHORT',   'nam2num'), 1,                      siz.lines    });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'BitsPerSample',    258  ,i_type('SHORT',   'nam2num'), 1,                      1.*8.*ByteCount_im }); tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'Compression',      259  ,i_type('SHORT',   'nam2num'), 1,                      1            });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'PhotometricInterp',262  ,i_type('SHORT',   'nam2num'), 1,                      1            });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'StripOffsets',     273  ,i_type('LONG',    'nam2num'), num_strip,              zeros(num_strip,1)});   tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1; % value is given later
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'SamplesPerPixel',  277  ,i_type('SHORT',   'nam2num'), 1,                      1            });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1; % not a necessary field for simple RGB images
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'RowsPerStrip',     278  ,i_type('SHORT',   'nam2num'), 1,                      RowsPerStrip });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'StripByteCounts',  279  ,i_type('LONG',    'nam2num'), num_strip,              [ones(num_strip-1,1).*RowsPerStrip.*siz.samples.*1.*ByteCount_im; last_strip.*siz.samples.*1.*ByteCount_im] }); tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'XResolution',      282  ,i_type('RATIONAL','nam2num'), 1,                      tif.res.x    });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'YResolution',      283  ,i_type('RATIONAL','nam2num'), 1,                      tif.res.y    });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'PlanarConfig',     284  ,i_type('SHORT',   'nam2num'), 1,                      1            });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'ResolutionUnit',   296  ,i_type('SHORT',   'nam2num'), 1,                      tif.res.unit });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'SampleFormat',     339  ,i_type('SHORT',   'nam2num'), ifd(c_ifd).value{tag(c_ifd).(strcat('t',num2str(277)))},repmat(SampleFormat,ifd(c_ifd).value{tag(c_ifd).(strcat('t',num2str(277)))})}); tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1; % this tag has len = SamplesPerPixel 
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'SMinSampleValue',  340  ,i_type(datatype,  'typ2num'), ifd(c_ifd).value{tag(c_ifd).(strcat('t',num2str(277)))},SMinSampleValue}); tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'SMaxSampleValue',  341  ,i_type(datatype,  'typ2num'), ifd(c_ifd).value{tag(c_ifd).(strcat('t',num2str(277)))},SMaxSampleValue}); tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1;
end
    c_ifd = flag_tumb + 1; % Added ONLY to the end of 1st ifd that comes after the thumbnail image (ie. the first image-cube page) 
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'WaveLength',     65000  ,i_type('FLOAT',   'nam2num'), length(lmbd),           lmbd         });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1; 
    if ~isempty(nvxml)
    ifd(c_ifd) = i_set_ifd(ifd(c_ifd),{c_tag, 'NVXML',          65111  ,i_type('ASCII',   'nam2num'), length(nvxml),          nvxml        });        tag(c_ifd).(strcat('t',num2str(ifd(c_ifd).tag(c_tag)))) = c_tag; c_tag = c_tag + 1; 
    end
%% --------------------------------------------------------------
% % Calculations ------------------------------------------------

% % calculating data sizes
siz.hdr = 8; % always fix
for c_ifd = 1:length(ifd)
    for c_tag = 1:length(ifd(c_ifd).tag)
        siz_tag_data = ifd(c_ifd).len(c_tag) * i_type(ifd(c_ifd).type(c_tag),'num2len');
        siz.d_tag(c_ifd,c_tag) = max(4 , siz_tag_data); % Note! siz.d_tag will be a squared matrix even if ifds do not have equal #-of-tags
    end
    siz.d_tag(siz.d_tag <= 4) = 0;
    siz.ifd(c_ifd,1) = 2 + length(ifd(c_ifd).tag) * 12 + 4;
    siz.d_ifd(c_ifd,1) = sum(siz.d_tag(c_ifd,:));
    
    % we find StripByteCounts from it's tag which is number 279
    for c_c = 1:length(ifd(c_ifd).tag)
        if ifd(c_ifd).tag(c_c) == 279; temp = ifd(c_ifd).value{c_c}; end
    end
    siz.strip(c_ifd,:) = temp;
    siz.im(c_ifd,1) = sum(siz.strip(c_ifd,:));
end

% % calculating data offsets
acc_sum_ifd = 0;
acc_sum_d_ifd = 0;
acc_sum_im = 0;
for c_ifd = 1:length(ifd)
    off.ifd(c_ifd) = siz.hdr + sum(siz.im) + acc_sum_ifd + acc_sum_d_ifd;
    acc_sum_d_tag = 0;
    for c_tag = 1:length(ifd(c_ifd).tag)
        off.d_tag(c_ifd,c_tag) = off.ifd(c_ifd) + siz.ifd(c_ifd) + acc_sum_d_tag;
        acc_sum_d_tag = acc_sum_d_tag + siz.d_tag(c_ifd,c_tag);
    end
    acc_sum_strip = 0;
    for c_strip = 1:num_strip
        off.strip(c_ifd,c_strip) = siz.hdr + acc_sum_im + acc_sum_strip;
        acc_sum_strip = acc_sum_strip + siz.strip(c_ifd,c_strip);
    end
    acc_sum_ifd   = acc_sum_ifd   + siz.ifd(c_ifd);
    acc_sum_d_ifd = acc_sum_d_ifd + siz.d_ifd(c_ifd);
    acc_sum_im    = acc_sum_im    + siz.im(c_ifd);
end
off.d_tag = off.d_tag .* logical(siz.d_tag);
% % giving the StripOffsets tag its values
for c_ifd = 1:length(ifd)
    for c_tag = 1:length(ifd(c_ifd).tag)
        if ifd(c_ifd).tag(c_tag) == 273 % StripOffsets tag
            ifd(c_ifd).value{c_tag} = nonzeros(off.strip(c_ifd,:));
        end
    end
end
%% --------------------------------------------------------------
% % Writing the file --------------------------------------------
machineformat = 'ieee-le'; % little endian --> the least significant BYTE of the word is written first and so on
if length(filename)<5
    filename = strcat(filename,'.tif');
elseif ~(strcmp(filename(end-3:end),'.tif') || strcmp(filename(end-4:end),'.tiff')) 
    filename = strcat(filename,'.tif');
end
fid = fopen(filename, 'w+', machineformat);
% % writitng Header -----------------------------
if strcmp(machineformat,'ieee-le')
    fwrite(fid, uint8(73), 'uint8');
    fwrite(fid, uint8(73), 'uint8');
elseif strcmp(machineformat,'ieee-be')
    fwrite(fid, uint8(77), 'uint8');
    fwrite(fid, uint8(77), 'uint8');
end
fwrite(fid, uint16(42), 'uint16');
fwrite(fid, uint32(off.ifd(1)), 'uint32');
% % writing Image data --------------------------
if ~isempty(rgb)
    c_write = 0;
    for c_strip = 1:num_strip
        c_write = c_write + fwrite(fid, data.tumb{c_strip}, 'uint8');
    end
end
c_write = 0;
for c_bands = 1:siz.bands
    for c_strip = 1:num_strip
        % % % Matlab Help: If you specify a bit precision, then fwrite saturates for all values outside the range.
        c_write = c_write + fwrite(fid, data.im{c_bands}{c_strip}, datatype);
    end
end

% % writing IFDs --------------------------------
for c_ifd = 1:length(ifd)
    fwrite(fid, uint16(length(ifd(c_ifd).tag)), 'uint16');
    for c_tag = 1:length(ifd(c_ifd).tag)
        fwrite(fid, uint16(ifd(c_ifd).tag(c_tag)), 'uint16');
        fwrite(fid, uint16(ifd(c_ifd).type(c_tag)), 'uint16');
        fwrite(fid, uint16(ifd(c_ifd).len(c_tag)), 'uint32');
        
        if siz.d_tag(c_ifd,c_tag) <= 4 
            fwrite(fid, uint32(ifd(c_ifd).value{c_tag}), 'uint32');
        else
            fwrite(fid, uint32(off.d_tag(c_ifd,c_tag)), 'uint32');
        end
    end
    if c_ifd < length(ifd) % not yet the last IFD
        fwrite(fid, uint32(off.ifd(c_ifd+1)), 'uint32'); % Offset to the next IFD
    else
        fwrite(fid, uint32(0), 'uint32'); % Last IFD, write 0
    end
    % % writing IFDs big data ------------------
    for c_tag = 1:length(ifd(c_ifd).tag)
        if siz.d_tag(c_ifd,c_tag) > 4
            type = i_type(ifd(c_ifd).type(c_tag),'num2typ');
            for c_data = 1:ifd(c_ifd).len(c_tag)
                if ~strcmpi(i_type(ifd(c_ifd).type(c_tag),'num2nam'),'RATIONAL') && ~strcmpi(i_type(ifd(c_ifd).type(c_tag),'num2nam'),'SRATIONAL')
                    fwrite(fid, ifd(c_ifd).value{c_tag}(c_data), type);
                elseif strcmpi(i_type(ifd(c_ifd).type(c_tag),'num2nam'),'RATIONAL')
                    fwrite(fid, ifd(c_ifd).value{c_tag}(2*c_data-1), 'uint32');
                    fwrite(fid, ifd(c_ifd).value{c_tag}(2*c_data  ), 'uint32');
                elseif strcmpi(i_type(ifd(c_ifd).type(c_tag),'num2nam'),'SRATIONAL')
                    fwrite(fid, ifd(c_ifd).value{c_tag}(2*c_data-1), 'int32');
                    fwrite(fid, ifd(c_ifd).value{c_tag}(2*c_data  ), 'int32');
                end
            end
        end
    end
end
fclose(fid);
end
    
%% --------------------------------------------------------------
% % internal functions ------------------------------------------
function out = i_config(in,config,num_strip,RowsPerStrip)
last_strip = size(in,1) - (num_strip-1)*RowsPerStrip;
if strcmpi(config,'BIP')
    out = cell(num_strip,1);
    % out = [zeros( (num_strip-1), size(in,3) * size(in,2) * RowsPerStrip, class(in)); ...
    %        zeros( 1,             size(in,3) * size(in,2) * last_strip  , class(in))];
    BIP = permute(in,[3,2,1]);
    for c_strip = 1:(num_strip-1)
        temp = BIP(:,:,1+(RowsPerStrip*(c_strip-1)):(RowsPerStrip*c_strip));
        out{c_strip,1} = temp(:)'; % each row is the (RGB,RGB,...) values of one row of original image cube
    end
    temp = BIP(:,:, (end-last_strip+1):end); % 1+(RowsPerStrip*(c_strip-1)):(RowsPerStrip*c_strip));
    out{end,1} = temp(:)'; % each row is the (RGB,RGB,...) values of one row of original image cube
elseif strcmpi(config,'BIL')
elseif strcmpi(config,'BSQ')
end
end
% % -------------------------------------------------------------
function ifd = i_set_ifd(ifd,data)
ifd.num(end+1,1)   = data{1};
ifd.name{end+1,1}  = data{2};
ifd.tag(end+1,1)   = data{3};
ifd.type(end+1,1)  = data{4};
ifd.len(end+1,1)   = data{5}; % this is #-of-values, not the memory space length of the tag
ifd.value{end+1,1} = data{6};
end
% % -------------------------------------------------------------
function out = i_type(in,request)
% % lookup table for type_num/type_length/type_name
% % request = 'num2len'/'num2typ'/'num2nam'/'nam2num'/'nam2len'/'nam2typ'/'typ2num'
LUT = {
    1, 1, 'uint8',  'BYTE';     % 8-bit unsigned integer
    2, 1, 'uint8',  'ASCII';    % 8-bit containing a 7-bit ASCII code + last BYTE must be NUL
    3, 2, 'uint16', 'SHORT';    % 16-bit (2-byte) unsigned integer
    4, 4, 'uint32', 'LONG';     % 32-bit (4-byte) unsigned integer
    5, 8, 'uint64', 'RATIONAL'; % two LONGs, 1st numerator, 2nd denominator
    6, 1, 'int8',   'SBYTE';    % 8-bit signed (2s-complement) integer
    7, 1, 'uint8',  'UNDEFINED';% 8-bit byte, contains anything, depends on the definition of the field
    8, 2, 'int16',  'SSHORT';   % 16-bit (2-byte) signed (2s-complement) integer
    9, 4, 'int32',  'SLONG';    % 32-bit (4-byte) signed (2s-complement) integer
    10,8, 'int64',  'SRATIONAL';% two SLONGs, 1st numerator, 2nd denominator
    11,4, 'single', 'FLOAT';    % Single precision (4-byte) IEEE format
    12,8, 'double', 'DOUBLE'    % Double precision (8-byte) IEEE format
    }; % num/len/typ/nam
if strcmp(request(1:3),'nam') % for nam2sth
    for c_type = 1:size(LUT,1)
        if strcmp(LUT{c_type,4},in)
            sel = c_type;
        end
    end
end
% % typ2num is only used to find the datatype of tags 339-340-341
% % in general typ2nam is not 1-to-1 (ex. for uint8) and here we only consider the numerical cases
% % with this code below for uint8 which is not 1-to-1, always num=1 which is the only numerical case, is chosen
if strcmp(request(1:3),'typ') % for typ2num
    for c_type = 1:size(LUT,1)
        if strcmp(LUT{c_type,3},in)
            sel = c_type;
            break % always only 1st uint8 is selected 
        end
    end
end
switch lower(request)
    case 'num2len'
        out = LUT{in,2};
    case 'num2typ'
        out = LUT{in,3};
    case 'num2nam'
        out = LUT{in,4};
    case 'nam2num'
        out = LUT{sel,1};
    case 'nam2len'
        out = LUT{sel,2};
    case 'nam2typ'
        out = LUT{sel,3};
    case 'typ2num' % This is not obviously 1-to-1, it is only used to find the datatype of tags 339-340-341
        out = LUT{sel,1};
end
end
% % -------------------------------------------------------------



