function [im, tumb, ifd, tif] = t_tif_read(filename, varargin)
% % function [im, tumb, ifd, tif] = t_tif_read(filename, bands)
% % --------------------------------------------------------------
% % 27.Mar.2017 (c) Arash Mirhashemi
% % A tiff reader for spectral image cubes (in combination with tif_write.m)
% % --------------------------------------------------------------
% % Inputs:
% % filename = the input tif file to be read
% % bands    = a vector that specifies which bands to load. It does not need to be ordered. 
% %            bands are counted starting by 1. Thumbnail image is not included here in counting, and it will be loaded anyway separately.
% %            if left empty, all the bands will load.
% %            if no valid number is given, no bands get loaded, only thumbnail image will load if it exists.
% %            (ex. if you need ONLY the thumbnail image you can set bands = [-1] for example, or any other invalid value)
% % Outputs:
% % im       = the output image cube
% % tumb     = the thumbnail image if exist 
% % ifd      = the image-file-directory(ifd) of the input tif file
% % tif      = the parsed ifd
% % --------------------------------------------------------------
% % usage example: [im, tumb, ifd, tif] = t_read_tif('example.tif');
% % --------------------------------------------------------------
% % Note-1 IMPORTANT! Although this tif-reader follows tif standards and can read any tif file, when it comes to preparing the final
% %        image cube, it follows (t_write_tif.m). This means it stack all ifd images (except 1st one if it is an thumbnail image) and makes a 3d image cube
% %        Any changes in the (t_write_tif.m) assumptions about the image-cube should be reflected here as well.
% % Note-2 Any tag in the ifds is read in this code, however tag names are recognized by TAG_NUM lookuptable in sub-function i_type
% %        if the tag is not in that list it will be given the name "Unknown" although its value is correctly read and outputed
% %        the current TAG_NUM contains TAG_NUM = [256,257,258,259,262,273,277,278,279,282,283,284,296,339,340,341,65000,65111];
% %        if any tag is added to (t_write_tif.m) it should be added here as well, otherwise it will be read as Unknown
% % Note-3 The wavelength range of the image-cube can be accessed like this:
% %        lmbd = tif(1-or-2 if tumb exist).WaveLength ( ex: lmbd = tif(1+~isempty(tumb)).WaveLength )
% % Note-4 The Minimum/Maximum values of the image-cube can be accessed like this:
% %        minimum = tif(1-or-2 if tumb exist).SMinSampleValue ( ex: lmbd = tif(1+~isempty(tumb)).SMinSampleValue )
% %        maximum = tif(1-or-2 if tumb exist).SMaxSampleValue ( ex: lmbd = tif(1+~isempty(tumb)).SMaxSampleValue )
% % Note-5 Natural Vision XML data --> If this exists in the file it can be accessed by: 
% %        NV_XML = char(ifd(1+~isempty(tumb)).value{ifd(1+~isempty(tumb)).tag==65111});

%% ---------------------------------------------------------------
% % Initializing -------------------------------------------------
bands = []; 

if mod(nargin-1,2)
    error('Name/Value pairs are not matching!');
end
for c_var = 1:2:nargin-2
    switch varargin{c_var}
        case 'bands'
            bands = varargin{c_var+1};
        otherwise
            error('Unknown given Name/Value pair!');
    end
end

TAG_NUM = [256,257,258,259,262,273,277,278,279,282,283,284,296,339,340,341,65000];

fid = fopen(filename,'r');
machineformat = fread(fid, 2, 'uint8'); % reading 7373 or 7777
if machineformat(1) == 73 && machineformat(2) == 73
    machineformat = 'ieee-le';
elseif machineformat(1) == 77 && machineformat(2) == 77
    machineformat = 'ieee-be';
end

if isempty(machineformat)
    fclose(fid);
    error('Not a valid TIFF file!!');
end
idnumber = fread(fid, 1, 'uint16', 0, machineformat); % reading 2A00
if idnumber ~= 42
    fclose(fid);
    error('Not a valid TIFF file!!');
end
% % Reading ifds -------------------------------------------------
offset = fread(fid, 1, 'uint32', 0, machineformat); % offset of 1sf ifd
fseek(fid,offset,'bof');
c_ifd = 0;
ifd = struct([]);
ifd(1).num = []; ifd(1).name = {}; ifd(1).tag = []; ifd(1).type = []; ifd(1).len = []; ifd(1).value = {};
while offset ~= 0
    c_ifd = c_ifd + 1;
    for c_tag = 1:fread(fid, 1, 'uint16', 0, machineformat); % #-of-tags in the current ifd
        ifd(c_ifd,1).num(c_tag,1)  = c_tag;
        ifd(c_ifd,1).tag(c_tag,1)  = fread(fid, 1, 'uint16', 0, machineformat);
        ifd(c_ifd,1).type(c_tag,1) = fread(fid, 1, 'uint16', 0, machineformat);
        ifd(c_ifd,1).len(c_tag,1)  = fread(fid, 1, 'uint32', 0, machineformat);
        ifd(c_ifd,1).value{c_tag,1}= fread(fid, 1, 'uint32', 0, machineformat);
        ifd(c_ifd,1).name{c_tag,1} = i_tag(ifd(c_ifd).tag(c_tag));
    end
    offset = fread(fid, 1, 'uint32', 0, machineformat); % offset of next ifd or 0 if it is the last
    if offset ~= 0
        fseek(fid,offset,'bof');
    end
end
% % Reading ifd large data ---------------------------------------
for c_ifd = 1:length(ifd)
    for c_tag = 1:length(ifd(c_ifd).num)
        if ifd(c_ifd).len(c_tag) * i_type(ifd(c_ifd).type(c_tag), 'num2len') > 4
            fseek(fid,ifd(c_ifd).value{c_tag},'bof');
            for c_len = 1:ifd(c_ifd).len(c_tag)
                if ifd(c_ifd).type(c_tag) == 5 
                    RATIONAL = fread(fid, 2, 'uint32', 0, machineformat);
                    ifd(c_ifd).value{c_tag}(c_len) = double(RATIONAL(2)) + 1/double(RATIONAL(1));
                elseif ifd(c_ifd).type(c_tag) == 10
                    SRATIONAL = fread(fid, 2, 'int32', 0, machineformat);
                    ifd(c_ifd).value{c_tag}(c_len) = double(SRATIONAL(2)) + 1/double(SRATIONAL(1));
                else
                    ifd(c_ifd).value{c_tag}(c_len) = fread(fid, 1, i_type(ifd(c_ifd).type(c_tag),'num2typ'), 0, machineformat);
                end
            end
        end
    end
end

% % Reading image data -------------------------------------------
flag_tumb = 0;
tif = struct([]);
for c_ifd = 1:length(ifd)
    % % Parsing each ifd's tags into tif (we need tif because tif field names are tag-names but ifd field names are the items in each tag)
    for c_tag = 1:length(TAG_NUM)
        tif(c_ifd,1).(i_tag(TAG_NUM(c_tag))) = i_pars(ifd(c_ifd),i_tag(TAG_NUM(c_tag)));
    end
    if isempty(tif(c_ifd).SampleFormat); tif(c_ifd).SampleFormat = 1; end
    if isempty(tif(c_ifd).SamplesPerPixel); 
        tif(c_ifd).SamplesPerPixel = 1; 
        % warning('SamplesPerPixel not specified - assumed to be 1!'); 
    end
    tif(c_ifd).ByteCount = tif(c_ifd).StripByteCounts(1)/tif(c_ifd).RowsPerStrip/tif(c_ifd).ImageWidth/tif(c_ifd).SamplesPerPixel;
    switch tif(c_ifd).SampleFormat
        case 1 % {'uint8','uint16','uint32','uint64'}
            tif(c_ifd).datatype = strcat('uint',num2str(tif(c_ifd).ByteCount * 8));
        case 2 % {'int8','int16','int32','int64'}
            tif(c_ifd).datatype = strcat('int',num2str(tif(c_ifd).ByteCount * 8));
        case 3 % {'single','double'}
            if tif(c_ifd).ByteCount == 4
                tif(c_ifd).datatype = 'single';
            elseif tif(c_ifd).ByteCount == 8
                tif(c_ifd).datatype = 'double';
            end
        otherwise % SampleFormat = 4
            error('Unkown data-type (SampleFormat = 4) !')
    end
    if c_ifd == 1 && tif(c_ifd).SamplesPerPixel == 3 && strcmpi(tif(c_ifd).datatype, 'uint8') && isempty(tif(1).WaveLength)
        flag_tumb = 1;
    end
end
if isempty(bands); bands = 1:length(tif); end;
bands(bands < 1) = []; 
if flag_tumb % the two for-loops (above and below) could be merged but then im size was not known beforehand and would change in each loop (=slow)
    bands = bands + 1; % 1st band is tumb so we add one to 'bands'
    bands(bands > length(ifd)) = []; bands = unique(bands);
    % im = zeros([tif(1).ImageLength, tif(c_ifd).ImageWidth, length(tif)-1], tif(2).datatype);
    im = zeros([tif(1).ImageLength, tif(c_ifd).ImageWidth, length(bands)], tif(2).datatype);
    if isempty(bands); bands = 1; end; % only tumb will be loaded
    tumb = zeros([tif(1).ImageLength, tif(c_ifd).ImageWidth, 3], 'uint8');
else
    bands(bands > length(ifd)) = []; bands = unique(bands);
    tumb = [];
    % im = zeros([tif(1).ImageLength, tif(c_ifd).ImageWidth, length(tif)], tif(1).datatype);
    im = zeros([tif(1).ImageLength, tif(c_ifd).ImageWidth, length(bands)], tif(1).datatype);
end

c_count = 0;
for c_ifd = 1:length(ifd)
        data = zeros([tif(c_ifd).ImageLength, tif(c_ifd).ImageWidth, tif(c_ifd).SamplesPerPixel],tif(c_ifd).datatype);
    for c_strip = 1:length(tif(c_ifd).StripOffsets)
        fseek(fid,tif(c_ifd).StripOffsets(c_strip),'bof');
        strip = fread(fid,tif(c_ifd).StripByteCounts(c_strip)./length(typecast(cast(0,tif(c_ifd).datatype),'uint8')), tif(c_ifd).datatype);
        strip_row_num = tif(c_ifd).StripByteCounts(c_strip)/tif(c_ifd).ImageWidth/tif(c_ifd).SamplesPerPixel/tif(c_ifd).ByteCount;
        strip = i_config(strip,'BIP',[strip_row_num, tif(c_ifd).ImageWidth, tif(c_ifd).SamplesPerPixel]); 
        data(((c_strip-1)*tif(c_ifd).RowsPerStrip+1):((c_strip-1)*tif(c_ifd).RowsPerStrip+strip_row_num),:,:) = cast(strip,tif(c_ifd).datatype);
    end
    if (c_ifd == 1) && flag_tumb
        tumb = data;
    elseif ismember(c_ifd,bands)
        c_count = c_count+1;
        im(:,:,(c_count)) = data;
    end
end
fclose(fid);
end
%% --------------------------------------------------------------
% % internal functions ------------------------------------------
function out = i_config(in,config,siz)
if strcmpi(config,'BIP')
    out = permute(reshape(in,[siz(3),siz(2),siz(1)]),[3,2,1]);
elseif strcmpi(config,'BIL')
elseif strcmpi(config,'BSQ')
end
end
% % -------------------------------------------------------------
function out = i_pars(ifd,tagname)
for c_tag = 1:length(ifd.name)
    if strcmpi(tagname,ifd.name{c_tag})
        out = ifd.value{c_tag};
        break
    else % at the end will happen if the tag does not exist
        out = [];
    end
end
end
% % -------------------------------------------------------------
function out = i_tag(in)
LUT = { 
    256,  'ImageWidth';
    257,  'ImageLength';
    258,  'BitsPerSample';
    259,  'Compression';
    262,  'PhotometricInterp';
    273,  'StripOffsets';
    277,  'SamplesPerPixel';
    278,  'RowsPerStrip';
    279,  'StripByteCounts';
    282,  'XResolution';
    283,  'YResolution';
    284,  'PlanarConfig';
    296,  'ResolutionUnit';
    339,  'SampleFormat';
    340,  'SMinSampleValue';
    341,  'SMaxSampleValue';
    % % reusable tags (65000-65535)
    65000,'WaveLength'
    65111,'NaturalVisionXML'
    };
for c_type = 1:size(LUT,1)
    if LUT{c_type,1} == in
        out = LUT{c_type,2};
        break
    else
        out = 'UnknownTag'; % not found in the above LUT
    end
end
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
        if strcmpi(LUT{c_type,3},in)
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
%% --------------------------------------------------------------




