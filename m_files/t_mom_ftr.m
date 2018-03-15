function out = t_mom_ftr(im, mom_ftr, varargin)
% out = t_mom_ftr(im, mom_ftr)
% % --------------------------------------------------------------
% % 27.Mar.2017 (c) Arash Mirhashemi
% % Moment Features for Spectral Data
% % --------------------------------------------------------------
% % Inputs:
% % im         3D spectral image (spectrum in 3rd dimension) or 2D matrix of spectral vectors (spectra in columns)
% %            input 'im' is assumed in [0,1] range, and for further calculations, it is force-casted to 'double' without any checking
% % mom_ftr    The cell array of the requested moment features, noted by the 4-character notation
% %            1st character is the moment type M/S
% %            2nd character is raw/central moment o/i
% %            3rd character is original or L1-normalized spectrum R/r
% %            4th character is the moment's order
% %            Example: {'Mor3','SiR2'}
% % Name-------Value------------------------(optional)
% % lmbd       the wavelength range of the spectral data from which the moments are calculated
% %            default value is 1:size(im,3) that is "Natural numbers"
% % 
% % Outputs:   
% % out        The gray-scale feature image of the input spectral image, 
% %            or the vector of the feature values for an input 2D spectral data (spectra in columns)
% %            
% % --------------------------------------------------------------
% % NOTE: The feature image can be visualized with this line of code: imshow( uint8((2^8-1)*mat2gray(ftr_im)) );
% % --------------------------------------------------------------

im = double(im);
[ysize, xsize, zsize] = size(im);
lmbd = permute(1:zsize,[1 3 2]); 
mom_ftr = cellstr(mom_ftr);
out = zeros(ysize, xsize, length(mom_ftr));

if mod(nargin,2)
    error('Name/Value pairs are not matching!');
end
for c_var = 1:2:nargin-2
    switch varargin{c_var}
        case 'lmbd'
            lmbd = varargin{c_var+1}(:);
            lmbd = permute(lmbd',[1 3 2]); 
        otherwise
            error('Unknown given Name/Value pair!');
    end
end
lmbd = repmat(lmbd,[ysize, xsize]);
M0 = sum(im,3);                  M0   = repmat(M0,  [1, 1, zsize]);
M1 = sum( lmbd.*im ,3);          M1   = repmat(M1,  [1, 1, zsize]);

for c_mtd = 1:length(mom_ftr);
    order = str2double(mom_ftr{c_mtd}(4));
    switch mom_ftr{c_mtd}(1:3)
        case 'MoR' % Raw
            S = sum( lmbd.^order.*im ,3);
        case 'MiR' % Ram
            S = sum( (lmbd-(M1./M0)).^order.*im ,3);
        case 'Mor' % paw
            S = sum( lmbd.^order.*im./M0 ,3);
        case 'Mir' % pam
            S = sum( (lmbd-(M1./M0)).^order.*im./M0 ,3);
        case 'SoR' % Rbw
            S = sum( im.^order ,3);
        case 'SiR' % Rbm
            S = sum( (im-(M0./zsize)).^order ,3);
        case 'Sor' % pbw
            S = sum( (im./M0).^order ,3);
        case 'Sir' % pbm
            S = sum( ((im./M0)-1/zsize).^order ,3);
        case 'NoR' % Rcw
            S = sum( im.^order .* lmbd ,3);
        case 'NiR' % Rcm
            S = sum( (im-M1/(zsize*(zsize+1)/zsize)).^order .* lmbd ,3);
        case 'Nor' % pcw
            S = sum( (im./M0).^order .* lmbd ,3);
        case 'Nir' % pcm
            S = sum( (im./M0 - M1./M0/(zsize*(zsize+1)/zsize) ).^order .* lmbd ,3);   
        otherwise
            error('Unknown Method Requested!');
    end
    out(:,:,c_mtd) = S;
    
end
%% --------------------------------------------------------------







