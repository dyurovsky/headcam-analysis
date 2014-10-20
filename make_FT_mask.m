function mask   = make_FT_mask(x,ind,sc)

mask            = zeros( size(x) );                             % create image mask, should be 4x image size
mask(ind)       = 1;                                            % set ring values to 1 in mask
mask            = imresize( mask, 1/sc, 'bicubic' );            % create anti-aliased mask by resizing to image size
mask            = max(mask,0);                                  % clamp mask min to zero
mask            = min(mask,1);                                  % clamp mask max to one