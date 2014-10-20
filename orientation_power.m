function [Xth,TH,Xsf,SF] = orientation_power(X)

X = double(X);

[d1,d2] = size(X);
if d1 ~= d2
    error('X must be square!')
end
if max(X(:))>1
    X = X./max(X(:));
end
X = X - mean(X(:));

SCALE = 4; % amount to scale image masks by for anti-aliasing

%%% handle image dimensions and image center
[ydim,xdim] = size(X);                                         % image dimensions
if( mod(xdim,2) == 0 );  
    midx = xdim/2 + 1;                                       % if xdim is even, x center is right of center line
else
    midx = ceil(xdim/2);                                     % if xdim is odd,  x center is center column
end
if( mod(ydim,2) == 0 );  
    midy = ydim/2 + 1;                                       % if ydim is even, y center is below center line
else
    midy = ceil(ydim/2);                                     % if ydim is odd,  y center is center row
end
[xramp,yramp] = meshgrid( [1:SCALE*xdim]-SCALE*midx, [1:SCALE*ydim]-SCALE*midy ); % centered coords at each point, supersize for antialiasing

%%% compute full spectrum
win = hanning(ydim) * hanning(xdim)';          % hanning window for FT
Fim = fftshift( abs( fft2( win.*X ) ) );       % FT of image
%Fim = abs(Fim).*abs(Fim);


%%% compute power spectrum as a function of spatial frequency
% NOTE: NEED TO FIX TO CONVERT TO CYCLES PER DEG BEFORE CALCULATION
% so that width and steps are in meaningful units
dist                = sqrt( xramp.^2 + yramp.^2 );      % radial distance from center at each point, will use this to select rings in FT  
steps               = 2*SCALE;                          % spacing between rings
wk                  = steps*2.5;                        % width of rings, rings overlap to slightly smooth data
c                   = 1;                                % counter

for k = 0 : steps : max(dist(:))   
    ind             = find( abs(dist-k) <= wk );                    % select a ring of values at k with width wk
    mask            = zeros( size(dist) );                          % create image mask, 4x image size
    mask(ind)       = 1;                                            % set ring values to 1 in mask
    mask            = imresize( mask, 1/SCALE, 'bicubic' );         % create anti-aliased mask by resizing to image size
    mask            = max(mask,0);                                  % clamp mask min to zero
    mask            = min(mask,1);                                  % clamp mask max to one
    Xsf(c)          = k;                                            % store sf value
    SF(c)           = sum( Fim(:).*mask(:) ) / sum( mask(:) );      % mean energy in image
    c               = c + 1;
end

%Xsf                = 30 * Xsf / max(Xsf);                          % convert distance to spatial frequency in cycles per deg
SF                  = SF/sum(SF);                                   % make probability disitribution -- sum to 1


%%% compute power spectrum as a function of orientation
theta               = 90 + 180/pi * atan( yramp./xramp );           % compute FT angle at each point
wk                  = 10;%2.5;                                          % width of wedges, overlap to slightly smooth data
c                   = 1;                                            % counter

for k = min(theta(:)) : max(theta(:))
    if( k < wk );           ind = find( (abs(theta-k) <= wk) | abs(theta-k) >= 180-wk );    % select a wedge at k with width wk, handle 0/180 transition
    elseif( k > 180-wk );   ind = find( (abs(theta-k) <= wk) |(abs(theta-k) >= 180-wk) );   % select a wedge at k with width wk, handle 0/180 transition
    else                    ind = find( abs(theta-k) <= wk );                               % select a wedge at k with width wk
    end
    
    mask            = zeros( size(theta) );                         % create image mask, 4x image size
    mask(ind)       = 1;                                            % set wedge values to 1 in mask
    mask            = imresize( mask, 1/SCALE, 'bicubic' );         % create anti-aliased mask by resizing to image size
    mask            = max(mask,0);                                  % clamp mask min to zero
    mask            = min(mask,1);                                  % clamp mask max to one
    Xth(c)          = k;                                            % store angle
    TH(c)           = sum( Fim(:).*mask(:) ) / sum( mask(:) );      % mean energy in image
    c               = c + 1;   
end

TH                  = TH/sum(TH);                                   % make probability disitribution -- sum to 1


% figure(); colormap(jet);
% subplot(221); imagesc(X); axis image off; colorbar; title('image');
% subplot(222); imagesc(Fim); axis image off; colorbar; title('fourier spectrum');
% 
% subplot(223);
% plot( Xsf, SF, 'k' ); hold on; xlabel('spatial frequency');
% 
% subplot(224);
% plot( Xth, TH, 'k' ); hold on; xlabel('orientation');
