function [Xth,TH,Xsf,SF] = orientation_power_new(im,maxcpd)

im = double(im);

[d1,d2] = size(im);
if d1 ~= d2
    error('im must be square!')
end
if max(im(:))>1
    im = im./max(im(:));
end
im = im - mean(im(:));

if nargin<2; maxcpd = 30; % default assumption of 1 arcmin per pixl
elseif isempty(maxcpd); maxcpd = 30; end


[ydim,xdim]                     = size(im);                                         % image dimensions
if( mod(xdim,2) == 0 );  midx   = xdim/2 + 1;                                       % if xdim is even, x center is right of center line
else                     midx   = ceil(xdim/2); end                                 % if xdim is odd,  x center is center column
if( mod(ydim,2) == 0 );  midy   = ydim/2 + 1;                                       % if ydim is even, y center is below center line
else                     midy   = ceil(ydim/2); end                                 % if ydim is odd,  y center is center row

sc                              = 4;                                                % up-sample ramps by this much (to create anti-aliased masks)
[xramp,yramp]                   = meshgrid( [1:sc*xdim]-sc*midx, [1:sc*ydim]-sc*midy ); % centered coords at each point, supersize for antialiasing

win                             = hanning(ydim) * hanning(xdim)';                   % hanning window for FT
im                              = im .* win;                                        %
im                              = im - mean(im(:));                                 % zero mean the image
Fim                             = fftshift( abs( fft2( fftshift(im) ) ) );     % FT of image

%%% compute power spectrum as a function of spatial frequency
% NOTE: NEED TO FIX TO CONVERT TO CYCLES PER DEG BEFORE CALCULATION
% so that width and steps are in meaningful units
dist                = sqrt( xramp.^2 + yramp.^2 );                  % radial distance from center at each point, will use this to select rings in FT
[ydim,xdim]         = size( dist );
maxdist             = min( xdim, ydim ) / 2;
maxdist             = maxcpd*maxdist/max(dist(:));

mindist             = (2*5*maxcpd)/size(im,1);          % min dist is 5 cycles per image
dist                = maxcpd*dist/max(dist(:));

step = 'log';

if strcmp(step,'lin');
    
    wk                  = 1;
    steps               = 1;
    c                   = 1;                                            % counter
    
    for k = mindist+wk : steps : maxdist - wk
        ind             = find( abs(dist-k) <= wk );                    % select a ring of values at k with width wk
        mask            = make_FT_mask(dist,ind,sc);                       % make ring mask for range of frequencies
        Xsf(c)          = k;                                            % store sf value
        SF(c)           = sum( Fim(:).*mask(:) ) / sum( mask(:) );      % mean energy in image
        c               = c + 1;
        %imagesc(mask); axis image; colormap gray; truesize; pause; 
    end
    
elseif strcmp(step,'log')
    
    edges = logspace(log10(mindist),log10(maxdist),16);
    
    c                   = 1;                                            % counter
    
    for k = 2:length(edges) - 1
        
        ind             = find( dist >= edges(k-1) & dist <= edges(k+1) );                    % select a ring of values at k with width wk
        
        %ind             = find( abs(dist-k) <= wk );                    % select a ring of values at k with width wk
        mask            = make_FT_mask(dist,ind,sc);                       % make ring mask for range of frequencies
        Xsf(c)          = edges(k);                                            % store sf value
        SF(c)           = sum( Fim(:).*mask(:) ) / sum( mask(:) );      % mean energy in image
        c               = c + 1;
        %display(num2str(k));
        %imagesc(mask); axis image; colormap gray; truesize; drawnow;
        
    end
    
end


%%% compute power spectrum as a function of orientation
theta               = 90 + 180/pi * atan( yramp./xramp );           % compute FT angle at each point
wk                  = 5; % width of wedges, overlap to slightly smooth data
steps               = 5;
c                   = 1;                                            % counter
dist2                = sqrt( xramp.^2 + yramp.^2 );                  % radial distance from center at each point, will use this to select rings in FT
[ydim,xdim]         = size( dist2 );

for k = min(theta(:)) : steps : max(theta(:))
    if( k < wk );           ind = find( ((dist >= mindist & dist2<=(min(xdim,ydim)/2)) & ((abs(theta-k) <= wk) | abs(theta-k) >= 180-wk)) );    % select a wedge at k with width wk, handle 0/180 transition
    elseif( k > 180-wk );   ind = find( ((dist >= mindist & dist2<=(min(xdim,ydim)/2)) & ((abs(theta-k) <= wk) | abs(theta-k) >= 180-wk)) );   % select a wedge at k with width wk, handle 0/180 transition
    else                    ind = find( ((dist >= mindist & dist2<=(min(xdim,ydim)/2)) & (abs(theta-k) <= wk)) );                               % select a wedge at k with width wk
    end
    
    mask            = make_FT_mask(dist,ind,sc);                       % make wedge mask for range of angles
    Xth(c)          = k;                                            % store angle
    TH(c)           = sum( Fim(:).*mask(:) ) / sum( mask(:) );      % mean energy in image
    c               = c + 1;
    %imagesc(mask); axis image; drawnow;
end