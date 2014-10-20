function [freq_vals,rad_power,theta_vals,orient_power] = power_analyses(X,freq_vals)

X = double(X);

[d1,d2] = size(X);
if d1 ~= d2
    error('X must be square!')
end

DIM = [d1,d2];

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]';
u = repmat(u,1,DIM(2));

v = u'; % because d1=d2

u = fftshift(u);
v = fftshift(v);

[theta,f] = cart2pol(u,v);

X = X - mean(X(:));
if max(X(:))>1
    X = X./max(X(:));
end
win = hanning(d2) * hanning(d1)';
xf = fftshift(fft2(X.*win));
xf = xf(:);
pf = abs(xf).*abs(xf);

if nargin<2
    freq_vals = 0:(d1/2);
end

for r = freq_vals
    ix{r + 1} = find(f == r);
end
rad_power = nan(size(freq_vals));
for r = freq_vals
    if ~isempty(ix{r+1})
        rad_power(r + 1) = mean( pf( ix{r+1} ) );
    end
end

theta = theta.*(180/pi);
theta_vals = 0:1:180;
wk = 2.5;
cnt = 1;
for k = theta_vals
    if( k < wk );           
        ind = find( (abs(theta-k) <= wk) | abs(theta-k) >= 180-wk );    % select a wedge at k with width wk, handle 0/180 transition
    elseif( k > 180-wk );   
        ind = find( (abs(theta-k) <= wk) |(abs(theta-k) >= 180-wk) );   % select a wedge at k with width wk, handle 0/180 transition
    else
        ind = find( abs(theta-k) <= wk );                               % select a wedge at k with width wk
    end
    mask = zeros(size(theta)); % create image mask
    mask(ind) = 1;
    mask = max(mask,0); % clamp mask min to zero
    mask = min(mask,1);
    orient_power(cnt) = sum( pf.*mask(:) ) / sum( mask(:) );
    cnt = cnt + 1;
end


