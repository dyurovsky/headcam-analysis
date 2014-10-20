function [freq_vals,rad_power] = radial_power(X,freq_vals)

X = double(X);

[d1,d2] = size(X);
if d1 ~= d2
    error('X must be square!')
end
if max(X(:))>1
    X = X./max(X(:));
end
X = X - mean(X(:));

DIM = [d1,d2];

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]';
u = repmat(u,1,DIM(2));

v = u'; % because d1=d2

u = fftshift(u);
v = fftshift(v);

[theta,f] = cart2pol(u,v);

win = hanning(d2) * hanning(d1)'; % hanning window for FT
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


