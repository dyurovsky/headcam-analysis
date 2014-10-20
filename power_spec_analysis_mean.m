auth = '/Volumes/Crop_Auth/';
fake = '/Volumes/Crop_Imit/';
a_im = dir([auth,'0*png']);
f_im = dir([fake,'*png']);
f_im = f_im(6:10);

for im_no = 1:13
    
    if im_no < 9
        filedir = auth;
        filename = a_im(im_no).name;
    else
        filedir = fake;
        filename = f_im(im_no-8).name;
        
    end
    im = imread([filedir,filename]);
    im = gray_hughes(im);
    [d1,d2] = size(im);
    if d1>d2
        im = im';
        [d1,d2] = size(im);
        needtotrans = true;
    else
        needtotrans = false;
    end
    im = imresize(im,[800 nan]);
    im = make_dims_even(im);
    [d1,d2] = size(im);
    % split image into two overlapping square images
    if d1<d2
        s1 = im(1:d1,1:d1);
        s2 = im(1:d1,end-d1+1:end);
    end
    
    if needtotrans
        s1 = s1';
        s2 = s2';
    end
    
    s1 = double(s1);
    s2 = double(s2);
    
    sw1 = whiten_hughes(s1);
    sw2 = whiten_hughes(s2);
    
    [fv_s1(im_no,:),rp_s1(im_no,:)] = radial_power(s1);
    [fv_s2(im_no,:),rp_s2(im_no,:)] = radial_power(s2);
    [fv_sw1(im_no,:),rp_sw1(im_no,:)] = radial_power(sw1);
    [fv_sw2(im_no,:),rp_sw2(im_no,:)] = radial_power(sw2);
    
    figure;
    loglog(fv_s1(im_no,:),rp_s1(im_no,:),'b')
    hold on;
    loglog(fv_s2(im_no,:),rp_s2(im_no,:),'m')
    
    loglog(fv_sw2(im_no,:),rp_sw2(im_no,:),'g')
    loglog(fv_sw1(im_no,:),rp_sw1(im_no,:),'c')
    ylim([1e10,1e20])
    title(filename);
    print('-djpeg100',sprintf('FreqAnalysis/800_%s',filename(1:end-4)));
    close;
    
end
%%
figure;
loglog(fv_s1(1,:),mean([rp_s1(1:8,:);rp_s2(1:8,:)]),'b');
hold on;
loglog(fv_s1(1,:),mean([rp_s1(9:13,:);rp_s2(9:13,:)]),'k');


%%
X = rp_s1;
y = [ones(8,1); zeros(5,1)];
b = glmfit(X,y,'binomial');

% 2. yhat = glmval(b,X,'logit')
% where X is the 2nd half of the data.

X2 = rp_s2;
yhat = glmval(b,X2,'logit');

c = yhat;
c(yhat>0.5) = 1;
c(yhat<=0.5) = 0;








