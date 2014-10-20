auth = '/kyb/agmbshare/image_databases/Bruegel/Crop_Auth/';
fake = '/kyb/agmbshare/image_databases/Bruegel/Crop_Imit/';
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
    im = imresize(im,0.5);
    im = make_dims_even(im);
    [d1,d2] = size(im);
    if d1>d2
        im = im';
        [d1,d2] = size(im);
        needtotrans = true;
    else
        needtotrans = false;
    end
    % split image into two overlapping square images
    if d1<d2
        s1 = im(1:d1,1:d1);
        s2 = im(1:d1,end-d1+1:end);
    end
    
    s1 = double(s1);
    s2 = double(s2);
    
    sw1 = whiten_hughes(s1);
    sw2 = whiten_hughes(s2);
    
    [fv_s1,rp_s1] = radial_power(s1);
    [fv_s2,rp_s2] = radial_power(s2);
    [fv_sw1,rp_sw1] = radial_power(sw1);
    [fv_sw2,rp_sw2] = radial_power(sw2);
    
    figure;
    loglog(fv_s1,rp_s1,'b')
    hold on;
    loglog(fv_s2,rp_s2,'m')
    
    loglog(fv_sw2,rp_sw2,'g')
    loglog(fv_sw1,rp_sw1,'c')
    ylim([1e10,1e20])
    title(filename);
    print('-djpeg100',sprintf('FreqAnalysis/%s',filename(1:end-4)));
    close;
    
end

    