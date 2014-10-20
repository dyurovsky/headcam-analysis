ages = [3 8 12];
nsets = 3;
colors = [0.25 0.25 0.25; 0.5 0.5 0.5; 0.75 0.75 0.75];

%% loop over images to get estimate of group's clutter exponent
N = [];
K = [250 500 750 1000 1500 2000];
figure; hold on;
for ca = ages  
    for cs = 1:nsets
        dirname = sprintf('samples/%dmo.%d',ca,cs);
        imlist = dir(sprintf('%s/*.jpg',dirname));
        for ci = 1:length(imlist)
            im = imread([dirname,'/',imlist(ci).name]);
            N = [N; SIclutter(im,-1.3)];
            plot(K,N(end,:),'k*')
        end
    end
end
Ktmp = repmat(K,[size(N,1),1]);
[B,BINT,R,RINT,STATS] = regress(log(Ktmp(:)),[log(N(:)) ones(numel(Ktmp),1)]);
err_var = STATS(end);
exponent_val = B(1);
prop_const = exp(B(2));
save('segment_info','K','N','err_var','exponent_val','prop_const')

%%
meas_clutter = false;
meas_phasecoher = false;
meas_power = true;
 
cnt = 1;
for ca = 12
    for cs = 2:nsets
        dirname = sprintf('samples/%dmo.%d',ca,cs);
        imlist = dir(sprintf('%s/*.jpg',dirname));
        
        for ci = 1:length(imlist)-9
            im = imread([dirname,'/',imlist(ci).name]);
            imgray = rgb2gray(im);
            [h,w] = size(imgray);
            imgray = imgray(:,11:end-10); % remove black vertical borders
            
            % analyses on the rectangular image
            bw = edge(imgray,'canny');
            edgedens(ci,cs,cnt) = sum(bw(:))/numel(bw);
            if meas_clutter
                [tmp,tmp,clutter_val(ci,cs,cnt)] = SIclutter(im,exponent_val);
            end
            
            if meas_phasecoher
                % sharpness index/GPC
                sharpness(ci,cs,cnt) = sharpness_index(imgray);
                gpc(ci,cs,cnt) = global_phase_coherence(imgray);
            end
            
            % analysis on square center image
            if meas_power
                sq_ix = (1:h) + (w-h)/2;
                imgraysq = imgray(:,sq_ix);
                % power exponent & power by orientation
                [freq_vals(ci,:),rad_power(ci,:)] = radial_power(imgraysq);
                
                x = freq_vals(ci,2:end-1); % ignore DC & highest s.f. (pixel boundaries)
                y = rad_power(ci,2:end-1);
                [B,BINT,R,RINT,STATS] = regress(log(y)',[log(x+eps)' ones(length(y),1)]);
                alph(ci,cs,cnt) = B(1);
                Rsqrd(ci,cs,cnt) = STATS(1);
                p_val(ci,cs,cnt) = STATS(3);
                %plot(log(x+eps),log(y))
                %hold on; plot(log(x+eps),B(1).*(log(x+eps))+B(2),'k-')
                
                % orientation analysis
                [orient_vals,energy_orient(ci,cs,cnt,:)] = orientation_power(imgraysq);
            end
            
        end
        all_rad_power{cs,ca} = rad_power;
    end
    cnt = cnt + 1;
end
        
%%
load power_spec_info
load clutter_info
load phaseco_info
load edgedens_info

nims = 10;
for ca = 1:length(ages)  
    figure;
    hold on;
    subplot(3,2,1)
    hold on;
    plot([45 45],[0 14e-3],'k--')
    plot([90 90],[0 14e-3],'k--')
    plot([135 135],[0 14e-3],'k--')
    for cs = 1:nsets
        datatmp = squeeze(energy_orient(:,cs,ca,:));
        datatmpm = mean(datatmp);
        SEM = std(datatmp)/sqrt(size(datatmp,1));
        ub = datatmpm + SEM;
        lb = datatmpm - SEM;
        h = fill([orient_vals,orient_vals(end:-1:1)],[lb ub(end:-1:1)],'k');
        set(h,'LineStyle','none','FaceAlpha',0.2,'FaceColor',[0.7 0.7 0.7]);
        plot(datatmpm,'k-','LineWidth',1.5)
    end
    title(sprintf('%d Mos',ages(ca)));
    ylabel('probability')
    xlabel('orientation')
    set(gca,'XTick',0:45:180)
    xlim([0 180])
    ylim([0 14e-3])
    
    subplot(3,2,2);
    for cs = 1:nsets
        datatmp = all_rad_power{cs,ages(ca)};
        datatmpm = mean(datatmp);
        SEM = std(datatmp)/sqrt(size(datatmp,1));
        ub = datatmpm + SEM;
        lb = datatmpm - SEM;
        %h = fill([freq_vals(1,:),freq_vals(1,end:-1:1)],[lb ub(end:-1:1)],'k');
        set(h,'LineStyle','none','FaceAlpha',0.2)
        loglog(freq_vals(1,:),datatmpm,'k-','LineWidth',1)
        hold on;
        %meanalph = alph(ci,cs,:);
    end
    ylim([0 10^8])
    xlim([0 10^2.5])
    title(sprintf('%d Mos',ages(ca)));
    xlabel('Freq (cycles/image)');
    ylabel('Power');
        
    subplot(3,2,3);
    hold on;
    datatmp = edgedens(:,:,ca);
    hist(datatmp(:),0:0.01:0.2)
    xlabel('edge density')
    ylabel('count')
    xlim([0 0.2])
    
    subplot(3,2,4);
    hold on;
    datatmp = alph(:,:,ca);
    ptmp = p_val(:,:,ca);
    datatmp = datatmp(ptmp<0.05);
    hist(datatmp(:),-4:0.1:-2);
    xlim([-4 -2]);
    xlabel('spectral slope')
    ylabel('count')
    
    subplot(3,2,5)
    hold on;
    datatmp = clutter_val(:,:,ca);
    hist(datatmp(:),0.5e4:1500:2e4)
    xlabel('clutter index')
    ylabel('count')
    
    subplot(3,2,6)
    hold on;
    datatmp = gpc(:,:,ca);
    hist(datatmp(:),0:3e4:3e5)
    xlabel('global phase consistency')
    ylabel('count')
    
    print('-dpsc2','Baby_cam_static_summary.ps','-append');
    
end
        
%loglog(freq_vals',rad_power','k-','Color',colors(cs,:),'LineWidth',0.5);
%hold on;
%loglog(freq_vals(1,:),mean(rad_power),'k-','LineWidth',1);
%xlabel('spatial freq in cyc/im')
%ylabel('power');
%title(sprintf('%dmo.%d',ca));
%print('-depsc2',sprintf('figs/%dmo.eps',ca))











