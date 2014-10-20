ages = [3 8 12];
nsets = [3 4 3];
colors = [0.25 0.25 0.25; 0.5 0.5 0.5; 0.75 0.75 0.75];

meas_clutter = false;
meas_phasecoher = false;
meas_power = true;
meas_orient = true;
meas_edgedens = false;


nFramesToJump = 30; % to sample every nFramesToJump
startFrame = 60;
std_thresh = 40; % ### ?? set intelligently!

%%

for ca = 1:length(ages)
    for cs = 1:nsets(ca)
        dirname = sprintf('../videoData/AslinBabyCam/Frames/wesse.%dmos.%d',ages(ca),cs);
        cliplist = dir(dirname);
        mkdir(sprintf('stats/%dmos.%d',ages(ca),cs));
        
        for cc = 1:length(cliplist)
            if ~strcmp(cliplist(cc).name(1),'.')
                imlist = dir(sprintf('%s/%s/*.jpg',dirname,cliplist(cc).name));
                % this is generally a very long list
                imsToCheck = startFrame:nFramesToJump:length(imlist);
                imcnt = 1;
                
                for ci = imsToCheck
                    
                    clear im;
                    im = imread([dirname,'/',cliplist(cc).name,'/',imlist(imcnt).name]);
                    imgray = double(rgb2gray(im));
                    imstd(imcnt) = std(imgray(:));
                    imstdIx(imcnt) = ci;
                    clipFrameStr = sprintf('stats/%dmos.%d/Clip%s.mat',ages(ca),cs,cliplist(cc).name);
                    if ci == imsToCheck(1)
                        save(clipFrameStr,'imstd','imstdIx','imsToCheck');
                    else
                        save(clipFrameStr,'imstd','imstdIx','imsToCheck','-append');
                    end
                    
                    if std(imgray(:)) > std_thresh  % should probably use a mean threshold as well or instead (more sensitive to dark bad frames)
                        
                        [h,w] = size(imgray);
                        imgray = imgray(:,11:end-10); % remove black vertical borders
                        
                        % analyses on the rectangular image
                        if meas_edgedens
                            bw = edge(imgray,'canny');
                            edgedens(imcnt) = sum(bw(:))/numel(bw);
                            save(clipFrameStr,'edgedens','ci','imsToCheck','-append');
                        end
                        
                        if meas_clutter
                            [tmp,tmp,clutter_val(imcnt)] = SIclutter(im,exponent_val);
                            save(clipFrameStr,'clutter_val','ci','imsToCheck','-append');
                        end
                        
                        if meas_phasecoher
                            % sharpness index/GPC
                            sharpness(imcnt) = sharpness_index(imgray);
                            gpc(imcnt) = global_phase_coherence(imgray);
                            save(clipFrameStr,'gpc','sharpness','ci','imsToCheck','-append');
                        end
                        
                        % analysis on square center image
                        if meas_power
                            sq_ix = (1:h) + (w-h)/2;
                            imgraysq = imgray(:,sq_ix);
                            % power exponent & power by orientation
                            [freq_vals(imcnt,:),rad_power(imcnt,:)] = radial_power(imgraysq);
                            
                            x = freq_vals(imcnt,2:end-1); % ignore DC & highest s.f. (pixel boundaries)
                            y = rad_power(imcnt,2:end-1); % ### MAYBE remove more?
                            [B,BINT,R,RINT,STATS] = regress(log(y)',[log(x+eps)' ones(length(y),1)]);
                            alph(imcnt) = B(1);
                            Rsqrd(imcnt) = STATS(1);
                            p_val(imcnt) = STATS(3);
                            %plot(log(x+eps),log(y))
                            %hold on; plot(log(x+eps),B(1).*(log(x+eps))+B(2),'k-')
                            save(clipFrameStr,'freq_vals','rad_power','alph','Rsqrd','p_val','ci','imsToCheck','-append');
                        end
                        
                        if meas_orient
                            % orientation analysis
                            [orient_vals(imcnt,:),energy_orient(imcnt,:)] = orientation_power_new(imgraysq);
                            save(clipFrameStr,'energy_orient','orient_vals','ci','imsToCheck','-append');
                        end
                        
                    end
                    imcnt = imcnt + 1;
                end
            end
        end
    end
end

%% Plot results

agecolors = {'r' 'g' 'b'};
facecolors = [0.5 0 0; 0 0.5 0; 0 0 0.5];

% For orientation analyses
figure;
hold on;
set(gca,'FontSize',12)
plot([45 45],[0 7e-2],'k--')
plot([90 90],[0 7e-2],'k--')
plot([135 135],[0 7e-2],'k--')

for ca = 1:length(ages)
    hold on;
    cnt = 1;
    
    for cs = 1:nsets(ca)
        
        datafiles = dir(sprintf('stats/%dmos.%d/*.mat',ages(ca),cs));
        
        for df = 1:length(datafiles)
            load(sprintf('stats/%dmos.%d/%s',ages(ca),cs,datafiles(df).name));
            orient_vals = orient_vals(1,:);
            goodIx = imstd>std_thresh;
            ener = energy_orient(goodIx,:);
            ener_tot_by_row = sum(ener,2);
            clip_orient_prob{ca,cnt} = ener./repmat(ener_tot_by_row,[1,size(ener,2)]);
            n_samp{ca,cnt} = sum(goodIx);
            
            all_means(cnt,:) = mean(clip_orient_prob{ca,cnt});
            all_means_N(cnt) = n_samp{ca,cnt};
            
            cnt = cnt + 1;
        end
    end
    
    all_means_N = all_means_N./sum(all_means_N);
    mean_val = sum(all_means.*repmat(all_means_N',[1 size(all_means,2)]));
    sum(mean_val)
    SEM = std(mean_val)/sqrt(size(all_means,1));
    ub = mean_val + SEM;
    lb = mean_val - SEM;
    
    h = fill([orient_vals,orient_vals(end:-1:1)],[lb ub(end:-1:1)],'k');
    set(h,'LineStyle','none','FaceAlpha',0.2,'FaceColor',facecolors(ca,:));
    plot(orient_vals,mean_val,'Color',agecolors{ca},'LineWidth',1.5)
    text(145, 0.06 - (ca-1).*0.005, sprintf('%d mos.',ages(ca)), 'Color', agecolors{ca},'FontSize',12);
    
end
ylabel('probability')
xlabel('orientation')
set(gca,'XTick',0:45:180)
xlim([0 180])
ylim([0 7e-2])