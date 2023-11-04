%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Probabilistic Seismic Risk Analysis Using Upper Tail Approximation of Site-specific Seismic Hazard Curves   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

NumSimu = 10^6;

% DATA_2905 : 0.1-degree grids
% DATA_18714: EA locations
load Data_EQrisk_Malawi
DATA = DATA_2905;
clear DATA_2905 DATA_18714

F = 1./[100 200 500 750 1000 2000 2500 5000 10000];

F_pick = round(NumSimu*(1-F));

site_pick = 1:2905

% [A-GI	A-LD A-DS B-GI B-LD	B-DS C-GI C-LDC-DS]
FragilityPara(1,:) = [0.159 0.127 0.152 0.179 0.152 0.177 0.278 0.236 0.272];
FragilityPara(2,:) = [0.401	0.351 0.365	0.439 0.381	0.410 0.412	0.389 0.399];
FragilityPara(3,:) = [0.333 0.334 0.333 0 0 0 0 0 0]; % User specified - maintain the weight sums to 1

count = 0;
for ii = 1:9
    FragilityIndex(count+1:round(sum(NumSimu*FragilityPara(3,1:ii))),1) = ii;
    count = round(sum(NumSimu*FragilityPara(3,1:ii)));
end

% Optional: Randomize the order of FragilityIndex.
FragilityIndex(:,2) = rand(NumSimu,1);
FragilityIndex = sortrows(FragilityIndex,2);

ExProb = 1-(1:NumSimu)/(NumSimu+1);

visual_check = 0

for ii = 1:length(site_pick)
    
    % Step 1: Fit four models and choose the best-fitting model based on the linear correlation coefficient.
    for jj = 1:4
        [~,~,RRtmp(jj,:)] = ModelFit_MalawiSeismicHazard(1-F',DATA(site_pick(ii),3:11)',jj);
    end
    [~,model_pick] = max(RRtmp(:,1));

    % Step 2: Obtain the model parameters for the best-fitting model.
    [para_im,~,RR_im] = ModelFit_MalawiSeismicHazard(1-F',DATA(site_pick(ii),3:11)',model_pick);

    % Step 3: Simulate annual maximum PGA values from the best-fitting model.
    if model_pick == 1
        disp(['Site ',num2str(site_pick(ii)),': Lognormal distribution with R2 = ',num2str(RR_im(1))]);
        pga = exp(para_im(1)+para_im(2)*randn(NumSimu,1));
    elseif model_pick == 2
        disp(['Site ',num2str(site_pick(ii)),': Gumbel distribution with R2 = ',num2str(RR_im(1))]);
        pga = para_im(1)-log(-log(rand(NumSimu,1)))/para_im(2);
    elseif model_pick == 3
        disp(['Site ',num2str(site_pick(ii)),': Frechet distribution with R2 = ',num2str(RR_im(1))]);
        pga = exp(log(para_im(1))-log(-log(rand(NumSimu,1)))/para_im(2));
    elseif model_pick == 4
        disp(['Site ',num2str(site_pick(ii)),': Weibull distribution with R2 = ',num2str(RR_im(1))]);
        pga = exp(log(para_im(1))+log(-log(1-rand(NumSimu,1)))/para_im(2));
    end
    
    % Step 4: Calculate collapse probability
    for ijk = 1:9
        index = find(FragilityIndex(:,1) == ijk);
        if ~isempty(index)
            pc(index,1) = normcdf(log(pga(index,1)/FragilityPara(1,ijk))/FragilityPara(2,ijk));
        end
        clear index
    end
    
    % Step 5: Visually compare the goodness-of-the-fit.
    pga = sort(pga);
    pc  = sort(pc);
    
    if visual_check == 1
        
        figure (1)
        loglog(DATA(site_pick(ii),3:11),F,'rs',pga,ExProb,'b-'); hold on;
        xlabel('PGA'); ylabel('Exceedance Probability'); grid on; axis square; legend('Data','Fit'); axis([0.01 10 10^-5 10^-1])

        figure (2)
        semilogy(DATA(site_pick(ii),15:23),F,'cs'); hold on; % Use equally weighted - permanent
        semilogy(DATA(site_pick(ii),24:32),F,'ms'); hold on; % Use equally weighted - semi-permanent
        semilogy(DATA(site_pick(ii),33:41),F,'rs'); hold on; % Use equally weighted - traditional
        semilogy(pc,ExProb,'b-'); hold on;  
        xlabel('PGA'); ylabel('Exceedance Probability'); grid on; axis square; legend('Data P','Data SP','Data T','Fit'); axis([0 1 10^-5 10^-1]);
        
        clf (1)
        clf (2)
    
    end
    
    % Step 6: Store relevant results
    EQRISK(ii,1) = model_pick;
    EQRISK(ii,2) = para_im(1);
    EQRISK(ii,3) = para_im(2);

    EQRISK(ii,4:12) = pga(F_pick);

    EQRISK(ii,13:21) = pc(F_pick);
    
    clear RRtmp model_pick para_im RR_im pga pc
    
end

%% Check

figure (101)
plot(EQRISK(:,1),DATA(site_pick,12),'b.'); axis square;

figure (102)
plot(EQRISK(:,2),DATA(site_pick,13),'b.'); axis square;

figure (103)
plot(EQRISK(:,3),DATA(site_pick,14),'b.'); axis square;

figure (104)
subplot(331); plot(EQRISK(:,4) ,DATA(site_pick,3) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(332); plot(EQRISK(:,5) ,DATA(site_pick,4) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(333); plot(EQRISK(:,6) ,DATA(site_pick,5) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(334); plot(EQRISK(:,7) ,DATA(site_pick,6) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(335); plot(EQRISK(:,8) ,DATA(site_pick,7) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(336); plot(EQRISK(:,9) ,DATA(site_pick,8) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(337); plot(EQRISK(:,10),DATA(site_pick,9) ,'b.',[0 1],[0 1],'r-'); axis square;
subplot(338); plot(EQRISK(:,11),DATA(site_pick,10),'b.',[0 1],[0 1],'r-'); axis square;
subplot(339); plot(EQRISK(:,12),DATA(site_pick,11),'b.',[0 1],[0 1],'r-'); axis square;

figure (105)
% index_start = 14
% index_start = 23
index_start = 32
subplot(331); plot(EQRISK(:,13),DATA(site_pick,index_start+1),'b.',[0 1],[0 1],'r-'); axis square;
subplot(332); plot(EQRISK(:,14),DATA(site_pick,index_start+2),'b.',[0 1],[0 1],'r-'); axis square;
subplot(333); plot(EQRISK(:,15),DATA(site_pick,index_start+3),'b.',[0 1],[0 1],'r-'); axis square;
subplot(334); plot(EQRISK(:,16),DATA(site_pick,index_start+4),'b.',[0 1],[0 1],'r-'); axis square;
subplot(335); plot(EQRISK(:,17),DATA(site_pick,index_start+5),'b.',[0 1],[0 1],'r-'); axis square;
subplot(336); plot(EQRISK(:,18),DATA(site_pick,index_start+6),'b.',[0 1],[0 1],'r-'); axis square;
subplot(337); plot(EQRISK(:,19),DATA(site_pick,index_start+7),'b.',[0 1],[0 1],'r-'); axis square;
subplot(338); plot(EQRISK(:,20),DATA(site_pick,index_start+8),'b.',[0 1],[0 1],'r-'); axis square;
subplot(339); plot(EQRISK(:,21),DATA(site_pick,index_start+9),'b.',[0 1],[0 1],'r-'); axis square;




