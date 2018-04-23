function [ mobr, nair, ttt4, ppt4 ] = fcn_vpbrox( r, ttt, ppt, dswt, mono, moo3, mobro )
% 
% Shaojie Song, 03/30/2018
% This is a matlab routine used to estimate the vertical profiles of BrOx and obtain monthly, diurnal profiles (Holmes et al. 2010; Legrand et al. 2016; Sherwen et al. 2016)

%% Br, pptv
    
    % calculate Br based on photochemical steady state (Holmes et al., 2010)
    
    % (R31) Br + O3 --> BrO + O2
    % (R32) BrO + hv --> Br + O (JBrO is the BrO photolysis frequency)
    % (R33) BrO + NO --> Br + NO2
    % [Br]/[BrO] = ( JBrO + k33*[NO] ) / ( k31*[O3] )
    
    % vertical profile of temperature (K)
    m1 = [0, 31,59, 90,120,151,181,212,243,273,304,334];
    m2 = [31,59,90,120,151,181,212,243,273,304,334,365];
    nlay = size(ttt,2);
    ttt2 = reshape(ttt(end+1-87600:end,:),[10,8760,nlay]);
    ttt3 = reshape(squeeze(mean(ttt2,1)),[24,365,nlay]);
    ttt4 = zeros(12,24,nlay); %initialize
    for cnt1=1:1:12
        ttt4(cnt1,:,:) = squeeze(mean(ttt3(:,m1(cnt1)+1:m2(cnt1),:),2));
    end
    
    % vertical profile of pressure (hPa)
    ppt2 = reshape(ppt(end+1-87600:end,:),[10,8760,nlay]);
    ppt3 = reshape(squeeze(mean(ppt2,1)),[24,365,nlay]);
    ppt4 = zeros(12,24,nlay); %initialize
    for cnt1=1:1:12
        ppt4(cnt1,:,:) = squeeze(mean(ppt3(:,m1(cnt1)+1:m2(cnt1),:),2));
    end
    
    % Downward shortwave radiation (W/m2)
    dsw2 = mean(reshape(dswt(end+1-10*24*365:end),[10,8760]),1);
    dsw3 = reshape(dsw2,[24,365]);
    dsw4 = zeros(24,12);
    for cnt=1:1:12
        dsw4(:,cnt) = mean(dsw3(:,m1(cnt)+1:m2(cnt)),2);
    end
    
    % local air number density [molec/cm3 LTP]
    nair = ppt4 .* 1d2 ./ ( 0.029d0 .* r .* ttt4 ) .* 6.02d23 ./ 1d6;
    
    % convert DSW to JBrO (s-1)
    % relationship between DSW and JBrO, see j_brox_hox.m
    jbro = 3.8947e-08 .* dsw4 .* dsw4 + 9.8059e-05 .* dsw4;
    
    % rate constants and concentrations
    k31 = 1.6d-11 .* exp( -780d0 ./ ttt4 ); %cm3 molecule-1 s-1 from JPL
    k33 = 8.8d-12 .* exp(  260d0 ./ ttt4 ); %cm3 molecule-1 s-1 from JPL
    
    % compute br mixing ratios, pptv
    mobr = zeros(size(mono))-999;
    for cnt1=1:1:12
        for cnt2=1:1:24
            for cnt3=1:1:nlay
                mobr(cnt1,cnt2,cnt3) = mobro(cnt1,cnt2,cnt3) * ...
                                       ( jbro(cnt2,cnt1) + k33(cnt1,cnt2,cnt3) * mono(cnt1,cnt2,cnt3) * nair(cnt1,cnt2,cnt3) / 1d12 ) ./ ...
                                       ( k31(cnt1,cnt2,cnt3) * moo3(cnt1,cnt2,cnt3) * nair(cnt1,cnt2,cnt3) / 1d9 );
            end
        end
    end