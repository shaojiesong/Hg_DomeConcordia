%% main.m
    
    % This is the main routine for the Hg Dome C box model
    % Shaojie SONG, 06/10/2016, initiated
    % Shaojie SONG, 03/30/2018, finalized
    
%% Set up
    
    close all; clear; clc
    
    global ROU EDPTH FACSTP2LTP HG0DVD HG2DVD DT TS BETA DD %define global variables
    
    % load data
    load( './data/mdatain_bro_monthly.mat' ) %BrO
    upHg = csvread( './data/free_troposphere_hg_conc.csv', 1, 0 ); %free trop Hg0 and HgII
    
    % select a combination of oxidation mechanism (0 not included; 1 lower bound; 2 upper bound)
    LBROX = 2;
    LOHOX = 0;
    LO3OX = 0;
    
    % select a combination of oxidants concentration levels (0 lower bound; 1 best estimate; 2 upper bound)
    LBRMR = 2;
    LOHMR = 1;
    LO3MR = 1;
    
    % Reduction of snow mercury
    TAU1 = 14; %lifetime of snow Hg2 against photo-reduction in the summer period, day
    TAU2 = 365; %lifetime of snow Hg2 against dark-reduction in the non-summer period, day
    
    % select MAR data
    LMAR = 1;
    
    % variables
    HG0DVD = 1d-6; %hg0 dry deposition velocity (m/s)
    HG2DVD = 1d-2; %hg2 dry deposition velocity (m/s)
    BETA = 3d-3; %scaling factor of turbulent diffusion in surface snow
    
    % integration time periods
    DT = 1; %time step [min]
    TS = (0+DT):DT:1440; %tspan
    
    % snow properties
    ROU = 300.; %snow density, kg/m3, ref: Dommergue et al. 2012 ACP
    EDPTH = 0.2; %surface snow 1~2 fold e-depth, m, ref: France et al. 2011 ACP
    
    % the uncertainties of O3, HOx and BrOx are 2%, 50%, and a factor of 2.5
    fac_br = 2.5;
    fac_oh = 0.5;
    fac_o3 = 0.02;
    
    % parameters
    R = 287d0; %molar gas constant
    N0 = 1013.25d2/(0.029d0*R*273.15d0)*6.02d23/1d6; %number air density STP [molec/cm3]
    
%% meteorological parameters
    
    % load MAR data
    switch LMAR
        case 0
            load( './data/mdatain_MARorg.mat' )
        case 1
            load( './data/mdatain_MARadj.mat' )
    end
    
    % adjust according to average surface TP
    tmpTT = mean(TTt(end+1-87600:end,1));
    tmpPP = mean(PPt(end+1-87600:end,1));
    Nair = double(tmpPP*1d2/(0.029d0*R*tmpTT)*6.02d23/1d6); %local air density LTP [molec cm-3]
    FACSTP2LTP = double((6.02d23/201d0/1d9/1d6)*(Nair/N0)); %ng/m3 STP Hg --> molec cm-3 LTP Hg
    
    % interpolate and extrapolate
    XX1 = (0.5:1:239.5)./10.;
    XXm = (DT/2:DT:1440-DT/2)./60;
    XXh = (0.5:1:23.5);
    nday = numel(DSWt)/numel(XX1);
    nlay = size(Kzt,2);
    
    % altitude [m]
    ZZ = double( mean(ZZt(end+1-87600:end,:)) )';
    DD = double( ZZ(1:end)-[0;ZZ(1:end-1)] ); %depth of each layer [m]
    
    % DownSW [W/m2]
    DSWt2 = reshape(DSWt,[numel(XX1),nday]);
    DSWtm = zeros(numel(XXm),nday);
    for cnt = 1:1:nday
        DSWtm(:,cnt) = interp1(XX1,DSWt2(:,cnt),XXm,'linear','extrap');
    end
    
    % effective Kz [m2/s] and TKE [m2/s2]
    tmpfac = (PPt./tmpPP)./(TTt./tmpTT);
    Kze    = double( Kzt .* tmpfac .* tmpfac );
    Kze2   = reshape(Kze,[numel(XX1),nday,nlay]);
    Kzem   = zeros(numel(XXm),nday,nlay);
    for cnt1 = 1:1:nlay
        for cnt2 = 1:1:nday
            Kzem(:,cnt2,cnt1) = interp1(XX1,Kze2(:,cnt2,cnt1),XXm,'linear','extrap');
        end
    end
    
    TKEt2 = reshape(TKEt,[numel(XX1),nday]);
    TKEtm = zeros(numel(XXm),nday);
    for cnt = 1:1:nday
        TKEtm(:,cnt) = interp1(XX1,TKEt2(:,cnt),XXm,'linear','extrap');
    end
    
    % avoid negative values
    Kzem(Kzem < 1d-5) = 1d-5; %[m2/s]
    TKEtm(TKEtm < 1d-2) = 1d-2; %[m2/s2]
    
%% chemical tracers' concentrations
    
    % load oxidants' concentrations
    switch LMAR
        case 0
            load( './data/mdatain_OXorg.mat' )
        case 1
            load( './data/mdatain_OXadj.mat' )
    end
    
    % apply uncertainty levels
    switch LOHMR
        case 0
            moOH  = moOH  .* (1-fac_oh);
            moHO2 = moHO2 .* (1-fac_oh);
        case 2
            moOH  = moOH  .* (1+fac_oh);
            moHO2 = moHO2 .* (1+fac_oh);
    end
    
    switch LO3MR
        case 0
            moO3 = moO3 .* (1-fac_o3);
        case 2
            moO3 = moO3 .* (1+fac_o3);
    end
    
    % BrOx = Br + BrO, pptv
    moBrO = zeros(size(moNO))-999; %BrO, pptv
    for cnt1=1:1:12
        moBrO(cnt1,:,:) = avg_bro(cnt1);
    end
    [ moBr, Nair3, TTt3, PPt3 ] = fcn_vpbrox( R, TTt, PPt, DSWt, moNO, moO3, moBrO );
    
    % apply uncertainty levels
    switch LBRMR
        case 0
            error('Error found NO case 0 for Br')
        case 2
            moBr  = moBr  .* fac_br;
            moBrO = moBrO .* fac_br;
    end
    
    % Convert units to molec cm-3 LTP
    moO3c  = moO3  .* Nair3 ./ 1d9; %ppbv --> molec cm-3 LTP
    moNOc  = moNO  .* Nair3 ./ 1d12; %pptv --> molec cm-3 LTP
    moNO2c = moNO2 .* Nair3 ./ 1d12; %pptv --> molec cm-3 LTP
    moBrc  = moBr  .* Nair3 ./ 1d12; %pptv --> molec cm-3 LTP
    moBrOc = moBrO .* Nair3 ./ 1d12; %pptv --> molec cm-3 LTP
    moOHc  = moOH  .* Nair3 ./ N0; %molec cm-3 STP --> molec cm-3 LTP
    moHO2c = moHO2 .* Nair3 ./ N0; %molec cm-3 STP --> molec cm-3 LTP
    
%% Free tropospheric Hg concentrations
    
    % adjust free tropospheric Hg concentration (ng/m3 STP)
    kohdiff = 8;    kbrdiff = 4; %correction factors for high oxidation rates
    if LBROX == 1
        upHg0 = upHg(:,1)-upHg(:,3);
        upHg2 = upHg(:,3);
    elseif LBROX == 2
        upHg0 = upHg(:,1).*((kbrdiff.*(upHg(:,3)./(upHg(:,1)-upHg(:,3)))+1).^-1);
        upHg2 = upHg(:,1)-upHg0;
    elseif LBROX == 0 && ( LOHOX == 1 || LO3OX == 1 )
        upHg0 = upHg(:,2)-upHg(:,4);
        upHg2 = upHg(:,4);
    elseif LBROX == 0 && ( LOHOX == 2 || LO3OX == 2 )
        upHg0 = upHg(:,2).*((kohdiff.*(upHg(:,4)./(upHg(:,2)-upHg(:,4)))+1).^-1);
        upHg2 = upHg(:,2)-upHg0;
    end
    
    % Convert units to molec cm-3 LTP
    upHg0m = FACSTP2LTP .* fcn_reshape2( upHg0, DSWtm );
    upHg2m = FACSTP2LTP .* fcn_reshape2( upHg2, DSWtm );
    
%% Hg0 gas-phase oxidation mechanism, rate constants
    
    % -------------------- two-step Br reaction --------------------
    % (R1) Hg0 + Br + ([M]) --> HgBr + ([M])          unit: cm3 molecule-1 s-1
    % (R2) HgBr + ([M]) --> Hg0 + Br + ([M])          unit: s-1 (k2=k1/keq12)
    % (R3) HgBr + Br --> Hg0 + Br2                    unit: cm3 molecule-1 s-1
    % (R4) HgBr + X --> HgII (X = OH/Br/HO2/NO2/BrO)  unit: cm3 molecule-1 s-1
    
    k1max = 3.2d-12 .* ( Nair3 ./ N0 ); %Ariya et al. (2002)
    k1min = 1.46d-32 .* ( ( TTt3 ./ 298d0 ) .^ (-1.86d0) ) .* Nair3; %Donohoue et al (2006)
    keq12 = 9.14d-24 .* exp( 7801 ./ TTt3 ); %Dibble et al. (2012)
    k3    = 3.9d-11; %Balabanov et al. (2005)
    k4_no2 = 8.6d-11; %Wang et al. (2014)
    k4_ho2 = 8.2d-11; %Wang et al. (2014)
    k4_oh  = 6.33d-11; %Wang et al. (2014)
    k4_br  = 6.33d-11; %Wang et al. (2014)
    k4_bro = 1.09d-10; %Wang et al. (2014)
    
    % ------------------------- O3 reaction -------------------------
    % (R5) Hg0 + O3 --> HgIIO + O2                    unit: cm3 molecule-1 s-1
    
    k5max = 1.7d-18 .* ( Nair3 ./ N0 ); %Iverfeldt and Lindqvist (1986)
    k5min = 3d-20 .* ( Nair3 ./ N0 ); %Hall (1995)
    
    % ------------------------- OH reaction ------------------------- 
    % (R6) Hg0 + OH --> HgII                          unit: cm3 molecule-1 s-1
    
    k6max = 3.2d-13 .* ( ( TTt3 ./ 298d0 ) .^ (-3.06d0) ) .* ( Nair3 ./ N0 ); %Goodsite et al. (2004)
    k6min = 8.7d-14 .* ( Nair3 ./ N0 ); %Sommar et al. (2001)
    
    switch LBROX
        case 0
            kbrox = zeros(size(TTt3));
        case 1
            tmp = k4_no2.*moNO2c + k4_ho2.*moHO2c + k4_oh.*moOHc + k4_br.*moBrc + k4_bro.*moBrOc;
            kbrox = k1min .* moBrc .* tmp ./ ( k1min ./ keq12 + k3 .* moBrc + tmp );
        case 2
            tmp = k4_no2.*moNO2c + k4_ho2.*moHO2c + k4_oh.*moOHc + k4_br.*moBrc + k4_bro.*moBrOc;
            kbrox = k1max .* moBrc .* tmp ./ ( k1max ./ keq12 + k3 .* moBrc + tmp );
    end
    
    switch LO3OX
        case 0
            ko3ox = zeros(size(TTt3));
        case 1
            ko3ox = k5min .* moO3c;
        case 2
            ko3ox = k5max .* moO3c;
    end
    
    switch LOHOX
        case 0
            kohox = zeros(size(TTt3));
        case 1
            kohox = k6min .* moOHc;
        case 2
            kohox = k6max .* moOHc;
    end
    
    % Total OX Rates [s-1]
    ktox = double( kbrox + ko3ox + kohox );
    ktoxm = fcn_reshape1( ktox, Kzem, XXm, XXh ); %reshape array
    ktoxm(ktoxm < 0) = 0d0; %avoid negative values in interpolation
    
%% processes in surface snowpack
    
    %----------------------------------------------------------------------
    % photoreduction of HgII in surface snowpack
    %----------------------------------------------------------------------
    
    jo1d = DSWtm .^ 1.8; %based on J(O1D), Kukui et al. 2014 ACP and GEOS-Chem
    ksn1m = ((TAU1*86400).^-1) .* jo1d ./ mean(mean(jo1d(:,1:120))); %red rates in snow [s-1]
    
    %----------------------------------------------------------------------
    % Dark reduction of HgII in surface snowpack
    %----------------------------------------------------------------------
    
    ksn2m = zeros( size(DSWtm) ) + (TAU2*86400).^-1; %red rates in snow [s-1]
    
    % scale reaction rates by NOx levels
    m1=[0, 31,59, 90,120,151,181,212,243,273,304,334]+15;
    tmp1=interp1(m1,(mean(moNO(:,:,1),2)+mean(moNO2(:,:,1),2))./43,(0.5:1:364.5),'linear','extrap');
    tmp2=[fliplr(tmp1(1:61)),tmp1]; %add Nov-Dec 2012
    for cnt1=1:1:size(ksn2m,1)
        ksn2m(cnt1,:)=ksn2m(cnt1,:).*tmp2;
    end
    clear tmp1 tmp2
    
    % apply to dark period
    ksn2m(DSWtm>=1e-4) = 0.;
    ksn2m(:,[1:(nday-365+53),(nday-365+296):end])=0.;
    
%% save input variables
    
    % plot monthly diel profiles
    junk = fcn_vplot( ZZt, moO3, moNO, moNO2, moOH, moHO2, moBr, moBrO, ktox, 0 );
    
    save( './data/input.mat', ...
          'ROU', 'EDPTH', 'FACSTP2LTP', 'HG0DVD', 'HG2DVD', 'DT', 'TS', 'BETA', 'DD', ...
          'Kzem', 'TKEtm', 'ktoxm', 'ksn1m', 'ksn2m','upHg0m', 'upHg2m', ...
          'nday', 'nlay', 'LMAR', 'LBROX', 'LBRMR', 'LOHOX', 'LOHMR', 'LO3OX', 'LO3MR', ...
          'TAU1', 'TAU2', 'upHg0', 'upHg2', 'moBr', 'kbrox', 'ko3ox', 'kohox', 'XXm', 'XXh', 'XX1' )
    
%% Prepare output variables: mass concentrations and flux
    
    svairHg0 = zeros([size(Kzem,1),nday,size(Kzem,3)]); %[ng/m3 STP]
    svairHg2 = zeros([size(Kzem,1),nday,size(Kzem,3)]); %[ng/m3 STP]
    svsnHg0  = zeros([size(Kzem,1),nday]); %[ng/m3 STP]
    svsnHg2m = zeros([size(Kzem,1),nday]); %[ng/m2 STP]
    svHg0dep = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    svHg2dep = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    svHg0exg = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    svHg2exg = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    svHg0sne = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    svHg2sn1 = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    svHg2sn2 = zeros([size(Kzem,1),nday]); %[ng/m2/TS STP]
    
%% Loop over all the days
    
    tStart=tic; fprintf('Loop begins: \n');
    
    for cnt=1:1:nday %loop starts
        
        boxin.Kz    = squeeze(Kzem(:,cnt,:)); %Kz [m2/s]
        boxin.TKE   = TKEtm(:,cnt); %TKE [m2/s2]
        boxin.ktox  = squeeze(ktoxm(:,cnt,:)); %tot ox rates [s-1]
        boxin.ksn1  = ksn1m(:,cnt); %snow photored rates [s-1]
        boxin.ksn2  = ksn2m(:,cnt); %snow darkred rates [s-1]
        boxin.upHg0 = upHg0m(:,cnt); %FT Hg0 conc [molec cm-3 LTP]
        boxin.upHg2 = upHg2m(:,cnt); %FT Hg2 conc [molec cm-3 LTP]
        
        % initial conditions
        if cnt == 1
            snHg0int = 0.0; %ng/m3 STP
            snHg2cint = 5.0; %initial Hg conc in surface snow, ng/L
            snHg2mint = snHg2cint * ROU * EDPTH; %initial total Hg mass in surface snow, ng/m2
            hg0int = 1.0; %ng/m3 STP
            hg2int = 0.0; %ng/m3 STP
            boxin.airHg0i = zeros(size(DD'))+hg0int; %[ng/m3 STP]
            boxin.airHg2i = zeros(size(DD'))+hg2int; %[ng/m3 STP]
            boxin.snHg0i  = snHg0int; %[ng/m3 STP]
            boxin.snHg2mi = snHg2mint; %[ng/m2 STP]
        else
            boxin.airHg0i = boxout.airHg0(end,:); %[ng/m3 STP]
            boxin.airHg2i = boxout.airHg2(end,:); %[ng/m3 STP]
            boxin.snHg0i  = boxout.snHg0(end); %[ng/m3 STP]
            boxin.snHg2mi = boxout.snHg2m(end); %[ng/m2 STP]
        end
        
        % solve continuity equations
        boxout = fcn1( boxin );
        
        fprintf('Finish day # %d\n', round(cnt));
        fprintf('Air surf  Hg0 conc: %d\n', mean(boxout.airHg0(:,1)));
        fprintf('Air surf  Hg2 conc: %d\n', mean(boxout.airHg2(:,1)));
        fprintf('Snow surf Hg0 conc: %d\n', mean(boxout.snHg0));
        fprintf('Snow surf Hg2 mass: %d\n', mean(boxout.snHg2m));
        
        % save variables
        svairHg0(:,cnt,:) = boxout.airHg0; %[ng/m3 STP]
        svairHg2(:,cnt,:) = boxout.airHg2; %[ng/m3 STP]
        svsnHg0(:,cnt)    = boxout.snHg0; %[ng/m3 STP]
        svsnHg2m(:,cnt)   = boxout.snHg2m; %[ng/m2 STP]
        svHg0dep(:,cnt)   = boxout.Hg0dep; %[ng/m2 STP]
        svHg2dep(:,cnt)   = boxout.Hg2dep; %[ng/m2 STP]
        svHg0exg(:,cnt)   = boxout.Hg0exg; %[ng/m2 STP]
        svHg2exg(:,cnt)   = boxout.Hg2exg; %[ng/m2 STP]
        svHg0sne(:,cnt)   = boxout.Hg0sne; %[ng/m2 STP]
        svHg2sn1(:,cnt)   = boxout.Hg2sn1; %[ng/m2 STP]
        svHg2sn2(:,cnt)   = boxout.Hg2sn2; %[ng/m2 STP]
        
    end %loop ends
    
    tElapsed = toc(tStart); 
    fprintf('Total Loop time in min %d\n', round(tElapsed/60));
    
%% Save model output variables
    
    save( './data/output.mat', ...
          'svairHg0', 'svairHg2', 'svsnHg0', 'svsnHg2m', 'svHg0dep', 'svHg2dep', ...
          'svHg0exg', 'svHg2exg', 'svHg0sne', 'svHg2sn1', 'svHg2sn2' )
    