function boxout = fcn1( boxin )
% 
% Shaojie SONG, 06/29/2016, created
% Shaojie SONG, 03/30/2018, modified
% 
    %----------------------------------------------------------------------
    % Initialize
    
    global ROU EDPTH FACSTP2LTP DT TS DD
    
    airHg0mr = zeros( size(boxin.Kz) ); %air Hg mass [m * molec cm-3 LTP]
    airHg2mr = zeros( size(boxin.Kz) ); %air Hg mass [m * molec cm-3 LTP]
    snHg0mr  = zeros( size(boxin.Kz,1), 1 ); %snow Hg0 mass [m * molec cm-3 LTP]
    snHg2mmr = zeros( size(boxin.Kz,1), 1 ); %snow Hg2 mass [m * molec cm-3 LTP]
    Hg0exgr  = zeros( size(boxin.Kz,1), 1 ); %flux [m * molec cm-3 LTP per time step]
    Hg2exgr  = zeros( size(boxin.Kz,1), 1 ); %below the same
    Hg0dfdr  = zeros( size(boxin.Kz,1), 1 );
    Hg2dfdr  = zeros( size(boxin.Kz,1), 1 );
    Hg0sner  = zeros( size(boxin.Kz,1), 1 );
    Hg2sn1r  = zeros( size(boxin.Kz,1), 1 );
    Hg2sn2r  = zeros( size(boxin.Kz,1), 1 );
    
    %----------------------------------------------------------------------
    % Loop over all time steps
    
    for cnt = 1:1:numel(TS) %loop begins
        
        if cnt == 1
            y0 = [ boxin.airHg0i' .* FACSTP2LTP .* DD; ...
                   boxin.airHg2i' .* FACSTP2LTP .* DD; ...
                   boxin.snHg0i   .* FACSTP2LTP .* EDPTH .* ( 1 - ROU./900 ); ...
                   boxin.snHg2mi  .* FACSTP2LTP; ...
                   0; 0; 0; 0; 0; 0; 0 ];
        else
            y0 = [ airHg0mr(cnt-1,:)'; airHg2mr(cnt-1,:)'; snHg0mr(cnt-1); snHg2mmr(cnt-1); ...
                   0; 0; 0; 0; 0; 0; 0 ];
        end
        
        if max(boxin.Kz(cnt,:)) < 1
            [~,y] = ode23( @(t,y) fcn2( t,y, ...
                    boxin.Kz(cnt,:), boxin.ktox(cnt,:), boxin.upHg0(cnt), boxin.upHg2(cnt), ...
                    boxin.ksn1(cnt), boxin.TKE(cnt), boxin.ksn2(cnt) ), [0 DT*60], y0 );
        else
            [~,y] = ode15s( @(t,y) fcn2( t,y, ...
                    boxin.Kz(cnt,:), boxin.ktox(cnt,:), boxin.upHg0(cnt), boxin.upHg2(cnt), ...
                    boxin.ksn1(cnt), boxin.TKE(cnt), boxin.ksn2(cnt) ), [0 DT*60], y0 );
        end
        
        airHg0mr(cnt,:) = y(end, 1:numel(DD)); %mass [m * molec cm-3 LTP]
        airHg2mr(cnt,:) = y(end, numel(DD)+1:numel(DD)*2);
        snHg0mr(cnt)    = y(end, numel(DD)*2+1);
        snHg2mmr(cnt)   = y(end, numel(DD)*2+2);
        Hg0exgr(cnt)    = y(end, numel(DD)*2+3); %flux [m * molec cm-3 LTP per time step]
        Hg2exgr(cnt)    = y(end, numel(DD)*2+4);
        Hg0dfdr(cnt)    = y(end, numel(DD)*2+5);
        Hg2dfdr(cnt)    = y(end, numel(DD)*2+6);
        Hg0sner(cnt)    = y(end, numel(DD)*2+7);
        Hg2sn1r(cnt)    = y(end, numel(DD)*2+8);
        Hg2sn2r(cnt)    = y(end, numel(DD)*2+9);
        
        if mod(cnt,240/DT) == 0
            fprintf('Finish hours # %d\n', cnt*DT/60);
        end
        
        clear y0
        
    end %loop ends
    
    %----------------------------------------------------------------------
    % Save output variables
    
    [c,nlay]=size(airHg0mr); tmp1=repmat(DD,c,1); tmp2=reshape(tmp1,[nlay,c])';
    
    boxout.airHg0 = airHg0mr ./ FACSTP2LTP ./ tmp2; %air Hg conc [ng/m3 STP]
    boxout.airHg2 = airHg2mr ./ FACSTP2LTP ./ tmp2; %air Hg conc [ng/m3 STP]
    boxout.snHg0  = snHg0mr  ./ FACSTP2LTP ./ EDPTH ./ ( 1 - ROU./900 ); %snow Hg conc [ng/m3 STP]
    boxout.snHg2m = snHg2mmr ./ FACSTP2LTP; %snow Hg mass [ng/m2 STP]
    boxout.Hg0dep = Hg0dfdr  ./ FACSTP2LTP; %ng m-2 STP per time step
    boxout.Hg2dep = Hg2dfdr  ./ FACSTP2LTP; %ng m-2 STP per time step
    boxout.Hg0exg = Hg0exgr  ./ FACSTP2LTP; %ng m-2 STP per time step
    boxout.Hg2exg = Hg2exgr  ./ FACSTP2LTP; %ng m-2 STP per time step
    boxout.Hg0sne = Hg0sner  ./ FACSTP2LTP; %ng m-2 STP per time step
    boxout.Hg2sn1 = Hg2sn1r  ./ FACSTP2LTP; %ng m-2 STP per time step
    boxout.Hg2sn2 = Hg2sn2r  ./ FACSTP2LTP; %ng m-2 STP per time step
    
    clear airHg0mr airHg2mr snHg0mr snHg2mmr Hg0dfdr Hg2dfdr Hg0exgr Hg2exgr Hg0sner Hg2sn1r Hg2sn2r
    clear c nlay tmp1 tmp2 cnt