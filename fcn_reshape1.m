function newdata = fcn_reshape1( olddata, kzem, xxm, xxh )
% 
% Shaojie Song, 03/30/2018
% This is a matlab routine used to reshape monthly to minutely data
% 
%% 
    
    % initialize
    nday = size( kzem, 2 );
    nlay = size( kzem, 3 );
    tmp1 = zeros( nday, 24, nlay );
    
    for cnt1=1:1:30 %2012 Nov
        tmp1(cnt1,:,:) = olddata(11,:,:);
    end
    
    for cnt1=31:1:61 %2012 Dec
        tmp1(cnt1,:,:) = olddata(12,:,:);
    end
    
    for cnt1=62:1:92 %2013 Jan
        tmp1(cnt1,:,:) = olddata(1,:,:);
    end
    
    for cnt1=93:1:120 %2013 Feb
        tmp1(cnt1,:,:) = olddata(2,:,:);
    end
    
    for cnt1=121:1:151 %2013 Mar
        tmp1(cnt1,:,:) = olddata(3,:,:);
    end
    
    for cnt1=152:1:181 %2013 Apr
        tmp1(cnt1,:,:) = olddata(4,:,:);
    end
    
    for cnt1=182:1:212 %2013 May
        tmp1(cnt1,:,:) = olddata(5,:,:);
    end
    
    for cnt1=213:1:242 %2013 Jun
        tmp1(cnt1,:,:) = olddata(6,:,:);
    end
    
    for cnt1=243:1:273 %2013 Jul
        tmp1(cnt1,:,:) = olddata(7,:,:);
    end
    
    for cnt1=274:1:304 %2013 Aug
        tmp1(cnt1,:,:) = olddata(8,:,:);
    end
    
    for cnt1=305:1:334 %2013 Sep
        tmp1(cnt1,:,:) = olddata(9,:,:);
    end
    
    for cnt1=335:1:365 %2013 Oct
        tmp1(cnt1,:,:) = olddata(10,:,:);
    end
    
    for cnt1=366:1:395 %2013 Nov
        tmp1(cnt1,:,:) = olddata(11,:,:);
    end
    
    for cnt1=396:1:426 %2013 Dec
        tmp1(cnt1,:,:) = olddata(12,:,:);
    end
    
%%
    
    % initialize
    tmp2 = zeros( size( kzem ) );
    
    for cnt1=1:1:nday
        for cnt2=1:1:nlay
            tmp3 = squeeze(tmp1(cnt1,:,cnt2));
            tmp4 = interp1(xxh,tmp3,xxm,'linear','extrap');
            tmp2(:,cnt1,cnt2) = tmp4;
            clear tmp3 tmp4
        end
    end
    
    newdata = tmp2;