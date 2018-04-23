function junk = fcn_vplot( zzt, moo3, mono, mono2, mooh, moho2, mobr, mobro, ktox, Lplot )
% 
% Shaojie Song, 03/30/2018
% This is a matlab routine used to plot the vertical profiles of chemical tracers, plot monthly, diurnal profiles
% 
%%
if Lplot == 1
    
    xx=0.5:1:23.5;
    moTLE = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
    m1 = [0, 31,59, 90,120,151,181,212,243,273,304,334];
    m2 = [31,59,90,120,151,181,212,243,273,304,334,365];
    nlay = size(zzt,2);
    zzt2 = reshape(zzt(end+1-87600:end,:),[10,8760,nlay]);
    zzt3 = reshape(squeeze(mean(zzt2,1)),[24,365,nlay]);
    
    % HO2, molec cm-3 STP
    figure(35)
    
    for cnt1=1:1:12
        
        subplot(4,3,cnt1)
        
        [X,Y] = meshgrid( xx, squeeze(mean(mean(zzt3(:,m1(cnt1)+1:m2(cnt1),:),1),2)) );
        
        surf(X,Y,squeeze(moho2(cnt1,:,:))','EdgeColor','None'); view(2)
        set(gca,'XLim',[0.4 23.6],'XTick',(4:4:20),'XTickLabel',{'4','8','12','16','20'})
        set(gca,'YLim',[2 400],'YTick',[2,10,100,400])
        set(gca,'YScale','log')
        set(gca,'TickDir','out')
        set(gca,'ticklength',2*get(gca,'ticklength'))
        
        if mod(cnt1,3) == 1
            ylabel('Height (m)')
        end
        if cnt1 >= 10
            xlabel('Local Time (h)')
        end
        
        caxis([0 2e8])
        if cnt1 == 12
            colorbar
            hcb=colorbar;
            set(get(hcb,'title'),'string','HO_2 level (cm^-^3)','Rotation',90.0);
        end
        title(moTLE(cnt1,:))
        
    end
    
    % BrO, pptv
    figure(36)
    
    for cnt1=1:1:12
        
        subplot(4,3,cnt1)
        
        [X,Y] = meshgrid( xx, squeeze(mean(mean(zzt3(:,m1(cnt1)+1:m2(cnt1),:),1),2)) );
        
        surf(X,Y,squeeze(mobro(cnt1,:,:))','EdgeColor','None'); view(2)
        set(gca,'XLim',[0.4 23.6],'XTick',(4:4:20),'XTickLabel',{'4','8','12','16','20'})
        set(gca,'YLim',[2 400],'YTick',[2,10,100,400])
        set(gca,'YScale','log')
        set(gca,'TickDir','out')
        set(gca,'ticklength',2*get(gca,'ticklength'))
        
        if mod(cnt1,3) == 1
            ylabel('Height (m)')
        end
        if cnt1 >= 10
            xlabel('Local Time (h)')
        end
        
        caxis([0 1])
        if cnt1 == 12
            colorbar
            hcb=colorbar;
            set(get(hcb,'title'),'string','BrO level (pptv)','Rotation',90.0);
        end
        title(moTLE(cnt1,:))
        
    end

    % ktox, s-1
    figure(39)
    
    for cnt1=1:1:12
        
        subplot(4,3,cnt1)
        
        [X,Y] = meshgrid( xx, squeeze(mean(mean(zzt3(:,m1(cnt1)+1:m2(cnt1),:),1),2)) );
        
        surf(X,Y,squeeze(ktox(cnt1,:,:))','EdgeColor','None'); view(2)
        set(gca,'XLim',[0.4 23.6],'XTick',(4:4:20),'XTickLabel',{'4','8','12','16','20'})
        set(gca,'YLim',[2 400],'YTick',[2,10,100,400])
        set(gca,'YScale','log')
        set(gca,'TickDir','out')
        set(gca,'ticklength',2*get(gca,'ticklength'))
        
        if mod(cnt1,3) == 1
            ylabel('Height (m)')
        end
        if cnt1 >= 10
            xlabel('Local Time (h)')
        end
        
        caxis([0 1d-5])
        if cnt1 == 12
            colorbar
            hcb=colorbar;
            set(get(hcb,'title'),'string','Hg^0 ox rate (s^-^1)','Rotation',90.0);
        end
        title(moTLE(cnt1,:))
        
    end
    
    % diel monthly profiles of all species
    xx2 = 0.5:1:12*24-0.5;
    zzt4 = zeros(24,12,nlay);
    for c1=1:1:24
        for c2=1:1:12
            zzt4(c1,c2,:) = mean(squeeze(zzt3(c1,m1(c2)+1:m2(c2),:)),1);
        end
    end
    zzt5 = reshape( zzt4, [12*24 nlay] );
    tmp1 = zeros( size(zzt5) );
    tmp1(:,2:nlay) = zzt5(:,1:nlay-1);
    zzt6 = tmp1./2 + zzt5./2; clear tmp1;
    
    % reshape
    [moo3_2, mooh_2, mobr_2, mono_2, mono2_2] = deal( zeros(24*12,nlay) );
    for c3=1:1:nlay
        for c2=1:1:24
            for c1=1:1:12
                mobr_2((c1-1)*24+c2,c3)=mobr(c1,c2,c3);
                moo3_2((c1-1)*24+c2,c3)=moo3(c1,c2,c3);
                mooh_2((c1-1)*24+c2,c3)=mooh(c1,c2,c3);
                mono_2((c1-1)*24+c2,c3)=mono(c1,c2,c3);
                mono2_2((c1-1)*24+c2,c3)=mono2(c1,c2,c3);
            end
        end
    end
    pmos = cat( 3, moo3_2, mooh_2, mobr_2 );
    pmos2 = cat( 3, mono_2, mono2_2 );
    
    [X,~] = meshgrid( xx2, mean(zzt6,1) );
    
    % X axis
    tmp1=(0:24:12*23)+4;
    tmp2=(0:24:12*23)+12;
    tmp3=(0:24:12*23)+20;
    tmp4=reshape(cat(1,tmp1,tmp2,tmp3),[12*3 1]);
    clear tmp1 tmp2 tmp3
    tmp1=zeros(3,12); tmp1(1,:)=4; tmp1(2,:)=12; tmp1(3,:)=20;
    tmp2=num2cell(reshape(tmp1,[12*3 1]));
    clear tmp1
    
    figure(37)
    
    for c1=1:1:3
        
        subplot(3,1,c1)
        surf(X,zzt6',pmos(:,:,c1)','EdgeColor','None'); view(2)
        hold on
        set(gca,'YLim',[1 500],'YTick',[1,10,100,500])
        set(gca,'YScale','log')
        set(gca,'TickDir','out')
        set(gca,'ticklength',1.0*get(gca,'ticklength'))
        ylabel('Height (m)')
        
        if c1==3
            xlabel('Local time (h)')
        end
        
        hcb = colorbar;
        set(hcb,'ticklength',0.1);
        set(hcb,'fontsize',11);
        
        for c2=1:1:13
            plot3([24*(c2-1) 24*(c2-1)],[1 500],[1e7 1e7],'color','k','linewidth',1.5)
        end
        plot3([0 288],[1 1],[1e7 1e7],'color','k','linewidth',1.5)
        plot3([0 288],[500 500],[1e7 1e7],'color','k','linewidth',1.5)
        
        set(gca,'XLim',[0 288],'XTick',tmp4-0.5,'XTickLabel',tmp2)
        
        switch c1
            case 1
                set(get(hcb,'title'),'string','O_3 (ppbv)','Rotation',90.0);
                caxis([0 40])
            case 2
                set(get(hcb,'title'),'string','OH (molecule cm^-^3)','Rotation',90.0);
                caxis([0 8e6])
            case 3
                set(get(hcb,'title'),'string','Br (pptv)','Rotation',90.0);
                caxis([0 0.25])
        end
        
        if c1==1
            for c2=1:1:12
                t = text(24*(c2-1)+12-0.5,1000,moTLE(c2,:),'HorizontalAlignment','center');
%                 s = t.FontSize;
                t.FontSize = 16;
            end
        end
        
    end
    
    figure(38)
    
    for c1=1:1:2
        
        subplot(2,1,c1)
        surf(X,zzt6',pmos2(:,:,c1)','EdgeColor','None'); view(2)
        hold on
        set(gca,'YLim',[1 500],'YTick',[1,10,100,500])
        set(gca,'YScale','log')
        set(gca,'TickDir','out')
        set(gca,'ticklength',1.0*get(gca,'ticklength'))
        ylabel('Height (m)')
        
        if c1==2
            xlabel('Local time (h)')
        end
        
        hcb = colorbar;
        set(hcb,'ticklength',0.1);
        set(hcb,'fontsize',11);
        
        for c2=1:1:13
            plot3([24*(c2-1) 24*(c2-1)],[1 500],[1e7 1e7],'color','k','linewidth',1.5)
        end
        plot3([0 288],[1 1],[1e7 1e7],'color','k','linewidth',1.5)
        plot3([0 288],[500 500],[1e7 1e7],'color','k','linewidth',1.5)
        
        set(gca,'XLim',[0 288],'XTick',tmp4-0.5,'XTickLabel',tmp2)
        
        switch c1
            case 1
                set(get(hcb,'title'),'string','NO (pptv)','Rotation',90.0);
                caxis([0 250])
            case 2
                set(get(hcb,'title'),'string','NO_2 (pptv)','Rotation',90.0);
                caxis([0 250])
        end
        
        if c1==1
            for c2=1:1:12
                t = text(24*(c2-1)+12-0.5,1000,moTLE(c2,:),'HorizontalAlignment','center');
%                 s = t.FontSize;
                t.FontSize = 16;
            end
        end
        
    end
    
end

junk = 1;