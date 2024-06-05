%% Estes el codgo bueno de INTERPLAY !!!%%
clear
close all
clc

read_EDGAR=0;
read_hysplit=0;
INTERPLAY_avg=0;
%%
% Tests: w/o rain or PBL, integrating in-situ BC full hour after BT R2 is
% 0.17 - count 19932

% Tests: w/o rain or PBL, integrating in-situ BC half hour before
%and after BT R2 is 0.168 - count 19935

% Tests: w/o PBL, adding rain condition (1mm), integrating in-situ BC half hour before
% and after BT R2 is 0.173 - count 14189

% Tests: w/o PBL, adding rain (1mm) condition , integrating
%in-situ BC full hour after BT BT R2 is 0.177 - count 14191

% Tests: w/o rain, adding PBL condition (+500m), integrating 
%in-situ BC full hour after BT BT R2 is 0.151 - n_points=376000

% Tests: adding rain (1mm) and PBL condition (+500m), integrating 
%in-situ BC full hour after BT BT R2 is 0.153 - n_points=376000

% Tests: w/o rain, adding PBL condition (+2000m), integrating 
%in-situ BC full hour after BT BT R2 is 0.167 - n_points=48000

% Tests: adding rain (1mm) and PBL condition (+2000m), integrating 
%in-situ BC full hour after BT BT R2 is 0.172 - n_points=48000

% Tests: w/o rain, adding PBL condition (+3000m), integrating 
%in-situ BC full hour after BT BT R2 is 0.168 - n_points=12000

% Tests: adding rain (1mm) and PBL condition (+3000m), integrating 
%in-situ BC full hour after BT BT R2 is 0.175 - n_points=12000

% Tests: w/o rain, adding PBL condition (+3000m), integrating 
%in-situ BC full hour after BT BT R2 is 0.17 - n_points=890

% Tests: adding rain (1mm) and PBL condition (+3000m), integrating 
%in-situ BC full hour after BT BT R2 is 0.176 - n_points=890

% Tests: w/o rain or PBL, integrating in-situ BC full hour after BT and 
%interpolating to 10min R2 is % 0.169 - count 19932

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT R2 is 0.176 - count 14191

% Tests: w/o PBL, adding rain condition (0.5), integrating 
%in-situ BC full hour after BT R2 is 0.180 - count 10903

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (360°) R2 is 0.176 - count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (270°) R2 is 0.171 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°) R2 is 0.165 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (90°) R2 is 0.155 - count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (60°) R2 is 0.147 - count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (30°) R2 is 0.139 - count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°), using 0.4°x0.4°
% the R2 is 0.165 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°), using 0.3°x0.3°
% the R2 is 0.164 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°), using 0.2°x0.2°
% the R2 is 0.161 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°), using 0.6°x0.6°
% the R2 is 0.165 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°), using 0.8°x0.8°
% the R2 is 0.164 count 14191

% Tests: w/o PBL, adding rain condition (1mm), integrating 
%in-situ BC full hour after BT, adding Lille condition (180°), using 1°x1°
% the R2 is 0.164 count 14191

%% %%%%%%%%%%%%%%%%%%%%%% INTERPLAY approach %%%%%%%%%%%%%%%%%%

%%% BC_inv=ncread('3_EDGAR\v50_BC_2015.0.1x0.1.nc','emi_bc');

A=122.98e6; %The area in m2 relative to each grid point

kg_2_Gg=1e-6;
kg_2_mg=1e6;
s_2_yr=3600.*24.*365;

sector={'AWB';'ENE';'IND';'RCO';'REF_TRF';'TNR_Ship';'TRO_noRES';'TRO_RES'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. IMPORT EDGAR inventory %%%%%%%%%%%%%%%%%%%%
if read_EDGAR==1
    
    for j=1:12
        for i=1:max(size(sector))
            fname=dir(['3_EDGAR\2018\*',num2str(j),'_',sector{i},'*.nc']);
            dummy_name=['3_EDGAR\2018\',fname(1).name];
            eval(['BC_',sector{i},'_',num2str(j),'=ncread(dummy_name,''emi_bc'');'])
            if i==1
                lat_inv=ncread(dummy_name,'lat');
                lon_inv=ncread(dummy_name,'lon')-180;
                eval(['BC_sum_',num2str(j),'=BC_',sector{i},'_',num2str(j),'.*0;'])
            end
            eval(['BC_',sector{i},'_',num2str(j),'=ncread(dummy_name,''emi_bc'').*A.*kg_2_Gg.*s_2_yr;'])
            eval(['dummy=BC_',sector{i},'_',num2str(j),';'])
            eval(['BC_',sector{i},'_',num2str(j),'(1:1800,:)=dummy(1801:end,:);'])
            eval(['BC_',sector{i},'_',num2str(j),'(1801:end,:)=dummy(1:1800,:);'])
            eval(['BC_sum_',num2str(j),'=BC_sum_',num2str(j),'+BC_',sector{i},'_',num2str(j),';'])
        end
        if j==1
            eval(['BC_sum=BC_sum_',num2str(j),';'])
        else
            eval(['BC_sum=BC_sum+BC_sum_',num2str(j),';'])
        end
    end
    save 3_EDGAR\2018_sectors.mat sector BC_* lat_inv lon_inv
else
    load 3_EDGAR\2018_sectors.mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. BC_GP - Hysplit+EDGAR %%%%%%%%%%%%%%%%%%%%%% 
if read_hysplit==1
    count=0;
    for i=datenum(2016,12,10):datenum(2020,01,01)
  %  for i=datenum(2016,10,01):datenum(2016,12,01)
        d0=day(i);
        m0=month(i);
        y0=year(i);
        y0_str=num2str(y0);
        y0_str=y0_str(3:4);
        
        for h=0:23
            [dummy_lat(:,h+1),dummy_lon(:,h+1),dummy_alt(:,h+1),dummy_rain(:,h+1),dummy_PBL(:,h+1)] = import_hysplit(['3_Hysplit\Lille\tdumpLIL',y0_str,num2str(m0,'%02d'),num2str(d0,'%02d'),num2str(h,'%02d')]);
            if max(isnan(dummy_lat(:,h+1)))==0
                count=count+1;
                time_BT(count)=i+h/24;
                lat(:,count)=dummy_lat(:,h+1);
                lon(:,count)=dummy_lon(:,h+1);
                alt(:,count)=dummy_alt(:,h+1);
                rain(:,count)=dummy_rain(:,h+1);
                PBL(:,count)=dummy_PBL(:,h+1);
                
              %  idx=find(height(:,count)<(PBL(:,count)+500));
                
                %%%%% 1. Angle between Lille & 2nd point BT (Local contribution)%%%%%
                %test_east=0;
                Li_vector(1,count)=atan2((lat(4,count)-lat(1,count)),(lon(4,count)-lon(1,count))).*180./pi;%before 4th point in 1hrs time
                if Li_vector(1,count)<0
                    Li_vector(1,count)=Li_vector(1,count)+360;
                end
                                
                
            end
        end
    end
    
    save 1_mat_files\Hysplit.mat lat lon time_BT rain alt PBL Li_vector
else
    load 1_mat_files\Hysplit.mat
end

load 1_mat_files\BC_AE33_ATOLL_flag.mat

lat=interp1(1:73,lat(:,1:end),1:1./6:73);
lon=interp1(1:73,lon(:,1:end),1:1./6:73);
alt=interp1(1:73,alt(:,1:end),1:1./6:73);
PBL=interp1(1:73,PBL(:,1:end),1:1./6:73);
rain=interp1(1:73,rain(:,1:end),1:1./6:73);

%deg=[0 5 30 60 90 120 150 175];
% deg=175;
% for i=1:max(size(deg))


d_Sec = 180;


d_Sec=d_Sec/2-d_Sec:5:d_Sec-d_Sec/2;

% Y_L=lat(1,1)+2.*sind(Li_vector-d_Sec/2);
% X_L=lon(1,1)+2.*cosd(Li_vector-d_Sec/2);
% Y_H=lat(1,1)+2.*sind(Li_vector+d_Sec/2);
% X_H=lon(1,1)+2.*cosd(Li_vector+d_Sec/2);



% figure
% geoscatter(Y_L(1,1),X_L(1,1))
% hold on
% geoscatter(Y_H(1,1),X_H(1,1))
% geoscatter(lat(1:20,1),lon(1:20,1))
% geoscatter(lat(1,1)-0.25,lon(1,1)-0.25)
% geoscatter(lat(1,1)-0.25,lon(1,1)+0.25)
% geoscatter(lat(1,1)+0.25,lon(1,1)-0.25)
% geoscatter(lat(1,1)+0.25,lon(1,1)+0.25)
% title(num2str(deg(i)))
% end

%[dummy_idx,~]=inpolygon(X,Y,[lon_lille X_L X_H lon_lille],[lat_lille Y_L Y_H lat_lille]);
                            
n_points=0;
count_EBC=0;
size_rec=0.5;% km?

BC_map=BC_sum.*0;
for k=1:max(size(sector))
    eval(['BC_map_',sector{k},' = BC_map;'])
    for i=1:12
        eval(['BC_map_',sector{k},'_',num2str(i),' = BC_map;'])
    end
end

if INTERPLAY_avg ==1
for i=1:max(size(time_BT))
% i=1;
% while count_EBC<1000
    %idx=find(time_AE33>=time_BT(i)-0.5./24&time_AE33<time_BT(i)+0.5./24&idx_flag==0);
    idx=find(time_AE33>=time_BT(i)&time_AE33<time_BT(i)+1./24&idx_flag==0);
    m0_str=num2str(month(time_BT(i)));
    
    if min(size(idx))>0
        EBC_corr(i,:)=nanmean(EBC(idx,:));
                
%         count_EBC=count_EBC+1;
        for j=1:size(lat,1)
            idx_lon=find(lon(j,i)>=(lon_inv-size_rec./2)&lon(j,i)<=(lon_inv+size_rec./2));
            idx_lat=find(lat(j,i)>=(lat_inv-size_rec./2)&lat(j,i)<=(lat_inv+size_rec./2));
            dist_point(j,i)=pos2dist(lat(1,1),lon(1,1),lat(j,i),lon(j,i),2);
            
            dummy_idx=ones(size(idx_lon,1),size(idx_lat,1));
            if dist_point(j,i)<40 %this is the distance where the rectangle can include downtown lille
                [X,Y] = meshgrid(lon_inv(idx_lon),lat_inv(idx_lat));
                lat_circle=lat(1,1)+2.*sind(Li_vector(i)+d_Sec);
                lon_circle=lon(1,1)+2.*cosd(Li_vector(i)+d_Sec);%200km selected
                [dummy_idx,~]=inpolygon(X,Y,[lon(1,1) lon_circle lon(1,1)],[lat(1,1) lat_circle lat(1,1)]);%0.5
                dummy_idx=dummy_idx';
            end
            eval(['BC_GP(j,i)=sum(sum(BC_sum_',m0_str,'(idx_lon,idx_lat).*dummy_idx));'])
            count_filter(j,i)=nansum(nansum(dummy_idx))./(size(dummy_idx,1).*size(dummy_idx,2));
            
            for k=1:max(size(sector))
                eval(['BC_GP_',sector{k},'(j,i)=sum(sum(BC_',sector{k},'_',m0_str,'(idx_lon,idx_lat).*dummy_idx));'])
            end
        
            %%%%%%%%% MAPS contribution %%%%%%%%%%
            eval(['BC_map(idx_lon,idx_lat)=BC_map(idx_lon,idx_lat)+ BC_sum_',m0_str,'(idx_lon,idx_lat).*dummy_idx;'])
            
            for k=1:max(size(sector))
                eval(['BC_map_',sector{k},'_',m0_str,'(idx_lon,idx_lat) = BC_map_',sector{k},'_',m0_str,'(idx_lon,idx_lat) + BC_',sector{k},'_',m0_str,'(idx_lon,idx_lat).*dummy_idx;'])
                eval(['BC_map_',sector{k},'(idx_lon,idx_lat) = BC_map_',sector{k},'(idx_lon,idx_lat) + BC_',sector{k},'_',m0_str,'(idx_lon,idx_lat).*dummy_idx;'])
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    else
        EBC_corr(i,:)=ones(7,1).*NaN;
        BC_GP(1:size(lat,1),i)=NaN;
        count_filter(1:size(lat,1),i)=NaN;
        
        for k=1:max(size(sector))
            eval(['BC_GP_',sector{k},'(1:size(lat,1),i)=NaN;'])
        end
        
    end
    %      BC_GP(:,i)=BC_GP(:,i).*(alt(:,i)<PBL(:,i)+5000);
    %      n_points=n_points+size(BC_GP(alt(:,i)>=PBL(:,i)+5000,i),1);
%     i=i+1;
end
    save 1_mat_files\INTERPLAY_BC.mat EBC_corr BC_GP BC_map BC_GP_* BC_map_* dist_point count_filter%
else
    load 1_mat_files\INTERPLAY_BC.mat
end

idx=find(~isnan(EBC_corr(:,6)));
size(idx)
LM=fitlm(EBC_corr(:,6),nansum(BC_GP))

idx=intersect(idx,find(nanmax(rain)<1));
size(idx)
LM=fitlm(EBC_corr(idx,6),nansum(BC_GP(:,idx)))
nanmean(count_filter(:,1))
 
% rain_val=0:0.1:5;
% for i=1:max(size(rain_val))
%     idx=find(nanmax(rain)<=rain_val(i));
%     LM=fitlm(EBC_corr(idx,6),nansum(BC_GP(:,idx)));
%     r2(i)=LM.Rsquared.Adjusted;
%     n_points(i)=size(intersect(idx,find(~isnan(EBC_corr(:,6)))),1);
% end
% 
% figure
% circ(rain_val,r2)
% yyaxis right
% circ(rain_val,n_points,'r')

%% conversion from hysplit lat/lon to edgar grid (0.1°)

clear fig1;
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,23,8.2])
axes1 = axes('Parent',fig1,...
    'Position',[0.0804889298892989 0.115167958656331 0.523062730627307 0.815]);
hold(axes1,'on');box on
%hold on;box on;
colororder(axes1,{'r','b'})
yyaxis right
circ(time_BT,sum(BC_GP))
yyaxis left
plot(time_BT,EBC_corr(:,6))
ylim([0 6])
dynamicDateTicks([],[],'dd-mm')

% 
% load 1_mat_files\BC_AE33_ATOLL_flag.mat
% plot(time_AE33,EBC(:,6))
% ylim([0 6])

figure
loglog(nansum(BC_GP),EBC_corr(:,6),'*')
grid on

for i=1:12
   idx=find(month(time_BT)==i);
    idx=intersect(idx,find(nanmax(rain)<1));
    size(idx);
    LM=fitlm(EBC_corr(idx,6),nansum(BC_GP(:,idx)));
    r2(i)=LM.Rsquared.Adjusted;
    n_points(i)=size(intersect(idx,find(~isnan(EBC_corr(:,6)))),1); 
end

%scatter(lon(:,1000),lat(:,1000),10,BC_corr(:,1000),'filled')


return
if plot_EDGAR==1
    
    latlim = [30 60];
    lonlim = [-30 50];
    [X,Y] = meshgrid(lon_inv(idx_lon),lat_inv(idx_lat));
    BC_sum(BC_sum<nanmean(prctile(BC_sum,15)))=NaN;
    figure
     worldmap(latlim,lonlim)
    geoshow(Y, X, (BC_sum(idx_lon,idx_lat).*dummy_idx)', 'DisplayType', 'surface')%pcolor(X,Y,BC_sum(idx_lon,idx_lat)')
%     shading interp
    caxis([nanmean(prctile(BC_sum,15)) max(prctile(BC_sum,95))])
    title('BC')
    colormap('hot')
    colormap(flipud(colormap))
    colorbar
    load coastlines
    plotm(coastlat,coastlon,'Color','k')
end  


