%% Countribution of regions & sectors to BC %%
clear all; close all; clc;

%%% 1. EDGAR inventory%%%
%load 3_EDGAR\2018_sectors.mat 

%%% 2. Hysplit %%%%
load 1_mat_files\Hysplit.mat %%%% Hysplit

lat=interp1(1:73,lat(:,1:end),1:1./6:73);
lon=interp1(1:73,lon(:,1:end),1:1./6:73);
alt=interp1(1:73,alt(:,1:end),1:1./6:73);
PBL=interp1(1:73,PBL(:,1:end),1:1./6:73);
rain=interp1(1:73,rain(:,1:end),1:1./6:73);


%%%% 3. INTERPLAY %%%
load 1_mat_files\INTERPLAY_BC.mat %%%% INTERPLAY BC

%%---
idx=find(nanmax(rain)>=1);
BC_GP(:,idx)=NaN;
size(find(isnan(sum(BC_GP))==0))%%14191
fitlm(EBC_corr(:,6)',sum(BC_GP))%0.165
%%---

%% %%%%%%%%%%%%%%%%%%%%%%% Calculation of contribution by countries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

read_BC_countries=0;

if read_BC_countries==1
    
    INTERPLAY_countries_calc_Seasons
    
    else
    load 1_mat_files\BC_countries_Cold.mat
    %load 1_mat_files\BC_countries_Warm.mat
    
end

%% %%%%%%%%%%% Seasons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cold season vs hot season %%
idx=find(month(time_BT)>=10|month(time_BT)<=3);%cold periods
%idx=find(month(time_BT)>=4&month(time_BT)<9);%warm periods

BC_GP = BC_GP(:,idx);
EBC_corr = EBC_corr(idx,:);
time_BT = time_BT(idx);
rain = rain(:,idx);
BC_countries = [BC_countries(:,1:8)];%%remove the rest of the countries for analysis


BC_GP_TNR_Ship_DJF = BC_GP_TNR_Ship(:,idx);
BC_countries = [BC_countries (sum(BC_GP_TNR_Ship_DJF',2))];

% BC_GP_TNR_Ship_JJA = BC_GP_TNR_Ship(:,idx);
% BC_countries = [BC_countries (sum(BC_GP_TNR_Ship_JJA',2))];




%% %%%%% Ploting pie chart countries %%%%%%%%%%%%

clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,10,13])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'Fontname','Garamond',...
    'FontSize',18);
box on;
b2=pie([nanmean(BC_countries)]);% nanmean(sum(BC_SHIPPING,2))]);
lgd=legend('Lille','Paris','France',...
'Belgium','Netherland','UK','Germany','Poland','Rest of land','Non-land (Shipping)');
set(lgd,...
    'Position',[0.102130113565601 0.0159307420761397 0.852754251179049 0.175941084098504],...
    'NumColumns',2,...
    'FontSize',10,...
    'FontName','Times New Roman',...
    'FontWeight','normal');
colormap(jet)
title('Relative Contribution to BC')
set(gca,'FontSize',12, 'FontWeight','Demi','FontName','times','LineWidth',1.0);
set(findobj(b2,'type','text'),'FontWeight','Demi','FontName','times','FontSize',12)
%export_fig RC_test1(12H_baseline) -png -r300 -transparent

% figure;
% %circ(time_BT,nansum(BC_sectors))%%BC_corr
% circ(time_BT,BC_GP)%%BC_corr
% ylabel('EDGAR_{inv}+HYSPLIT','FontSize',18, 'FontWeight','Demi');
% yyaxis right
% plot(time_BT,BC6,'r')
% ylabel('In situ','FontSize',18, 'FontWeight','Demi')
% ylim([0 10])
% set(gca,'Fontname','Garamond','Fontsize',15,'linewidth',2);
% %dynamicDateTicks
% %xlim([datenum(2018,12,09,07,00,00) datenum(2019,01,07,13,00,00)])
% xlabel('Time')
% %dynamicDateTicks([],[],'dd-mm')
% grid on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Baseline (INTERPLAY vs In situ)%%%%%%%%%
Countries_baseline_test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% BC wb BC ff & AAE analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend_var = {'AAE_{370-880nm}','BC wb [µg/m^{-3}]','BC ff [µg/m^{-3}]'};%,...
%     'Abs BrC_{470nm} [Mm/m^{-1}]','SSA_{450nm}','SSA_{525nm}','SSA_{635nm}',...
%     'SAE_{450-635nm}','GMD [nm/m^{-3}]','Scat_{525nm}','Abs_{520nm}'};
% Country = {'Lille','Paris','France',...
%     'Belgium','Netherland','UK','Germany','Poland','Shipping'};
Country = {'Lille','Paris','FR',...
    'BE','NL','UK','DE','PL','Shipping'};
var={'AAE','BCWB','BCFF'};

Absorption_calc

[BC_ff,BC_wb,alpha] = sandradewi(EBC_avg(:,1),EBC_avg(:,2),EBC_avg(:,3),EBC_avg(:,4),EBC_avg(:,5),EBC_avg(:,6),EBC_avg(:,7));

%[BC_ff,BC_wb,alpha] = sandradewi(EBC_corr(:,1),EBC_corr(:,2),EBC_corr(:,3),EBC_corr(:,4),EBC_corr(:,5),EBC_corr(:,6),EBC_corr(:,7));


for j=1:max(size(var))
    for i=1:max(size(Country))
    MAT_countries.(var{j}).(Country{i}) = [];
    end
end

BC_INTERPLAY = sum(BC_GP_avg',2);

for j=1:9%10→Remove rest of Europe
    
    %%-- BTs for the total BC → updated: BC_GP(i,:) ---%%%%
    
    idx=find(BC_countries_avg(:,j)./BC_INTERPLAY>0.2);%0.2
    
    size_AOPs_idx(:,j)=(size(idx))
    
    %%%% AAE %%%
    LogAbsorption = log10(abs(Abs_BC1(idx))./Abs_BC6(idx));
    Wavelenght = log10(370/880);
    AAE_370_880 = -(LogAbsorption/Wavelenght);
    
    MAT_countries.AAE.(Country{j}) = AAE_370_880;
    
    MAT.AAE.median(:,j)=nanmedian(AAE_370_880);
    MAT.AAE.prc25(:,j)=prctile(AAE_370_880,25);
    MAT.AAE.prc75(:,j)=prctile(AAE_370_880,75);
    
    %%% BC wb %%%
    MAT.BCWB.median(:,j) = nanmedian(BC_wb(idx));
    MAT.BCWB.prc25(:,j) = prctile(BC_wb(idx),25);
    MAT.BCWB.prc75(:,j) = prctile(BC_wb(idx),75);
    
    MAT_countries.BCWB.(Country{j}) = BC_wb(idx);
    
    %%% BC ff %%%
    MAT.BCFF.median(:,j) = nanmedian(BC_ff(idx));
    MAT.BCFF.prc25(:,j) = prctile(BC_ff(idx),25);
    MAT.BCFF.prc75(:,j) = prctile(BC_ff(idx),75);
    
    MAT_countries.BCFF.(Country{j}) = BC_ff(idx);
end

return
%save 1_mat_files/MAT_DJF_AOPs.mat MAT MAT_countries size_AOPs_idx
save 1_mat_files/MAT_JJA_AOPs.mat MAT MAT_countries size_AOPs_idx

plot_AOPs_seasons
