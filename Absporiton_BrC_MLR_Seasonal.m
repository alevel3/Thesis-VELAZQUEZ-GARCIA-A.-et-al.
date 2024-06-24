clear all; close all; clc;

%load 3_EDGAR\2018_sectors.mat %%% EDGAR inventory

load 1_mat_files\Hysplit.mat %%%% Hysplit

lat=interp1(1:73,lat(:,1:end),1:1./6:73);
lon=interp1(1:73,lon(:,1:end),1:1./6:73);
alt=interp1(1:73,alt(:,1:end),1:1./6:73);
PBL=interp1(1:73,PBL(:,1:end),1:1./6:73);
rain=interp1(1:73,rain(:,1:end),1:1./6:73);


load 1_mat_files\INTERPLAY_BC.mat %%%% INTERPLAY BC

%%---
idx=find(nanmax(rain)>=1);
BC_GP(:,idx)=NaN;
size(find(isnan(sum(BC_GP))==0))%%14191
fitlm(EBC_corr(:,6)',sum(BC_GP))%0.165
%%---
%% %%%%%%%%%%% Seasons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cold season vs hot season %%
idx=find(month(time_BT)>=10|month(time_BT)<=3);%cold periods
%idx=find(month(time_BT)>=4&month(time_BT)<9);%warm periods

BC_GP = BC_GP(:,idx);
BC_GP_RCO = BC_GP_RCO(:,idx);% EBC_avg = EBC_avg (idx,:);
BC_GP_TRO_noRES = BC_GP_TRO_noRES(:,idx);
BC_GP_TRO_RES = BC_GP_TRO_RES(:,idx);
time_BT = time_BT(idx);
rain = rain(:,idx);
%% %%%%%%%%%% Aged vs Fresh aerosols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC_RESIDENTIAL = BC_GP_RCO;
BC_TRANSPORT = BC_GP_TRO_noRES+BC_GP_TRO_RES;

Hrs_fresh = 1:144;
Hrs_aged = 145:433;

Var_sect = {'BC_RESIDENTIAL_fresh','BC_RESIDENTIAL_aged','BC_TRANSPORT_fresh','BC_TRANSPORT_aged'};
Var_BC = {'BC_RESIDENTIAL','BC_TRANSPORT'};

MAT.BC_RESIDENTIAL=BC_RESIDENTIAL';
MAT.BC_TRANSPORT=BC_TRANSPORT';

for j=1:max(size(Var_sect))
    MAT_MLR.(Var_sect{j}) = [];
end
% size(find(sum(BC_TRANSPORT)==0))

%%%%%%% RESIDENTIAL %%%%%
MAT_MLR.(Var_sect{1}) = sum(MAT.(Var_BC{1})(:,Hrs_fresh),2);
MAT_MLR.(Var_sect{2}) = sum(MAT.(Var_BC{1})(:,Hrs_aged),2);
%%%%%%%% TRANSPORT %%%%%%%
MAT_MLR.(Var_sect{3}) = sum(MAT.(Var_BC{2})(:,Hrs_fresh),2);
MAT_MLR.(Var_sect{4}) = sum(MAT.(Var_BC{2})(:,Hrs_aged),2);


%% %%%%%%%%%%%%%%%%%%%%% Baseline  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BrC_MLR_baseline


%% %%% Absportion and BrC calc %%%%%

%%% → BrC_contribution_abs %%%figures

Absorption_calc

[AbsBrC_370] = BrC_abs_calc(time_BT_avg,Abs_BC6,Abs_BC1,1,880,370);
[AbsBrC_470] = BrC_abs_calc(time_BT_avg,Abs_BC6,Abs_BC2,1,880,470);

% figure;
% plot([370 470 520 590 660 880 950],[nanmean(Abs_BC1) nanmean(Abs_BC2) nanmean(Abs_BC3) nanmean(Abs_BC4) nanmean(Abs_BC5) nanmean(Abs_BC6) nanmean(Abs_BC7)]);
% hold on
% plot([370 470],[nanmean(AbsBrC_370) nanmean(AbsBrC_470)]);
%% %%%%%%% BrC lifetime %%%%%%
% size(find(MAT_MLR.BC_RESIDENTIAL_aged==0))%141

Var_sect = {'BC_RESIDENTIAL_fresh_avg','BC_RESIDENTIAL_aged_avg','BC_TRANSPORT_fresh_avg','BC_TRANSPORT_aged_avg'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MLR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -- lsqlin ---%%%
[MAEs,Abs_retrived,RC_species,CI_370,CI_470] = Absorption2_4s(MAT_MLR.(Var_sect{1}),...
    MAT_MLR.(Var_sect{2}),MAT_MLR.(Var_sect{3}),MAT_MLR.(Var_sect{4}),...
    AbsBrC_370,AbsBrC_470,Abs_BC3,Abs_BC4,Abs_BC5,Abs_BC6,Abs_BC7);

%%% -- fitlm ---%%
% [MAEs,Abs_retrived,RC_species,CI_370,CI_470,CI_880,mdl_370,mdl_470,mdl_880] = Absorption_fitlm_4s(MAT_MLR.(Var_sect{1}),...
%     MAT_MLR.(Var_sect{2}),MAT_MLR.(Var_sect{3}),MAT_MLR.(Var_sect{4}),...
%    AbsBrC_370,AbsBrC_470,Abs_BC6,NO2);

%pie_chart_CI%%370nm
%pie_chart_CI_470nm

Figure_thesis_MLR

clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,10,20])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'Fontname','Garamond',...
    'FontSize',18);
box on;
%s1=subplot(2,2,[1,2]);%%%Mass
s1=subplot(3,1,1);%%%Mass
b2=pie([nanmean(MAT_MLR.(Var_sect{1})) nanmean(MAT_MLR.(Var_sect{2})) nanmean(MAT_MLR.(Var_sect{3})) nanmean(MAT_MLR.(Var_sect{4}))])
colormap([0.7020    0.5373    0.1804;...%%R-fresh
    0.3490    0.2588    0.0667;...%%R-aged
    0.8000    0.8000    0.8000;...%%T-fresh
    0.5020    0.5020    0.5020]);%%%T-aged
labels ={'BC RESIDENTIAL fresh','BC RESIDENTIAL aged','BC TRAFFIC fresh','BC TRAFFIC aged'};
title('Relative contribution to mass (Cold season)')
lgd = legend(labels);%,'Location','bestoutside','FontWeight','Demi','FontName','times','FontSize',9);
%legend1 = legend(subplot1,'show');
set(lgd,...
    'Position',[0.0851809610232283 0.0278884372044271 0.852754251179049 0.0459154122950441],...
    'NumColumns',2);
set(gca,'FontSize',9, 'FontWeight','Demi','FontName','times','LineWidth',1.0);
set(findobj(b2,'type','text'),'FontWeight','Demi','FontName','times','FontSize',9)

% s2=subplot(2,2,3);%%RC at 370nm
s2=subplot(3,1,2);%%RC at 370nm
b2=pie(RC_species.W470)
colormap([0.7020    0.5373    0.1804;...%%R-fresh
    0.3490    0.2588    0.0667;...%%R-aged
    0.8000    0.8000    0.8000;...%%T-fresh
    0.5020    0.5020    0.5020]);%%%T-aged
title('Relative contribution to Abs BrC_{370nm} (Cold season)')
set(gca,'FontSize',9, 'FontWeight','Demi','FontName','times','LineWidth',1.0);
set(findobj(b2,'type','text'),'FontWeight','Demi','FontName','times','FontSize',9)

% s3=subplot(2,2,4);
s3=subplot(3,1,3);
b2=pie(RC_species.W880)
colormap([0.7020    0.5373    0.1804;...%%R-fresh
    0.3490    0.2588    0.0667;...%%R-aged
    0.8000    0.8000    0.8000;...%%T-fresh
    0.5020    0.5020    0.5020]);%%%T-aged

title('Relative contribution to Abs_{880nm} (Cold season)')
set(gca,'FontSize',9, 'FontWeight','Demi','FontName','times','LineWidth',1.0);
set(findobj(b2,'type','text'),'FontWeight','Demi','FontName','times','FontSize',9)
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\INTERPLAY\Plots\Seasons\MLR_Abs(Warm)');
% savefig(namesauve)
% export_fig MLR_Abs(Warm) -png -r300 -transparent




Abs_BrC = AbsBrC_370;
isgood =  ~(isnan(Abs_BrC)|isnan(Abs_retrived.W470));
fitlm(Abs_BrC(find(isgood==1)),Abs_retrived.W470(find(isgood==1)))

clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,10,8.2])
hold on
plot([0 50],[0 50],'--r','LineWidth',1);hold on
plot(Abs_BrC(find(isgood==1)),Abs_retrived.W370(find(isgood==1)),'.k');
xlim([0 20])
ylim([0 20])
P = polyfit(Abs_BrC(find(isgood==1)),Abs_retrived.W370(find(isgood==1)),1);
yfit = P(1)*Abs_BrC+P(2);
box on;
hold on;
plot(Abs_BrC(find(isgood==1)),yfit(find(isgood==1)),'b')%,'-p','LineWidth',1,'MarkerSize',3,'MarkerIndices',1:500:length(ExtG_Yao(find(isgood==1))),'Color',[0 1 1]);
legend({'1:1','data','linear fit'},'Location','northwest')
ylabel('Absorption BrC_{370nm} reconstructed (Mm^{-1})');
xlabel('Absorption BrC_{370nm} observed (Mm^{-1})');
set(gca,'FontSize',10, 'FontWeight','Demi','FontName','times','LineWidth',1.0);
% 
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\INTERPLAY\Plots\Thesis figures\Line_fit_370(Cold)');
% savefig(namesauve);
% export_fig Line_fit_370_(Cold) -png -r300 -transparent

