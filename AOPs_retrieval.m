%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTE CODIGO ERA EL MAS ACTUALIZADO 14 MARCH 2022 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figuras para el paper %%
%% Retriving the Absoption, Extinction & Scattering by wavelength%%
clear all;clc;close all
%clearvars -except WS_corr WD_corr MTO Rayleigh_G Rayleigh_B Rayleigh_R
read_Rayleigh=0;
read_MTO = 0;
%% %%%%%%%%%%%%%%%%%%%%% Previous files %%%%%%%%%%%%%%%%%%%%%
% %% 1.1 Averaged Data Variables 2017 - 2019 -> Time stamp 29min (ACSM)
% load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\Data_avaraged_2017-2019');
% %% 1.2 Averaged Data Variables 2017 - 2019 -> Time stamp 29min (ACSM) SMPS
% load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\Data_avaraged_2017-2019_EBAS');
%%load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\Data_avaraged_2017-2019_EBAS_07Jun2021');
%%load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\Data_avaraged_2017-2019_EBAS_17Sep2021');
%% 1.3 Averaged Data Variables 2017 - 2019 -> Time stamp 29min (ACSM) Latest EBAS data
%%% MLR paper subm → load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\LOA_EBASdata_avg_2017_2019_insitu17NOV2021');
%MLR_paper = load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\LOA_EBASdata_avg_2017_2019_insitu17NOV2021');
load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\LOA_EBASdata_avg_2017_2019_insitu_14MARCH2022');
%% Seasonal matrix averagaed 2017-2019
Season = {'DJF' 'MAM' 'JJA' 'SON'};
%%Seasons
idx_MAM=find(month(Time_ACSM)>=3&month(Time_ACSM)<6);%%spring
idx_DJF = find(month(Time_ACSM)>=12|month(Time_ACSM)<3);
idx_JJA=find(month(Time_ACSM)>=6&month(Time_ACSM)<9);%summer
idx_SON=find(month(Time_ACSM)>=9&month(Time_ACSM)<12);%fall
%% 1.2 Meteorological variables
if read_MTO==1
    Annee_deb = 2017;
    Mois_deb = 07;
    Jour_deb = 29;
    Annee_fin = 2019;
    Mois_fin = 12;
    Jour_fin = 31;
    Instrument = 'MTO_';
    %%Time
    [MTO] = Concatenate_matrix(Annee_deb,Mois_deb,Jour_deb,Annee_fin,Mois_fin,Jour_fin,Instrument);
    %idx_PM1 = min(find(MTO.MTO.Time_MTO>datenum(2017,7,28,23,59,0)));%(idx_PM1:end)%-->to the end point of the matrix
    %Time_MTO = MTO.Time_MTO(idx_PM1:end);
    %%1.Convertion of Orientation
    MTO.MTO.Orientation(find(MTO.MTO.Orientation == -1))= NaN;
    WD = MTO.MTO.Orientation*22.5;
    %%2. %% Converting Data to vector
    Y=MTO.MTO.W_speed.*sin(WD*pi./180);
    X=MTO.MTO.W_speed.*cos(WD*pi./180);
    %% % Here you sync X and Y coordinates with wind time stamp being time_Anem
    %%%[time_corr,X_corr,Y_corr,p1_corr,p2_corr,p3_corr,p4_corr,p5_corr,p6_corr]=time_sync(X,time_Anem,Y,time_Anem,p1,time_PMF_Data,p2,time_PMF_Data,p3,time_PMF_Data,p4,time_PMF_Data,p5,time_PMF_Data,p6,time_PMF_Data);
    %% X and Y coordinates with wind time stamp being time_ACSM
    for i=1:max(size(Time_ACSM))
        %X and Y sync
        idx=find(MTO.MTO.Time>(Time_ACSM(i)-29./(24.*60))& MTO.MTO.Time<=(Time_ACSM(i)));
        
        if size(idx,1)>0 %threshold for data coverage before averaging
            X_corr(i,1)=nanmean(X(idx));
            Y_corr(i,1)=nanmean(Y(idx));
        else
            X_corr(i,1)=NaN;
            Y_corr(i,1)=NaN;
        end
    end
    WS_corr=(X_corr.^2+Y_corr.^2).^0.5;
    WD_corr=atan2(Y_corr,X_corr).*180./pi;
    WD_corr(WD_corr<0)=WD_corr(WD_corr<0)+360;
    
    save('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\MAtrix_season_2017-2019_AVG\MTO_TmpACSM','WD_corr','-v7.3');
else
    load ('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\MAtrix_season_2017-2019_AVG\MTO_TmpACSM');
end
    %% %%%%%%%%%%%%%%%%%%%%% 2. Filters %%%%%%%%%%%%%%%%%%%%%%%%
%% 2.1 AE33 %%%
lambda = [370	470	525	590	660	880	940];
SG = [18.47 ;14.54 ;13.14 ; 11.58 ;10.35;7.77 ; 7.19];
%lim = 600;%% <= peaks in AE33
lim=0;
% %%% Removing %%
% % Abs_BC1(find(Abs_BC1>lim))=NaN;
% % Abs_BC2(find(Abs_BC2>lim))=NaN;
% % Abs_BC3(find(Abs_BC3>lim))=NaN;
% % Abs_BC4(find(Abs_BC4>lim))=NaN;
% % Abs_BC5(find(Abs_BC5>lim))=NaN;
% % Abs_BC6(find(Abs_BC6>lim))=NaN;
% % Abs_BC7(find(Abs_BC7>lim))=NaN;
% % size(find(isnan(Abs_BC3)==0))%%94931
idx=find(Abs_BC1==lim|Abs_BC2==lim|Abs_BC3==lim|Abs_BC4==lim|Abs_BC5==lim|Abs_BC6==lim|Abs_BC7==lim);
%%% Removing %%
Abs_BC1(idx)=NaN;
Abs_BC2(idx)=NaN;
Abs_BC3(idx)=NaN;
Abs_BC4(idx)=NaN;
Abs_BC5(idx)=NaN;
Abs_BC6(idx)=NaN;
Abs_BC7(idx)=NaN;
Abs_BC6(idx)=NaN;
BC6(idx)=NaN;
% % size(find(isnan(Abs_BC3)==0));%%24925
%%Zeros in AE33
% idx=find(Time_ACSM>datenum(2017,08,21,09,57,59)& Time_ACSM<=datenum(2017,08,21,16,21,28));%% ALL PERIOD
% Abs_BC1(idx)=NaN;
% Abs_BC2(idx)=NaN;
% Abs_BC3(idx)=NaN;
% Abs_BC4(idx)=NaN;
% Abs_BC5(idx)=NaN;
% Abs_BC6(idx)=NaN;
% Abs_BC7(idx)=NaN;
% Abs_BC6(idx)=NaN;
% BC6(idx)=NaN;
%% 2.2 Neph %%%
%%%%%%%%%%%%%%%%% NOTA(05 Oct 2021) %%%%%%%%%%%%%%%%%%%%%%%
%%% The chanels R and B and inverse !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% corrections of EBAS database hasn't been corrected here yet!! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lim = 600; %%<= 1peak of 1000
% %%%%%%%%% Removing %%%%%%%%%%%%
% % Scat_R(find(Scat_R>lim))=NaN;
% % Scat_G(find(Scat_G>lim))=NaN;
% % Scat_B(find(Scat_B>lim))=NaN;
% % size(find(isnan(Scat_G)==0));%%17827
%idx=find(Scat_R>lim|Scat_G>lim|Scat_B>lim);
% %%%%%%%%% Removing %%%%%%%%%%%%
% Scat_R(idx)=NaN;
% Scat_G(idx)=NaN;
% Scat_B(idx)=NaN;
% size(find(isnan(Scat_G)==0));%%17827
%% 2.3 ACSM %%%
%%limite de detection :
lim_NH4 = 0.284;% ?g/m3%   
lim_Org = 0.148;% ?g/m3%               
lim_SO4 = 0.024;% ?g/m3%               
lim_NO3 = 0.012;% ?g/m3%              
lim_Cl = 0.011;% ?g/m%                 
%%Molarities
M_molaire_SO4 = 96.06; %g/mol
M_molaire_NH4 =18.03;%g/mol
M_molaire_NO3 =62.00;%g/mol
M_molaire_Cl =35.45; %g/mol
%%%%%%%% Removing %%%%%%%%%%%
%%Doubling LD
% SO4(find(SO4<2.*lim_SO4))=NaN;
% NH4(find(NH4<2.*lim_NH4))=NaN;
% NO3(find(NO3<2.*lim_NO3))=NaN;
% Org(find(Org<2.*lim_Org))=NaN;
% Chl(find(Chl<2.*lim_Cl))=NaN;
%%Not doubling
SO4(find(SO4<lim_SO4))=NaN;
NH4(find(NH4<lim_NH4))=NaN;
NO3(find(NO3<lim_NO3))=NaN;
Org(find(Org<lim_Org))=NaN;
Chl(find(Chl<lim_Cl))=NaN;
% size(find(isnan(Org)==0));
%%%%%%%% Removing with the "magical filter" %%%%%%%%%%%
%idx=find(SO4<2.*lim_SO4|NH4<2.*lim_NH4|NO3<2.*lim_NO3|Org<2.*lim_Org|Chl<2.*lim_Cl);%...
% SO4(idx)=NaN;
% NH4(idx)=NaN;
% NO3(idx)=NaN;
% Org(idx)=NaN;
% Chl(idx)=NaN;
% size(find(isnan(Org)==0));
%% %%%%%%%%%%%%%%%%%%%%%%%%%% 3. Calculus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.1 PM1 measured & PM1 species fraction
PM1 = BC6+Org+NO3+SO4+NH4+Chl;
size(find(isnan(PM1)==0));
FIno = (NO3+SO4+NH4+Chl)./PM1;
FOrg = Org./PM1;
%% 3.2 PM1 retrieved
NH42SO4 = 2*SO4.*M_molaire_NH4./M_molaire_SO4 +SO4;
NH4NO3 =  NO3.*M_molaire_NH4./M_molaire_NO3 +NO3;
NH4Cl =  Chl.*M_molaire_NH4./M_molaire_Cl + Chl;
PM1_tot_calc = NH42SO4+NH4NO3+NH4Cl+Org+BC6;
%% 3.3 AAE 370nm-660nm (Chylek et al. 2019)
% A=log10(abs(Abs_BC1)./Abs_BC5);
% B=log10(370/660);
% AAE_370_660=-(A/B);
% figure;
% histogram(abs(AAE_370_660));
% xlim([0.5 2.5]);
% xlabel('AAE_{370-660nm}');
% set(gca,'Fontname','Garamond','Fontsize',14,'linewidth',2);
%%3.3.1 AEE 370nm-520nm
A=log10(abs(Abs_BC1)./Abs_BC3);
B=log10(370/520);
AAE_370_520=-(A/B);
%% 3.4 SAE 450-550 (Cappa et al. 2016)
%lambda2 = 525;%Scat_G
lambda3 = 635;%Scat_R
lambda2 = 525;%%SCat_G
lambda_xx = 550;
%%Scat 550
xx=log10(Scat_G./Scat_R);
yy=log10(lambda2/lambda3);
ASE1=-(xx/yy);
[Scat_550] = change_wavelength(Scat_R,ASE1,lambda3,lambda_xx);
%%SAE 450_550
a=log10(Scat_B./Scat_550);
b=log10(450/lambda_xx);
SAE_450_550=-(a/b);
%% 3.5 Absoprtion at 450nm
%%1.AEE 370-470 to get Abs_450nm
lambda1 = 370;% La longitud de onda del primer punto (Prop. Opt. Medida) 470;%AE16
lambda2 = 470;% La longitud de onda del punto final (Prop. Opt. Medida)880;%AE33
lambda_x = 450;%% La onda en la que quieres tener la propiedad optica (ej.lambda Neph)
x=log10(Abs_BC1./Abs_BC2);
y=log10(lambda1/lambda2);
AAE=-(x/y);
[Abs_450] = change_wavelength(Abs_BC2,AAE,lambda2,lambda_x);
%% 3.6 Absportion at 635
%%1.AEE 590-660 to get Abs_635nm
lambda4 = 590;% La longitud de onda del primer punto (Prop. Opt. Medida) 470;%AE16
lambda5 = 660;% La longitud de onda del punto final (Prop. Opt. Medida)880;%AE33
lambda_x = 635;%% La onda en la que quieres tener la propiedad optica (ej.lambda Neph)
x=log10(Abs_BC4./Abs_BC5);
y=log10(lambda4/lambda5);
AAE=-(x/y);
[Abs_635] = change_wavelength(Abs_BC4,AAE,lambda5,lambda_x);
%% 3.7 SSA
%% 450 nm => Blue ultraviolet range
EXT_450 = Scat_B+abs(Abs_450);
SSA_B = Scat_B./EXT_450;
%% 520-525nm => Green visible range
EXT_525 = Scat_G+Abs_BC3;%%520nm
SSA_G = Scat_G./EXT_525;
%% 635nm 
EXT_635 = Scat_R+abs(Abs_635);%%Abs_BC4(590) & Abs_BC5(660)
SSA_R = Scat_R./EXT_635;
%% 3.8 SSA vs BrC (Suzane function)
[Contrib_BrC] = BrC_calculations(Abs_BC1,Abs_BC2,Abs_BC3,Abs_BC4,Abs_BC5,Abs_BC6,Abs_BC7);
Contrib_BrC.AAE880_940(Contrib_BrC.AAE880_940<0)=NaN;
% figure;
% plot(SSA_G,Contrib_BrC.AAE880_940./PM1,'.k');
% % hold on
% % plot(SSA_G,abs(Contrib_BrC.AAE880_940)./PM1,'.r');%%...???
% ylabel('BrC_{370nm}');
% xlabel('SSA_{PM1} [520-525nm]');
% grid on
% set(gca,'Fontname','Garamond','Fontsize',15,'linewidth',2);
%% 3.9 Mean_DPg calculus 
diam=Dia./1e7;%%cm
for i=1:max(size(diam))
    surf_part(i,:)=SD(:,i).*diam(i).^2.*pi; %dS in cm2/cm3
    vol_part(i,:)=SD(:,i).*diam(i).^3.*pi./6; %dV in cm3/cm3
end

SMPS_vol=zeros(1,size(SD,1));
SMPS_surf=zeros(1,size(SD,1));

% diam_T = ones(size(SMPS_EBAS.SD));
% Try_1=diam_T.*diam;
% Try_2 = diam';

for j=1:size(SD,2)-1
    corr_bin(j)=log10(diam(j+1)/diam(j));
end

for i=1:max(size(Time_ACSM))
    %Calculating mean diameter and volume
    mean_DPg(i)=SD(i,:)*diam'./(sum(SD(i,:))).*1e7;%%%*1e7=>Back to nm3
    mean_vol_DPg(i)=diam*vol_part(:,i)./(sum(vol_part(:,i))).*1e7;

    SMPS_surf(i)=corr_bin*surf_part(1:end-1,i);
    SMPS_vol(i)=corr_bin*vol_part(1:end-1,i);
%     for j=1:size(SMPS_EBAS.SD,2)-1
%         SMPS_vol(i)=SMPS_vol(i)+vol_part(j,i).*log10(diam(j+1)/diam(j)); %Should be in grams now %Total vol
%         SMPS_surf(i)=SMPS_surf(i)+surf_part(j,i).*log10(diam(j+1)/diam(j));%%Tot Surf
%     end
end

diam=diam.*1e7; %Back to nm3

dens=1.7; %Density is 1.7 g/cm3 here
PM_SMPS_2=SMPS_vol.*1e6.*dens; %conversion to micrograms/cm3
PM_SMPS_2=PM_SMPS_2.*1e6; %conversion to micrograms/m3

[DPg_range_final,DPg_freq,DPg_data_mean,DPg_data_std,DPg_ind_test1] = frequency([40:0.2:120],mean_DPg');

Eff_Radius_Ale = (SMPS_vol')./(SMPS_surf');
%% %%%%%%%%%%%%%%%%%% 4. Retrieval of the Optical properties %%%%%%%%%%%%%
%% 4.1Rayleigh by wavelength% CO2=407ppm %change every year
if read_Rayleigh==1
Rayleigh_G = Rayleigh(Press_abs,Temperature,525,407);
Rayleigh_B = Rayleigh(Press_abs,Temperature,450,407);
Rayleigh_R = Rayleigh(Press_abs,Temperature,635,407);
save Joel/Rayleigh.mat Rayleigh_G Rayleigh_B Rayleigh_R
else
load ('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Joel\Rayleigh');%.mat 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%
%plot_season_diurnal(Time_ACSM,Abs_BC3)
%plot_season_diurnal(Time_ACSM,Scat_G)
%plot_abstract_figure
%% %%%%%%%%%%%%%%%%%% 4.2 Scattering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From 450nm to 635nm
% [mdl_SB,MSE_B,CI_SB,mdl_SG,MSE_G,CI_SG,mdl_SR,MSE_R,CI_SR,Scat_CB,Scat_CG,Scat_CR] = Scattering_4s(NH42SO4,NH4NO3,Org,BC6,Scat_B,Scat_G,Scat_R,Rayleigh_B,Rayleigh_G,Rayleigh_R);
%% NO BC
%[mdl_SB,MSE_B,CI_SB,mdl_SG,MSE_G,CI_SG,mdl_SR,MSE_R,CI_SR,Scat_CB,Scat_CG,Scat_CR] = Scattering_3s(NH42SO4,NH4NO3,Org,Scat_B,Scat_G,Scat_R,Rayleigh_B,Rayleigh_G,Rayleigh_R);
[MSE_B,CI_SB,MSE_G,CI_SG,MSE_R,CI_SR,Scat_CB,Scat_CG,Scat_CR] = Scattering_3s(NH42SO4,NH4NO3,Org,Scat_B,Scat_G,Scat_R,Rayleigh_B,Rayleigh_G,Rayleigh_R);
Scat_CB(find(isinf(Scat_CB)==1))= NaN;
Scat_CG(find(isinf(Scat_CG)==1))= NaN;
Scat_CR(find(isinf(Scat_CR)==1))= NaN;
Scat_CB(find(Scat_CB==0))= NaN;
Scat_CG(find(Scat_CG==0))= NaN;
Scat_CR(find(Scat_CR==0))= NaN;
ScatB_MLR = Scat_B./Scat_CB;
ScatG_MLR = Scat_G./Scat_CG;
ScatR_MLR = Scat_R./Scat_CR;
% ScatB_MLR = Scat_CB./Scat_B;
% ScatG_MLR = Scat_CG./Scat_B;
% ScatR_MLR = Scat_CR./Scat_B;
%% %%%%%%%%%%%%%%%%%% 4.3 Absportion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From 370nm to 880nm
[MAEs,Abs370,Abs470,Abs520,Abs590,Abs660,Abs880,Abs940] = Absorption_4s(NH42SO4,NH4NO3,Org,BC6,Abs_BC1,Abs_BC2,Abs_BC3,Abs_BC4,Abs_BC5,Abs_BC6,Abs_BC7,NO2);
Abs370(find(isinf(Abs370)==1))= NaN;
Abs470(find(isinf(Abs470)==1))= NaN;
Abs520(find(isinf(Abs520)==1))= NaN;
Abs590(find(isinf(Abs590)==1))= NaN;
Abs660(find(isinf(Abs660)==1))= NaN;
Abs880(find(isinf(Abs880)==1))= NaN;
Abs940(find(isinf(Abs940)==1))= NaN;
Abs370(find(Abs370==0))= NaN;
Abs470(find(Abs470==0))= NaN;
Abs520(find(Abs520==0))= NaN;
Abs590(find(Abs590==0))= NaN;
Abs660(find(Abs660==0))= NaN;
Abs880(find(Abs880==0))= NaN;
Abs940(find(Abs940==0))= NaN;
Abs370_MLR = Abs_BC1./Abs370;
Abs470_MLR = Abs_BC2./Abs470;
Abs520_MLR = Abs_BC3./Abs520;
Abs590_MLR = Abs_BC4./Abs590;
Abs660_MLR = Abs_BC5./Abs660;
Abs880_MLR = Abs_BC6./Abs880;
Abs940_MLR = Abs_BC7./Abs940;
%% %%%%%%%%%%%%%%%%% 4.4 Extinction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4.4.1 IMPROVE 
%% 525nm
IMPROVE_B = IMPROVE_4s(NH42SO4,NH4NO3,Org,BC6,NO2,Rayleigh_B);
IMPROVE_B(find(isinf(IMPROVE_B)==1))= NaN;
IMPROVE_B(find(IMPROVE_B==0))= NaN;
Ext_IMPROVEB = EXT_450./IMPROVE_B;
%% 525nm
IMPROVE_G = IMPROVE_4s(NH42SO4,NH4NO3,Org,BC6,NO2,Rayleigh_G);
IMPROVE_G(find(isinf(IMPROVE_G)==1))= NaN;
IMPROVE_G(find(IMPROVE_G==0))= NaN;
Ratio_IMPROVE = EXT_525./IMPROVE_G; %=> Pitford et al. 2007
%% 635nm
IMPROVE_R = IMPROVE_4s(NH42SO4,NH4NO3,Org,BC6,NO2,Rayleigh_R);
IMPROVE_R(find(isinf(IMPROVE_R)==1))= NaN;
IMPROVE_R(find(IMPROVE_R==0))= NaN;
Ext_IMPROVER = EXT_635./IMPROVE_R;
%% 4.4.2 MLR
%% EXT 450 nm =<fitlm mdl_EB
[MEE_B,CI_EB,Ext_retB] = Extinction_4s(NH42SO4,NH4NO3,Org,BC6,EXT_450,NO2,Rayleigh_B);
Ext_retB(isinf(Ext_retB)==1)= NaN;
Ext_retB(find(Ext_retB==0))= NaN;
ExtB_MLR = EXT_450./Ext_retB;
%% EXT 525nm  mdl_EG
[MEE_G,CI_EG,Ext_retG] = Extinction_4s(NH42SO4,NH4NO3,Org,BC6,EXT_525,NO2,Rayleigh_G);
Ext_retG(find(isinf(Ext_retG)==1))= NaN;
Ext_retG(find(Ext_retG==0))= NaN;
Ratio_MLR = EXT_525./Ext_retG;
%% EXT 635nm mdl_ER
[MEE_R,CI_ER,Ext_retR] = Extinction_4s(NH42SO4,NH4NO3,Org,BC6,EXT_635,NO2,Rayleigh_R);
Ext_retR(find(isinf(Ext_retR)==1))= NaN;
Ext_retR(find(Ext_retR==0))= NaN;
ExtR_MLR = EXT_635./Ext_retR;
%% OTHER AUTHORS: Yao et al 2010; Wang et al. 2015; Groblicki et al. 1981;
%%EXT 525nm or 550nm?
%% PM1
%%1.Yao et al. 2010
[MEE_Yao,ExtG_Yao] = Extinction_Yao(NH42SO4,NH4NO3,Org,BC6,EXT_525,NO2,Rayleigh_G);
ExtG_Yao(find(isinf(ExtG_Yao)==1))= NaN;
ExtG_Yao(find(ExtG_Yao==0))= NaN;
Ratio_Yao = EXT_525./ExtG_Yao;
%%2.Wang et al. 2015
[MEE_Wang,ExtG_Wang] = Extinction_Wang(NH42SO4,NH4NO3,Org,BC6,EXT_525,NO2,Rayleigh_G);
ExtG_Wang(find(isinf(ExtG_Wang)==1))= NaN;
ExtG_Wang(find(ExtG_Wang==0))= NaN;
Ratio_Wang = EXT_525./ExtG_Wang;
%%3.Groblicki et al. 1981
[MEE_Gro,ExtG_Gro] = Extinction_Groblicki(NH42SO4,NH4NO3,Org,BC6,EXT_525,NO2,Rayleigh_G);
ExtG_Gro(find(isinf(ExtG_Gro)==1))= NaN;
ExtG_Gro(find(ExtG_Gro==0))= NaN;
Ratio_Gro = EXT_525./ExtG_Gro;
%%4.Chen et al. 2008
[MEE_Cheng,ExtG_Cheng] = Extinction_Cheng(NH42SO4,NH4NO3,Org,BC6,EXT_525,NO2,Rayleigh_G);
ExtG_Cheng(find(isinf(ExtG_Cheng)==1))= NaN;
ExtG_Cheng(find(ExtG_Cheng==0))= NaN;
Ratio_Cheng = EXT_525./ExtG_Cheng;
%% PM2.5
%%Valentini et a1 2018
[MEE_Val,ExtG_Val] = Extinction_Valentini(NH42SO4,NH4NO3,Org,Abs_BC3,EXT_525,NO2,Rayleigh_G);
ExtG_Val(find(isinf(ExtG_Val)==1))= NaN;
ExtG_Val(find(ExtG_Val==0))= NaN;
Ratio_Val = EXT_525./ExtG_Val;
%%Yuan et al. 2006
[MEE_Yuan,ExtG_Yuan] = Extinction_Yuan(NH42SO4,NH4NO3,Org,BC6,EXT_525,NO2,Rayleigh_G);
ExtG_Yuan(find(isinf(ExtG_Yuan)==1))= NaN;
ExtG_Yuan(find(ExtG_Yuan==0))= NaN;
Ratio_Yuan = EXT_525./ExtG_Yuan;
%% %%%%%%%%%%%%%%%%%%%%%% 5. Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SMPS & PM1 fraction n = 2911 / Only PM1 Fraction n=6061
% isgood =  ~(isnan(Scat_B)|isnan(Scat_G)|isnan(Scat_R)|...
%     isnan(Scat_CB)|isnan(Scat_CG)|isnan(Scat_CR)|...
%     isnan(ScatB_MLR)|isnan(ScatG_MLR)|isnan(ScatR_MLR)|...
%     isnan(Abs_BC1)|isnan(Abs_BC2)|isnan(Abs_BC3)|isnan(Abs_BC4)|isnan(Abs_BC5)|isnan(Abs_BC6)|isnan(Abs_BC7)|...
%     isnan(Abs370)|isnan(Abs470)|isnan(Abs520)|isnan(Abs590)|isnan(Abs660)|isnan(Abs880)|isnan(Abs940)|...
%     isnan(Abs370_MLR)|isnan(Abs470_MLR)|isnan(Abs520_MLR)|isnan(Abs590_MLR)|isnan(Abs660_MLR)|isnan(Abs880_MLR)|isnan(Abs940_MLR)|...
%     isnan(EXT_450)|isnan(ExtB_MLR)|isnan(SSA_B)|isnan(Ext_retB)|isnan(Ext_IMPROVEB)|...
%     isnan(EXT_525)|isnan(ExtG_MLR)|isnan(SSA_G) |isnan(Ext_retG)|isnan(Ext_IMPROVEG)|...
%     isnan(EXT_635)|isnan(ExtR_MLR)|isnan(SSA_R) |isnan(Ext_retR)|isnan(Ext_IMPROVER)|...
%     isnan(FIno)|isnan(FOrg)|isnan(Eff_radius)|isnan(mean_DPg'));
%% Only optcial properties before(n = 15686) now( n = 18120)
Ext_text = EXT_525+Rayleigh_G+(NO2*0.33);
Ratio_EXT = Ext_text./Ext_retG;
isgood =  ~(isnan(Scat_B)|isnan(Scat_G)|isnan(Scat_R)|...
    isnan(Scat_CB)|isnan(Scat_CG)|isnan(Scat_CR)|...
    isnan(ScatB_MLR)|isnan(ScatG_MLR)|isnan(ScatR_MLR)|...
    isnan(Abs_BC1)|isnan(Abs_BC2)|isnan(Abs_BC3)|isnan(Abs_BC4)|isnan(Abs_BC5)|isnan(Abs_BC6)|isnan(Abs_BC7)|...
    isnan(Abs370)|isnan(Abs470)|isnan(Abs520)|isnan(Abs590)|isnan(Abs660)|isnan(Abs880)|isnan(Abs940)|...
    isnan(Abs370_MLR)|isnan(Abs470_MLR)|isnan(Abs520_MLR)|isnan(Abs590_MLR)|isnan(Abs660_MLR)|isnan(Abs880_MLR)|isnan(Abs940_MLR)|...
    isnan(EXT_450)|isnan(ExtB_MLR)|isnan(SSA_B)|isnan(Ext_retB)|isnan(Ext_IMPROVEB)|...
    isnan(EXT_525)|isnan(Ratio_MLR)|isnan(SSA_G) |isnan(Ext_retG)|isnan(Ratio_IMPROVE)|isnan(Ext_text)|...
    isnan(EXT_635)|isnan(ExtR_MLR)|isnan(SSA_R) |isnan(Ext_retR)|isnan(Ext_IMPROVER));
% % isgoodG =  ~(isnan(EXT_525)|isnan(ExtG_MLR)|isnan(SSA_G) |isnan(Ext_retG)|isnan(Ext_IMPROVEG));
% % isgoodB =  ~(isnan(EXT_450)|isnan(ExtB_MLR)|isnan(SSA_B) |isnan(Ext_retB));
% % isgoodR =  ~(isnan(EXT_635)|isnan(ExtR_MLR)|isnan(SSA_R) |isnan(Ext_retR));
%% save data
% save('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\ATOLL_2017-2019_30minAVG_paper',...
%     'BC6','Abs_BC3','Abs_450','Abs_635','Abs_BC1','Abs_BC2','Abs_BC4','Abs_BC5','Abs_BC6','Abs_BC7',...
%     'Scat_G','Scat_B','Scat_R',...
%     'Time_ACSM','Org','NH42SO4','NH4NO3',...
%     'Rayleigh_B','Rayleigh_G','Rayleigh_R',...
%     'NO2');
%% %%%%%%% ------ plots for paper -----%%%%%
%Main_figures_P1%%→code for plotting the paper ES&T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% 4.1 Scatterplot IMPROVE vs MLR %%%%%%%%%%%%%%%%%
mdl_MLR=fitlm(Ext_text(find(isgood==1)),Ext_retG(find(isgood==1)));%,'Intercept',false);%'RobustOpts','on');
mdl_IMPROVE=fitlm(EXT_525,IMPROVE_G);%,'Intercept',false);%'RobustOpts','on');
mdl=fitlm(Ext_text(find(isgood==1)),Ext_retG(find(isgood==1)),'Intercept',false)
mdl.Rsquared.Ordinary
round( mdl.Rsquared.Ordinary,2)
%mdl_MLR=fitlm(Ext_text(find(isgood==1)),Ext_retG(find(isgood==1)));%,'Intercept',false);%'RobustOpts','on');
%mdl_IMPROVE=fitlm(Ext_text,IMPROVE_G);%,'Intercept',false);%'RobustOpts','on');

% Create figure
clear fig1
fig1=figure;
%set(gcf,'units','centimeters','position',[0,0,40,25])%%[0,0,8.4,8.4] / [0.15 0.15 0.75 0.75]
% Create axes
axes1 = axes('Parent',fig1,...
    'Position',[[0.150745833333334 0.195311341260709 0.50776799756437 0.7]],...
    'LineWidth',1.5,...
    'FontWeight','Demi','FontName','times',...
    'FontSize',25);%'Fontname','Garamond',...
hold(axes1,'on');box on;
%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Dots of MLR %%%%%%%%%%%%%%%%
% MLR
plot2=plot(Ext_text(find(isgood==1)),Ext_retG(find(isgood==1)),'.','Color',[1 0.75 0.75]);%gris:0.8 0.8 0.8/%redlight:
% IMPROVE
%plot2=plot(Ext_text(find(isgood==1)),IMPROVE_G(find(isgood==1)),'.','Color',[0.729411780834198 0.831372559070587 0.95686274766922]);%0.75 0.75 1/0 0.447 0.741
hold on
%text(250,40,{'#15686'},'FontSize',25)
%%%%%%%%%%%%%%%%%%%%%%%%%% 2. Line 1:1 %%%%%%%%%%%%%%%%%%%
plot([0 300],[0 300],'--k','LineWidth',1);
xlim([0 300])
ylim([0 300])
% plot([0 250],[0 250],'--k','LineWidth',2.5);
xlim([15 300])
ylim([15 300])
%%%%%%%%%%%%%%%%%%%%%%%%%% 2.1 Line 2:2 %%%%%%%%%%%%%%%%
% hold on
% plot([0 500],[0 250],'--k','LineWidth',2.5);
% plot([0 250],[0 500],'--k','LineWidth',2.5);
%%%%%%%%%%%%%%%%%%%%%%%%% 3. Groblicki et al. 1981 %%%%%%%%%%%%%%%%%%%%%%%%% 
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Gro(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-*','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Gro(find(isgood==1))),'Color',[0.929 0.694 0.125]);%0 0.8 0.8/0.85 0.32 0.098
P = polyfit(Ext_text(find(isgood==1)),ExtG_Gro(find(isgood==1)),1);
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-*','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Gro(find(isgood==1))),'Color',[0.929 0.694 0.125]);%0 0.8 0.8/0.85 0.32 0.098
%%%%%%%%%%%%%%%%%%%%%%%%% 4. Yuan et al. 2006 %%%%%%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Yuan(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-d','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Yuan(find(isgood==1))),'Color',[0.5 0.6 0.6]);%0.2 0.8 1
% P = polyfit(Ext_text(find(isgood==1)),ExtG_Yuan(find(isgood==1)),1);
% yfit = P(1)*Ext_text+P(2);
% hold on;
% plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-d','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Yuan(find(isgood==1))),'Color',[0.5 0.6 0.6]);%0.2 0.8 1
%%%%%%%%%%%%%%%%%%%%%%%%% 5. Pitchford et al. 2007 %%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),IMPROVE_G(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-s','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(IMPROVE_G(find(isgood==1))),'Color',[0 0.447 0.741]);%0.6 0.8 1
P = polyfit(Ext_text(find(isgood==1)),IMPROVE_G(find(isgood==1)),1);
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-s','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(IMPROVE_G(find(isgood==1))),'Color',[0 0.447 0.741]);%0.6 0.8 1
%%%%%%%%%%%%%%%%%%%%%%%% 6. Chen et al. 2008 %%%%%%%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Cheng(find(isgood==1)),1);%0.4940 0.1840 0.5560
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-+','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Cheng(find(isgood==1))),'Color',[0.494 0.184 0.556]);%0.8 0.7 0.65
P = polyfit(Ext_text(find(isgood==1)),ExtG_Cheng(find(isgood==1)),1);%0.4940 0.1840 0.5560
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-+','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Cheng(isgood==1)),'Color',[0.494 0.184 0.556]);%0.8 0.7 0.65
%%%%%%%%%%%%%%%%%%%%%%%% 7. Yao et al. 2010 %%%%%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Yao(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-p','LineWidth',2,'MarkerSize',6,'MarkerIndices',1:500:length(ExtG_Yao(find(isgood==1))),'Color',[0 1 1]);%0.5 0.6 0.6
P = polyfit(Ext_text(find(isgood==1)),ExtG_Yao(find(isgood==1)),1);
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-p','LineWidth',2,'MarkerSize',6,'MarkerIndices',1:500:length(ExtG_Yao(find(isgood==1))),'Color',[0 1 1]);%0.5 0.6 0.6
%%%%%%%%%%%%%%%%%%%%% 8. Wang et al. 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Wang(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-^','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Wang(find(isgood==1))),'Color',[1 0.5 1]);%1 0.7 1
P = polyfit(Ext_text(find(isgood==1)),ExtG_Wang(find(isgood==1)),1);
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-^','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Wang(find(isgood==1))),'Color',[1 0.5 1]);%1 0.7 1
%%%%%%%%%%%%%%%%%%%% 9. Valentini et al 2018 %%%%%%%%%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Val(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-o','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Val(find(isgood==1))),'Color',[0.5 0.7 0.2]);%%0.5 0.7 0.2/0.8 0.9 0.8
P = polyfit(Ext_text(find(isgood==1)),ExtG_Val(find(isgood==1)),1);
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-o','LineWidth',2,'MarkerSize',8,'MarkerIndices',1:500:length(ExtG_Val(find(isgood==1))),'Color',[0.5 0.7 0.2]);%%0.5 0.7 0.2/0.8 0.9 0.8
%%%%%%%%%%%%%%%%%%%% 10. Velazquez-Garcia et al. 2022 %%%%%%%%%%%%
% p2=polyfit(EXT_525(find(isgood==1)),Ext_retG(find(isgood==1)),1);
% yfit = p2(1)*EXT_525+p2(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'DisplayName','linear','Tag','linear','Parent',axes1,'LineWidth',10,'Color',[1 0 0]);
p2=polyfit(Ext_text(find(isgood==1)),Ext_retG(find(isgood==1)),1);
yfit = p2(1)*Ext_text+p2(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(isgood==1),'DisplayName','linear','Tag','linear','LineWidth',3,'Color',[1 0 0]);%'Parent',axes1
%%%%%%%%%%%%%%%%%%%%%%%% 7. Yao et al. 2010 %%%%%%%%%%%%%%%%%%%%
% P = polyfit(EXT_525(find(isgood==1)),ExtG_Yao(find(isgood==1)),1);
% yfit = P(1)*EXT_525+P(2);
% plot(EXT_525(find(isgood==1)),yfit(find(isgood==1)),'-p','LineWidth',2,'MarkerSize',6,'MarkerIndices',1:500:length(ExtG_Yao(find(isgood==1))),'Color',[0 1 1]);%/1 0.9 0.55
P = polyfit(Ext_text(find(isgood==1)),ExtG_Yao(find(isgood==1)),1);
yfit = P(1)*Ext_text+P(2);
hold on;
plot(Ext_text(find(isgood==1)),yfit(find(isgood==1)),'-p','LineWidth',2,'MarkerSize',6,'MarkerIndices',1:500:length(ExtG_Yao(find(isgood==1))),'Color',[0 1 1]);%/1 0.9 0.55
%%%%%%%%%%%%%%%%%%% Labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylabel('Reconstructed extinction_{525nm} (Mm^{-1})','FontWeight','Demi','FontName','times','Fontsize',25)
xlabel('Measured extinction_{525nm} (Mm^{-1})','FontWeight','Demi','FontName','times','Fontsize',25)%%10-12
legend({'','1:1','Groblicki et al. 1981','Pitchford et al. 2007','Cheng et al. 2008',...
    'Yao et al. 2010','Wang et al. 2015','Valentini et al. 2018','This work (MLR)',''},...
    'Location','best','FontWeight','Demi','FontName','times','FontSize',6);% 'Yuan et al. 2006'/RMSE:0.83%'y=0.87x+0.08 r^{2}=0.72'

set(gca,'FontSize',30, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\ObservedvsMLRvsLiterature_525(30NOV2021)');
% export_fig ObservedvsMLRvsLiterature_525(28FEB2021) -png -r600 -transparent%%-png
% savefig(namesauve);
% %print( '-djpeg', '-r650', namesauve);
%% Ext measured vs retrieved
% legend({'','1:1','Pitchford et al. 2007'},...
%     'Location','bestoutside','FontWeight','Demi','FontName','times','FontSize',35);
% legend({'','1:1','Groblicki et al. 1981','Yuan et al. 2006','Pitchford et al. 2007','Cheng et al. 2008',...
%     'Yao et al. 2010','Wang et al. 2015','Valentini et al. 2018'},...
%     'Location','bestoutside','FontWeight','Demi','FontName','times','FontSize',35);% RMSE:0.83%'y=0.87x+0.08 r^{2}=0.72'
% legend({'','1:1','Yao et al. 2010','This work (MLR)',''},...
%     'Location','bestoutside','FontWeight','Demi','FontName','times','FontSize',35);% RMSE:0.83%'y=0.87x+0.08 r^{2}=0.72'
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\Literature(new)');
% export_fig Literature(new) -png -r600 -transparent%%-png
% savefig(namesauve);
% print( '-djpeg', '-r650', namesauve);
%% %%%%%%%%%% 5.2 Contribution of the chemical composition to the AOPs %%%%%%%
labels = {'NH_{4}(2)SO_{4}','NH_{4}NO_{3}','Org','eBC'};
%%SPECIES
MAT.NH42SO4 = NH42SO4;
MAT.NH4NO3 = NH4NO3;
MAT.Org = Org;
MAT.BC6 = BC6;
MAT.NO3 = NO3;
MAT.NH4 = NH4;
MAT.SO4 = SO4;
MAT.Chl = Chl;
MAT.PM1_calc = NH42SO4+NH4NO3+Org+BC6;
MAT.PM1_meas = Org+BC6+NH4+NO3+SO4+Chl;
MAT.ExtG = EXT_525;
MAT.ExtG_MLR = Ext_retG;
MAT.ExtG_IMPROVE = IMPROVE_G;
MAT.RH = humidity;
%%MEAN PM1_CALCULADA
Var  = {'NH42SO4' 'NH4NO3' 'Org' 'BC6'};
for j=1:max(size(Var))
    MAT.PM1_4s(j) = nanmean(MAT.(Var{j}));
    MAT.PM1_4s_Porc(j) = 100.*nanmean(MAT.(Var{j})./MAT.PM1_calc);
end
%%MEAN PM1 MEDIDA
Var  = {'Org' 'NO3' 'NH4' 'SO4' 'Chl' 'BC6'};
for j=1:max(size(Var))
    MAT.PM1_6s(j) = nanmean(MAT.(Var{j}));
    MAT.PM1_avg(j) = 100.*nanmean(MAT.(Var{j})./MAT.PM1_meas);
end
%% 5.2.1 Scattering (3 wavelenght, bar graph)
Var  = {'NH42SO4' 'NH4NO3' 'Org'};
labels = {'NH_{4}(2)SO_{4}','NH_{4}NO_{3}','Org'};
%%BC
% y = [MAT.PM1_4s(1)*MSE_B(1) MAT.PM1_4s(2)*MSE_B(2) MAT.PM1_4s(3)*MSE_B(3) MAT.PM1_4s(4)*MSE_B(4);...
%     MAT.PM1_4s(1)*MSE_G(1) MAT.PM1_4s(2)*MSE_G(2) MAT.PM1_4s(3)*MSE_G(3) MAT.PM1_4s(4)*MSE_G(4);...
%     MAT.PM1_4s(1)*MSE_R(1) MAT.PM1_4s(2)*MSE_R(2) MAT.PM1_4s(3)*MSE_R(3) MAT.PM1_4s(4)*MSE_R(4)]';
% x=[nansum(y(:,1)) nansum(y(:,2)) nansum(y(:,3))]%
% Cont_Scat = [y(1)*100./x(1,1) y(2,1)*100./x(1) y(3,1)*100./x(1) y(4,1)*100./x(1);...
%     y(1,2)*100./x(2) y(2,2)*100./x(2) y(3,2)*100./x(2) y(4,2)*100./x(2);...
%     y(1,3)*100./x(3) y(2,3)*100./x(3) y(3,3)*100./x(3) y(4,3)*100./x(3)]';
%%NO BC
y = [MAT.PM1_4s(1)*MSE_B(1) MAT.PM1_4s(2)*MSE_B(2) MAT.PM1_4s(3)*MSE_B(3);...
    MAT.PM1_4s(1)*MSE_G(1) MAT.PM1_4s(2)*MSE_G(2) MAT.PM1_4s(3)*MSE_G(3);...
    MAT.PM1_4s(1)*MSE_R(1) MAT.PM1_4s(2)*MSE_R(2) MAT.PM1_4s(3)*MSE_R(3)]';
x=[nansum(y(:,1)) nansum(y(:,2)) nansum(y(:,3))]%
Cont_Scat = [y(1)*100./x(1,1) y(2,1)*100./x(1) y(3,1)*100./x(1);...
             y(1,2)*100./x(2) y(2,2)*100./x(2) y(3,2)*100./x(2);...
             y(1,3)*100./x(3) y(2,3)*100./x(3) y(3,3)*100./x(3)];
%Wav = [450; 525; 635]';
%%%PLOT
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(gcf,'units','centimeters','position',[0,0,100,100])
% Create multiple lines using matrix input to bar
bar1 = bar(Cont_Scat,'stacked');
set(bar1(1),'DisplayName','NH_{4}(2)SO_{4}',...
    'FaceColor',[0.9 0 0]);
set(bar1(2),'DisplayName','NH_{4}NO_{3}','FaceColor',[0 0.2 0.9]);
set(bar1(3),'DisplayName','Org',...
    'FaceColor',[0.200000002980232 0.800000011920929 0]);
%set(bar1(4),'DisplayName','eBC','FaceColor',[0 0 0]);
%legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northeastoutside','FontSize',15);
%legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org'},'Location','northoutside','Orientation','horizontal','FontSize',35);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','FontSize',35
legend({'Ammonium sulfate','Ammonium nitrate','Organics','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',35);%
%title({'Contribution of the chemical composition to the σ_{scatt}'});
xticklabels({'','','450nm','','525nm','','635nm'});box on;
%ylabel('m^{2}g^{-1}');ylabel('µg^-3');
ylabel('Contribution to Scatt Coeff (%)');%xlabel('Wavelength');
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_scattering(17NOV2021)');
% export_fig PM1_scattering(17NOV2021) -png -r600 -transparent
% savefig(namesauve);
%% 5.2.2 Absorption (7 wavelenght, area graph)
Var  = {'NH42SO4' 'NH4NO3' 'Org' 'BC6'};
labels = {'NH_{4}(2)SO_{4}','NH_{4}NO_{3}','Org','eBC'};
% Y=[2 3 4;7 6 5]';
% area(1:3,Y)
% Cont_Abs 
y = [MAT.PM1_4s(1)*MAEs(1) MAT.PM1_4s(2)*MAEs(2,1) MAT.PM1_4s(3)*MAEs(3,1) MAT.PM1_4s(4)*MAEs(4,1);...
    MAT.PM1_4s(1)*MAEs(1,2) MAT.PM1_4s(2)*MAEs(2,2) MAT.PM1_4s(3)*MAEs(3,2) MAT.PM1_4s(4)*MAEs(4,2);...
    MAT.PM1_4s(1)*MAEs(1,3) MAT.PM1_4s(2)*MAEs(2,3) MAT.PM1_4s(3)*MAEs(3,3) MAT.PM1_4s(4)*MAEs(4,3);...
    MAT.PM1_4s(1)*MAEs(1,4) MAT.PM1_4s(2)*MAEs(2,4) MAT.PM1_4s(3)*MAEs(3,4) MAT.PM1_4s(4)*MAEs(4,4);...
    MAT.PM1_4s(1)*MAEs(1,5) MAT.PM1_4s(2)*MAEs(2,5) MAT.PM1_4s(3)*MAEs(3,5) MAT.PM1_4s(4)*MAEs(4,5);...
    MAT.PM1_4s(1)*MAEs(1,6) MAT.PM1_4s(2)*MAEs(2,6) MAT.PM1_4s(3)*MAEs(3,6) MAT.PM1_4s(4)*MAEs(4,6);...
    MAT.PM1_4s(1)*MAEs(1,7) MAT.PM1_4s(2)*MAEs(2,7) MAT.PM1_4s(3)*MAEs(3,7) MAT.PM1_4s(4)*MAEs(4,7)]';

% y = [MAT.PM1_4s(1)*round(MAEs(1),1) MAT.PM1_4s(2)*round(MAEs(2,1),1) MAT.PM1_4s(3)*round(MAEs(3,1),1) MAT.PM1_4s(4)*round(MAEs(4,1),1);...
%     MAT.PM1_4s(1)*round(MAEs(1,2),1) MAT.PM1_4s(2)*round(MAEs(2,2),1) MAT.PM1_4s(3)*round(MAEs(3,2),1) MAT.PM1_4s(4)*round(MAEs(4,2),1);...
%     MAT.PM1_4s(1)*round(MAEs(1,3),1) MAT.PM1_4s(2)*round(MAEs(2,3),1) MAT.PM1_4s(3)*round(MAEs(3,3),1) MAT.PM1_4s(4)*round(MAEs(4,3),1);...
%     MAT.PM1_4s(1)*round(MAEs(1,4),1) MAT.PM1_4s(2)*round(MAEs(2,4),1) MAT.PM1_4s(3)*round(MAEs(3,4),1) MAT.PM1_4s(4)*round(MAEs(4,4),1);...
%     MAT.PM1_4s(1)*round(MAEs(1,5),1) MAT.PM1_4s(2)*round(MAEs(2,5),1) MAT.PM1_4s(3)*round(MAEs(3,5),1) MAT.PM1_4s(4)*round(MAEs(4,5),1);...
%     MAT.PM1_4s(1)*round(MAEs(1,6),1) MAT.PM1_4s(2)*round(MAEs(2,6),1) MAT.PM1_4s(3)*round(MAEs(3,6),1) MAT.PM1_4s(4)*round(MAEs(4,6),1);...
%     MAT.PM1_4s(1)*round(MAEs(1,7),1) MAT.PM1_4s(2)*round(MAEs(2,7),1) MAT.PM1_4s(3)*round(MAEs(3,7),1) MAT.PM1_4s(4)*round(MAEs(4,7),1)]';

x=[nansum(y(:,1)) nansum(y(:,2)) nansum(y(:,3)) nansum(y(:,4)) nansum(y(:,5)) nansum(y(:,6)) nansum(y(:,7))];
Cont_Abs = [y(1,1)*100./x(1,1) y(2,1)*100./x(1) y(3,1)*100./x(1) y(4,1)*100./x(1);...
    y(1,2)*100./x(2) y(2,2)*100./x(2) y(3,2)*100./x(2) y(4,2)*100./x(2);...
    y(1,3)*100./x(3) y(2,3)*100./x(3) y(3,3)*100./x(3) y(4,3)*100./x(3);...
    y(1,4)*100./x(4) y(2,4)*100./x(4) y(3,4)*100./x(4) y(4,4)*100./x(4);...
    y(1,5)*100./x(5) y(2,5)*100./x(5) y(3,5)*100./x(5) y(4,5)*100./x(5);...
    y(1,6)*100./x(6) y(2,6)*100./x(6) y(3,6)*100./x(6) y(4,6)*100./x(6);...
    y(1,7)*100./x(7) y(2,7)*100./x(7) y(3,7)*100./x(7) y(4,7)*100./x(7)];
%Cont_Abs(find(Cont_Abs<=0.2))=NaN;
clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,40,25])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'Fontname','Garamond',...
    'FontSize',18);
box on;
plot(1:7,Cont_Abs,'LineWidth',5);
%shadedplot(1:7,MAT_mean22.DJF.(Var{lol}).prct25,MAT_mean22.DJF.(Var{lol}).prct75,[0.7 0.7 1],[0.7 0.7 1]);
%area(1:7,Cont_Abs);
% plot(1:7,Cont_Abs)
ylim([0 100]);
newcolors = [0.200000002980232 0.800000011920929 0;0 0 0];
colororder(newcolors)
ylabel('Contribution to Abs Coeff (%)','FontSize',45);
%xlabel('Wavelength')
xticklabels({'370nm','470nm','520nm','590nm','660nm','880nm','940nm'});box on;
legend({'Organics','eBC'},'Location','northeastoutside','FontSize',35);
%title({'Contribution of the chemical composition to the σ_{abs}'});
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_absportion(18NOV2021)');
% export_fig PM1_absportion(01MARCH2022) -png -r600 -transparent
% savefig(namesauve);
%% Extinction by wavelength
y = [MAT.PM1_4s(1)*MEE_B(1) MAT.PM1_4s(2)*MEE_B(2) MAT.PM1_4s(3)*MEE_B(3) MAT.PM1_4s(4)*MEE_B(4);...
    MAT.PM1_4s(1)*MEE_G(1) MAT.PM1_4s(2)*MEE_G(2) MAT.PM1_4s(3)*MEE_G(3) MAT.PM1_4s(4)*MEE_G(4);...
    MAT.PM1_4s(1)*MEE_R(1) MAT.PM1_4s(2)*MEE_R(2) MAT.PM1_4s(3)*MEE_R(3) MAT.PM1_4s(4)*MEE_R(4)]';
x=[nansum(y(:,1)) nansum(y(:,2)) nansum(y(:,3))];
Cont_Ext = [y(1)*100./x(1,1) y(2,1)*100./x(1) y(3,1)*100./x(1) y(4,1)*100./x(1);...
    y(1,2)*100./x(2) y(2,2)*100./x(2) y(3,2)*100./x(2) y(4,2)*100./x(2);...
    y(1,3)*100./x(3) y(2,3)*100./x(3) y(3,3)*100./x(3) y(4,3)*100./x(3)];
Ext_520_525 = [y(1,2)*100./x(2) y(2,2)*100./x(2) y(3,2)*100./x(2) y(4,2)*100./x(2)];
Ext_Total = [MAT.PM1_4s(1)*MEE_G(1) MAT.PM1_4s(2)*MEE_G(2) MAT.PM1_4s(3)*MEE_G(3) MAT.PM1_4s(4)*MEE_G(4) nanmean(Rayleigh_G) (0.33*nanmean(NO2))];
%pie(Cont_Ext(2,:))
%Wav = [450; 525; 635]';
%%%PLOT
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(gcf,'units','centimeters','position',[0,0,100,100])
% Create multiple lines using matrix input to bar
bar1 = bar(Cont_Ext,'stacked');
set(bar1(1),'DisplayName','NH_{4}(2)SO_{4}',...
    'FaceColor',[0.9 0 0]);
set(bar1(2),'DisplayName','NH_{4}NO_{3}','FaceColor',[0 0.2 0.9]);
set(bar1(3),'DisplayName','Org',...
    'FaceColor',[0.200000002980232 0.800000011920929 0]);
set(bar1(4),'DisplayName','eBC','FaceColor',[0 0 0]);
%set(bar1(4),'DisplayName','eBC','FaceColor',[0 0 0]);
%legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northeastoutside','FontSize',15);
%legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',50);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','FontSize',35
legend({'Ammonium sulfate','Ammonium nitrate','Organics','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',35);%
%title({'Contribution of the chemical composition to the σ_{scatt}'});
xticklabels({'','','450nm','','525nm','','635nm'});box on;
%ylabel('m^{2}g^{-1}');ylabel('µg^-3');
ylabel('Contribution to Ext Coeff (%)');%xlabel('Wavelength');
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\Ext_by_Wavelengths(17NOV2021)');
% export_fig Ext_by_Wavelengths(17NOV2021) -png -r600 -transparent
% savefig(namesauve);
%% %%%%%%%%%% 5.3 Relative contribution of each species acording literature %%%
MEE_Pit = [2.2 2.4 2.8 10];%%=>Pitchofrt et al 2007 
Cont_Pitchf = [MAT.PM1_4s(1)*MEE_Pit(1),MAT.PM1_4s(2)*MEE_Pit(2),MAT.PM1_4s(3)*MEE_Pit(3),MAT.PM1_4s(4)*MEE_Pit(4)]; 
Cont_Yao = [MAT.PM1_4s(1)*MEE_Yao(1),MAT.PM1_4s(2)*MEE_Yao(2),MAT.PM1_4s(3)*MEE_Yao(3),MAT.PM1_4s(4)*MEE_Yao(4)];
Cont_Wang = [MAT.PM1_4s(1)*MEE_Wang(1),MAT.PM1_4s(2)*MEE_Wang(2),MAT.PM1_4s(3)*MEE_Wang(3),MAT.PM1_4s(4)*MEE_Wang(4)];
Cont_Gro = [MAT.PM1_4s(1)*MEE_Gro(1),MAT.PM1_4s(2)*MEE_Gro(2),MAT.PM1_4s(3)*MEE_Gro(3),MAT.PM1_4s(4)*MEE_Gro(4)];
Cont_Cheng = [MAT.PM1_4s(1)*MEE_Cheng(1),MAT.PM1_4s(2)*MEE_Cheng(2),MAT.PM1_4s(3)*MEE_Cheng(3),MAT.PM1_4s(4)*MEE_Cheng(4)];
Cont_Valentini = [MAT.PM1_4s(1)*MEE_Val(1),MAT.PM1_4s(2)*MEE_Val(2),MAT.PM1_4s(3)*MEE_Val(3),MAT.PM1_4s(4)*13.14];%7.7
Cont_Yuan = [MAT.PM1_4s(1)*MEE_Yuan(1),MAT.PM1_4s(2)*MEE_Yuan(2),MAT.PM1_4s(3)*MEE_Yuan(3),MAT.PM1_4s(4)*MEE_Yuan(4)];
Cont_MLR = [MAT.PM1_4s(1)*MEE_G(1),MAT.PM1_4s(2)*MEE_G(2),MAT.PM1_4s(3)*MEE_G(3),MAT.PM1_4s(4)*MEE_G(4)];
Cont_MLR = 100.*(Cont_MLR./nansum(37.9270))
%%%%%%%%%%%% Cont_Pitchf*100/dataSum %%%%%%%%%%%%%%%%%
dataSum = [nansum(Cont_Gro) nansum(Cont_Yuan) nansum(Cont_Pitchf) nansum(Cont_Cheng) nansum(Cont_Yao) nansum(Cont_Wang) nansum(Cont_Valentini)];
PCont_Pitchf = Cont_Pitchf*100/dataSum(:,1);
PCont_Yao = Cont_Yao*100/dataSum(:,2);
PCont_Wang = Cont_Wang*100/dataSum(:,3);
PCont_Gro = Cont_Gro*100/dataSum(:,4);
PCont_Cheng = Cont_Cheng*100/dataSum(:,5);
PCont_Valentini = Cont_Valentini*100/dataSum(:,6);
PCont_Yuan = Cont_Yuan*100/dataSum(:,7);
%%%%%% Total Concentrations in (%) %%%%%%%%%%%%%%%%%
Cont_NH42SO4 = [PCont_Gro(:,1),PCont_Yuan(:,1),PCont_Pitchf(:,1),PCont_Cheng(:,1),PCont_Yao(:,1),PCont_Wang(:,1),PCont_Valentini(:,1)]';
Cont_NH4NO3 = [PCont_Gro(:,2),PCont_Yuan(:,2),PCont_Pitchf(:,2),PCont_Cheng(:,2),PCont_Yao(:,2),PCont_Wang(:,2),PCont_Valentini(:,2)]';
Cont_Org = [PCont_Gro(:,3),PCont_Yuan(:,3),PCont_Pitchf(:,3),PCont_Cheng(:,3),PCont_Yao(:,3),PCont_Wang(:,3),PCont_Valentini(:,3)]';
Cont_eBC = [PCont_Gro(:,4),PCont_Yuan(:,4),PCont_Pitchf(:,4),PCont_Cheng(:,4),PCont_Yao(:,4),PCont_Wang(:,4),PCont_Valentini(:,4)]';
Cont_PM1 = [Cont_NH42SO4 Cont_NH4NO3 Cont_Org Cont_eBC];
%%% BOX PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = Cont_PM1;
x = 1:4;
colors = [0.9 0 0;0 0.2 0.9;0.200000002980232 0.800000011920929 0;0 0 0];%%4 filas* 3 columnas
figure();
ax = axes();
hold(ax);
for i=1:4
    boxchart(x(i)*ones(size(data(:,i))), data(:,i), 'BoxFaceColor', colors(i,:),'orientation','horizontal')
end
box on
% boxplot(Cont_PM1,'Labels',{'Ammonium sulfate','Ammonium nitrate','Organics','eBC'},'orientation','horizontal')
% legend({'Ammonium sulfate','Ammonium nitrate','Organics','eBC','This work','','',''},'Location','bestoutside','Orientation','vertical','FontSize',35);%
xlabel('Contribution to Ext Coeff_{525nm} (%)')
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
xlim([0 100]);
set(gcf,'units','centimeters','position',[0,0,40,25]);
yticklabels({'','','',''})
%%%%%%%% MLR markers %%%%%%%%%%
hold on
plot(15,1,'p','MarkerEdgeColor','red','MarkerFaceColor',[1 0 1],'MarkerSize',20);
plot(34,2,'p','MarkerEdgeColor','red','MarkerFaceColor',[1 0 1],'MarkerSize',20);
plot(25,3,'p','MarkerEdgeColor','red','MarkerFaceColor',[1 0 1],'MarkerSize',20);%0.8500 0.3250 0.0980
plot(26,4,'p','MarkerEdgeColor','red','MarkerFaceColor',[1 0 1],'MarkerSize',20);%1 .6 .6
legend({'Ammonium sulfate','Ammonium nitrate','Organics','eBC','This work','','',''},'Location','bestoutside','Orientation','vertical','FontSize',35);%
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\Uncertanities_speciesMLR(17NOV2021)');
% export_fig Uncertanities_speciesMLR(17NOV2021) -png -r600 -transparent
% savefig(namesauve);
%% Contributions Literature Pie charts (%)
%close all
Pie_charts_authors

%% 5.3.3 Extinction (montly, area graph) MONTHS!!!!!!
%% A. Matix in months
MAT.Time = Time_ACSM;
MAT_Month = sort_month(MAT.Time,MAT);
monthofyear = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'} ;
% MEE_DJF = [10.48 7.61 5.30 18.62; 11.52 5.30 3.89 15.38; 12.12 3.96 2.42 15.41]';
% MEE_MAM = [10.68 0.28 3.83 18.44; 8.60 0.27 2.37 15.51; 6.34 0.17 1.25 17.51]';
% MEE_JJA = [8.29 5.53 2.11 21.97; 6.44 3.99 1.42 19.04; 5.12 2.96 0.86 17.25]';
% MEE_SON =[8.31 7.44 2.50 20.28; 6.96 5.35 1.53 17.59; 6.68 3.56 1.00 16.73]';
%%Seasons
% MEE_DJF = [11.52 5.30 3.89 15.38]';
% MEE_MAM = [8.60 0.27 2.37 15.51]';
% MEE_JJA = [6.44 3.99 1.42 19.04]';
% MEE_SON =[6.96 5.35 1.53 17.59]';
%%All
MEE_DJF = [MEE_G(1) MEE_G(2) MEE_G(3) MEE_G(4)]';
MEE_MAM = [MEE_G(1) MEE_G(2) MEE_G(3) MEE_G(4)]';
MEE_JJA = [MEE_G(1) MEE_G(2) MEE_G(3) MEE_G(4)]';
MEE_SON =[MEE_G(1) MEE_G(2) MEE_G(3) MEE_G(4)]';
%%Average PM1 calc
Var  = {'NH42SO4' 'NH4NO3' 'Org' 'BC6'};
for i=1:max(size(monthofyear))
    for j=1:max(size(Var))
    MAT_Month.(monthofyear{i}).PM1_4s(j) = nanmean(MAT_Month.(monthofyear{i}).(Var{j}));
    MAT_Month.(monthofyear{i}).PM1_AVG = nanmean(MAT_Month.(monthofyear{i}).PM1_calc);
    end
end
%%Average PM1 meas
Var  = {'Org','NO3','NH4','SO4','Chl','BC6'};
for i=1:max(size(monthofyear))
    for j=1:max(size(Var))
    MAT_Month.(monthofyear{i}).PM1_AVG = nanmean(MAT_Month.(monthofyear{i}).PM1_meas);
    MAT_Month.(monthofyear{i}).PM1_ACSM(j) = nanmean(MAT_Month.(monthofyear{i}).(Var{j}));
    MAT_Month.(monthofyear{i}).PM1_ACSM_std(j) = nanstd(MAT_Month.(monthofyear{i}).(Var{j}));
    
    %MAT_Month.(monthofyear{i}).PM1_ACSM(j) = nanmedian(MAT_Month.(monthofyear{i}).(Var{j}));
    
    end
end
% %%Average Optical properties
% Var = {'ExtG','ExtG_MLR','ExtG_IMPROVE'};
% for i=1:max(size(monthofyear))
%     for j=1:max(size(Var))
%     %MAT_Month.(monthofyear{i}).PM1_AVG = nanmean(MAT_Month.(monthofyear{i}).PM1_meas);
%     MAT_Month.(monthofyear{i}).Extinction(j) = nanmean(MAT_Month.(monthofyear{i}).(Var{j}));
%     end
% end
%% C. Calculating contribution in %
%%wINTER
monthofyear = {'Jan' 'Feb' 'Dec'};
for i=1:max(size(monthofyear))
    %for j=1:max(size(Var))
    MAT_Month.(monthofyear{i}).y(1) = MAT_Month.(monthofyear{i}).PM1_4s(1)*MEE_DJF(1);
    MAT_Month.(monthofyear{i}).y(2) = MAT_Month.(monthofyear{i}).PM1_4s(2)*MEE_DJF(2);
    MAT_Month.(monthofyear{i}).y(3) = MAT_Month.(monthofyear{i}).PM1_4s(3)*MEE_DJF(3);
    MAT_Month.(monthofyear{i}).y(4) = MAT_Month.(monthofyear{i}).PM1_4s(4)*MEE_DJF(4);
    MAT_Month.(monthofyear{i}).x = [nansum(MAT_Month.(monthofyear{i}).y)];
    MAT_Month.(monthofyear{i}).Cont_Ext(1) = MAT_Month.(monthofyear{i}).y(1)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(2) = MAT_Month.(monthofyear{i}).y(2)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(3) = MAT_Month.(monthofyear{i}).y(3)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(4) = MAT_Month.(monthofyear{i}).y(4)*100/MAT_Month.(monthofyear{i}).x;
    %end
end
%%SPRING
monthofyear = {'Mar' 'Apr' 'May'};
for i=1:max(size(monthofyear))
    %for j=1:max(size(Var))
    MAT_Month.(monthofyear{i}).y(1) = MAT_Month.(monthofyear{i}).PM1_4s(1)*MEE_MAM(1);
    MAT_Month.(monthofyear{i}).y(2) = MAT_Month.(monthofyear{i}).PM1_4s(2)*MEE_MAM(2);
    MAT_Month.(monthofyear{i}).y(3) = MAT_Month.(monthofyear{i}).PM1_4s(3)*MEE_MAM(3);
    MAT_Month.(monthofyear{i}).y(4) = MAT_Month.(monthofyear{i}).PM1_4s(4)*MEE_MAM(4);
    MAT_Month.(monthofyear{i}).x = [nansum(MAT_Month.(monthofyear{i}).y)];
    MAT_Month.(monthofyear{i}).Cont_Ext(1) = MAT_Month.(monthofyear{i}).y(1)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(2) = MAT_Month.(monthofyear{i}).y(2)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(3) = MAT_Month.(monthofyear{i}).y(3)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(4) = MAT_Month.(monthofyear{i}).y(4)*100/MAT_Month.(monthofyear{i}).x;
    %end
end
%%SUMMER
monthofyear = {'Jun' 'Jul' 'Aug'};
for i=1:max(size(monthofyear))
    %for j=1:max(size(Var))
    MAT_Month.(monthofyear{i}).y(1) = MAT_Month.(monthofyear{i}).PM1_4s(1)*MEE_JJA(1);
    MAT_Month.(monthofyear{i}).y(2) = MAT_Month.(monthofyear{i}).PM1_4s(2)*MEE_JJA(2);
    MAT_Month.(monthofyear{i}).y(3) = MAT_Month.(monthofyear{i}).PM1_4s(3)*MEE_JJA(3);
    MAT_Month.(monthofyear{i}).y(4) = MAT_Month.(monthofyear{i}).PM1_4s(4)*MEE_JJA(4);
    MAT_Month.(monthofyear{i}).x = [nansum(MAT_Month.(monthofyear{i}).y)];
    MAT_Month.(monthofyear{i}).Cont_Ext(1) = MAT_Month.(monthofyear{i}).y(1)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(2) = MAT_Month.(monthofyear{i}).y(2)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(3) = MAT_Month.(monthofyear{i}).y(3)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(4) = MAT_Month.(monthofyear{i}).y(4)*100/MAT_Month.(monthofyear{i}).x;

    %end
end
%%FALL
monthofyear = {'Sep' 'Oct' 'Nov'};
for i=1:max(size(monthofyear))
    %for j=1:max(size(Var))
    MAT_Month.(monthofyear{i}).y(1) = MAT_Month.(monthofyear{i}).PM1_4s(1)*MEE_SON(1);
    MAT_Month.(monthofyear{i}).y(2) = MAT_Month.(monthofyear{i}).PM1_4s(2)*MEE_SON(2);
    MAT_Month.(monthofyear{i}).y(3) = MAT_Month.(monthofyear{i}).PM1_4s(3)*MEE_SON(3);
    MAT_Month.(monthofyear{i}).y(4) = MAT_Month.(monthofyear{i}).PM1_4s(4)*MEE_SON(4);
    MAT_Month.(monthofyear{i}).x = [nansum(MAT_Month.(monthofyear{i}).y)];
    MAT_Month.(monthofyear{i}).Cont_Ext(1) = MAT_Month.(monthofyear{i}).y(1)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(2) = MAT_Month.(monthofyear{i}).y(2)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(3) = MAT_Month.(monthofyear{i}).y(3)*100/MAT_Month.(monthofyear{i}).x;
    MAT_Month.(monthofyear{i}).Cont_Ext(4) = MAT_Month.(monthofyear{i}).y(4)*100/MAT_Month.(monthofyear{i}).x;
    %end
end
Cont_Ext = [MAT_Month.Jan.Cont_Ext; MAT_Month.Feb.Cont_Ext; MAT_Month.Mar.Cont_Ext; MAT_Month.Apr.Cont_Ext;...
    MAT_Month.May.Cont_Ext; MAT_Month.Jun.Cont_Ext; MAT_Month.Jul.Cont_Ext; MAT_Month.Aug.Cont_Ext;...
    MAT_Month.Sep.Cont_Ext; MAT_Month.Oct.Cont_Ext; MAT_Month.Nov.Cont_Ext; MAT_Month.Dec.Cont_Ext];
%% D. Plotting in '%'
clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,150,100])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'FontSize',18);
box on;
area(1:12,Cont_Ext);
ylim([0 100]);
xlim([1 12]);
newcolors = [0.9 0 0;0 0.2 0.9;0.200000002980232 0.800000011920929 0;0 0 0];
colororder(newcolors)
ylabel('Contribution Ext Coeff_{525nm} (%)');
xlabel('Months')
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
%xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});box on;
legend({'Ammonium sulfate','Ammonium nitrate','Organics','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',35);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
%legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',50);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
%title({'Contribution of the chemical composition to the σ_{ext 525-520nm}'});
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_year(17NOV2021)');
% export_fig Abstract_graph -png -r600 -transparent%%-png r600
% savefig(namesauve);
%% D. Plotting in Mm-1
Cont_Ext_M = [MAT_Month.Jan.y; MAT_Month.Feb.y; MAT_Month.Mar.y; MAT_Month.Apr.y;...
    MAT_Month.May.y; MAT_Month.Jun.y; MAT_Month.Jul.y; MAT_Month.Aug.y;...
    MAT_Month.Sep.y; MAT_Month.Oct.y; MAT_Month.Nov.y; MAT_Month.Dec.y];
clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,150,100])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'FontSize',18);
box on;
area(1:12,Cont_Ext_M);
ylim([0 80]);
xlim([1 12]);
newcolors = [0.9 0 0;0 0.2 0.9;0.200000002980232 0.800000011920929 0;0 0 0];
colororder(newcolors)
ylabel('Extinction_{525nm} (Mm^{-1})');
%xlabel('Months')
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
%xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});box on;
%legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',50);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
legend({'AS','AN','Org','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',35);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
%title({'Contribution of the chemical composition to the σ_{ext 525-520nm}'});
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_year_Mm(17NOV2021)');
% export_fig PM1_year_Mm(17NOV2021) -png -r600 -transparent
% savefig(namesauve);
%% 5.4 PM1 mass concentration
Var  = {'Org','NO3','NH4','SO4','Chl','eBC'};
labels = {'Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','Cl^{-}','eBC'};
PM_tot = [nansum(MAT_Month.Jan.PM1_ACSM); nansum(MAT_Month.Feb.PM1_ACSM); nansum(MAT_Month.Mar.PM1_ACSM); nansum(MAT_Month.Apr.PM1_ACSM);...
    nansum(MAT_Month.May.PM1_ACSM); nansum(MAT_Month.Jun.PM1_ACSM); nansum(MAT_Month.Jul.PM1_ACSM); nansum(MAT_Month.Aug.PM1_ACSM);...
    nansum(MAT_Month.Sep.PM1_ACSM); nansum(MAT_Month.Oct.PM1_ACSM); nansum(MAT_Month.Nov.PM1_ACSM); nansum(MAT_Month.Dec.PM1_ACSM)];
PM1_M = [MAT_Month.Jan.PM1_ACSM; MAT_Month.Feb.PM1_ACSM; MAT_Month.Mar.PM1_ACSM; MAT_Month.Apr.PM1_ACSM;...
    MAT_Month.May.PM1_ACSM; MAT_Month.Jun.PM1_ACSM; MAT_Month.Jul.PM1_ACSM; MAT_Month.Aug.PM1_ACSM;...
    MAT_Month.Sep.PM1_ACSM; MAT_Month.Oct.PM1_ACSM; MAT_Month.Nov.PM1_ACSM; MAT_Month.Dec.PM1_ACSM];

PM1_sdt = [MAT_Month.Jan.PM1_ACSM_std; MAT_Month.Feb.PM1_ACSM_std; MAT_Month.Mar.PM1_ACSM_std; MAT_Month.Apr.PM1_ACSM_std;...
    MAT_Month.May.PM1_ACSM_std; MAT_Month.Jun.PM1_ACSM_std; MAT_Month.Jul.PM1_ACSM_std; MAT_Month.Aug.PM1_ACSM_std;...
    MAT_Month.Sep.PM1_ACSM_std; MAT_Month.Oct.PM1_ACSM_std; MAT_Month.Nov.PM1_ACSM_std; MAT_Month.Dec.PM1_ACSM_std];

clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,150,100])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'FontSize',18);
box on;
area(1:12,PM1_M);
% hold on
% plot(1:12,PM_tot,'--','Color',[0.2 0.2 0.2],'LineWidth',8)
xlim([1 12]);
newcolors = [0.200000002980232 0.800000011920929 0;...%Org
    0 0.600000023841858 1;...%NO3
    0.929411768913269 0.694117665290833 0.125490203499794;%NH4
    0.9 0 0;...%SO4
    1 0 1;...%Cl
    0 0 0];%eBC
colororder(newcolors)
ylabel('PM_1 [µg/m^3]');
xlabel('Months')
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
%xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});box on;
legend(labels,'Location','northoutside','Orientation','horizontal','FontSize',35);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
%title({'Contribution of the chemical composition to the σ_{ext 525-520nm}'});
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_year_ACSM(17NOV2021)');
% export_fig PM1_year_ACSM(17NOV2021) -png -r600 -transparent
% savefig(namesauve);
%% 5.5 CONTRIBUTION TO THE SSA
Var  = {'NH42SO4' 'NH4NO3' 'Org' 'BC6'};
[range_final,freq,data_mean,data_std,ind_test1] = frequency([0.3:0.1:1],SSA_G);
h=[0.3:0.1:1];%--PERFECT
for ii=1:max((size(h)-1))
    %idx=find(h(ii)>=MAT.SSA_G & MAT.SSA_G<=h(ii+1));%%PERFECT
    idx=find(SSA_G>=h(ii) & SSA_G<=h(ii+1));%%NOT PERFECT
    for j=1:max(size(Var))
        MAT_SSA.(Var{j})(ii)=100.*nanmean(MAT.(Var{j})(idx)./MAT.PM1_calc(idx));
    end
end

X=[0.3 0.4 0.5 0.6 0.7 0.8 0.9];
species = [1:4];
size_SSA = [101; 408; 1006; 2722; 5666; 6637; 1579];
%size_SSA = [127; 457; 1152; 2855; 5782; 6340; 1410];
for j=1:max(size(Var))%-1
        P_PM1_calc(:,j) = MAT_SSA.(Var{j});
end
%%%PLOT
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(gcf,'units','centimeters','position',[0,0,100,100])
% Create multiple lines using matrix input to bar
yyaxis left
bar1 = bar(X,P_PM1_calc,'BarLayout','stacked','Parent',axes1);
set(bar1(1),'DisplayName','Ammonium sulfate',...
    'FaceColor',[0.9 0 0]);%1 0 0 NH_{4}(2)SO_{4}
set(bar1(2),'DisplayName','Ammonium nitrate','FaceColor',[0 0.2 0.9]);
set(bar1(3),'DisplayName','Organics',...
    'FaceColor',[0.200000002980232 0.800000011920929 0]);
set(bar1(4),'DisplayName','eBC','FaceColor',[0 0 0]);
ylim([0 100]);
xlim([0.2 1])
ylabel('Contribution (%)')
yyaxis right
plot(X,size_SSA,'o','MarkerFaceColor','black','MarkerSize',18);box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
xlabel('SSA in PM_{1}');
ylabel('Number of points')
%annotation('textbox',[.9 .5 .1 .2],'String','Number of points','EdgeColor','none','Fontname','Garamond','FontSize',50);
legend({'Ammonium sulfate','Ammonium nitrate','Organics','eBC','#'},'Location','northoutside','Orientation','horizontal','FontSize',35);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','FontSize',35
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
%title({'Contribution of the chemical composition to the SSA'})
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_SSA(17NOV2021)');
% export_fig PM1_SSA(17NOV2021) -png -r600 -transparent
% savefig(namesauve);

return

%% PM1 and ext
%% a. Extinction by month %%Feb =>75% of the data: 1051*100/1390.3448(Feb total coverage)
Var  = {'ExtG','ExtG_MLR','ExtG_IMPROVE'};
monthofyear = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'} ;
%% b. Extinction measured
for i=1:max(size(monthofyear))
    new_meas.Ext_median(i) = nanmedian(MAT_Month.(monthofyear{i}).ExtG);
    new_meas.Ext_mean(i) = nanmean(MAT_Month.(monthofyear{i}).ExtG);
    new_meas.RH(i) = nanmedian(MAT_Month.(monthofyear{i}).RH);
    %new_2017.SSA_std(i) = prctile(MAT_Month_2017.(monthofyear{i}).SSA_G,prc);
    new_meas.p25(i) = prctile(MAT_Month.(monthofyear{i}).ExtG,25);
    new_meas.p75(i) = prctile(MAT_Month.(monthofyear{i}).ExtG,75);
end
labels = {'Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','Cl^{-}','eBC','Ext Coeff'};
clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,150,100])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'FontSize',18);
box on;
colororder({'k','k'})
yyaxis left
area(1:12,PM1_M);
%hold on
xlim([1 12]);
newcolors = [0.200000002980232 0.800000011920929 0;...%Org
    0 0.600000023841858 1;...%NO3
    0.929411768913269 0.694117665290833 0.125490203499794;%NH4
    0.9 0 0;...%SO4
    1 0 1;...%Cl
    0 0 0];%eBC
colororder(newcolors)
ylabel('PM_1 [ug/m^3]');
yyaxis right
plot(1:12,new_meas.Ext_median,'o','Color',[0.2 0.2 0.2],'LineWidth',8)%,'MarkerEdgeColor','white');
ylim([0 70])
ylabel('Mm^{-1}');
xlabel('Months')
xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
%xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});box on;
legend(labels,'Location','northoutside','Orientation','horizontal','FontSize',35);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
%title({'Contribution of the chemical composition to the σ_{ext 525-520nm}'});
set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\600dp_resolution\PM1_year_ACSM_withEXT(17NOV2021)');
% export_fig PM1_year_ACSM_withEXT(17NOV2021) -png -r600 -transparent
% savefig(namesauve);

return

%% Time series NR-PM1 & Ext
%%%%%%%%%%% 1. Bar per month/per year
%MAT_Month = sort_month(MAT.Time,MAT);
MAT.PM1 = Org+NO3+NH4+SO4+BC6;
idx_2017=find(Time_ACSM>datenum(2017,8,1,0,0,0) & Time_ACSM<datenum(2018,1,1,0,0,0));
idx_2018=find(Time_ACSM>datenum(2018,1,1,0,0,0) & Time_ACSM<datenum(2019,1,1,0,0,0));
idx_2019=find(Time_ACSM>datenum(2019,1,1,0,0,0) & Time_ACSM<datenum(2020,1,1,0,0,0));
monthofyear = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'} ;
Annee = {'2017' '2018' '2019'};
Var_Annee = {'Y2017' 'Y2018' 'Y2019'};
Var = {'Time','Org','NO3','NH4','SO4','Chl','BC6','PM1'};
%%%%%% Splitting by year %%%%
for j=1:max(size(Annee))
        for ii =1:max(size(Var))
            MAT_years.(Var_Annee{j}).(Var{ii}) =[];
        end        
end

%for j=1:max(size(Annee))
%for i = 1:12
%Mois = num2str(i);
for ii =1:max(size(Var))
    %         MAT_years.(Var_Annee{j}).(Var{ii}) = [MAT_years.(Var_Annee{j}).(Var{ii});MAT.(Var{ii})];
    MAT_years.Y2017.(Var{ii}) = [MAT_years.Y2017.(Var{ii});MAT.(Var{ii})(idx_2017)];
    MAT_years.Y2018.(Var{ii}) = [MAT_years.Y2018.(Var{ii});MAT.(Var{ii})(idx_2018)];
    MAT_years.Y2019.(Var{ii}) = [MAT_years.Y2019.(Var{ii});MAT.(Var{ii})(idx_2019)];
end
%end
%end
%%%%%%%% Splitting by month %%%%%%%%%%%
MAT_Month_2017 = sort_month(MAT_years.Y2017.Time,MAT_years.Y2017);
MAT_Month_2018 = sort_month(MAT_years.Y2018.Time,MAT_years.Y2018);
MAT_Month_2019 = sort_month(MAT_years.Y2019.Time,MAT_years.Y2019);
%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%
%Xdata = [datenum(2017,8,1,0,0,0):datenum(0,1,0,0,0,0):datenum(2020,1,1,0,0,0)];
Data = [nanmean(MAT_Month_2017.Aug.Org) nanmean(MAT_Month_2017.Aug.NO3) nanmean(MAT_Month_2017.Aug.NH4) nanmean(MAT_Month_2017.Aug.SO4) nanmean(MAT_Month_2017.Aug.BC6);...
    nanmean(MAT_Month_2017.Sep.Org) nanmean(MAT_Month_2017.Sep.NO3) nanmean(MAT_Month_2017.Sep.NH4) nanmean(MAT_Month_2017.Sep.SO4) nanmean(MAT_Month_2017.Sep.BC6);...
    nanmean(MAT_Month_2017.Oct.Org) nanmean(MAT_Month_2017.Oct.NO3) nanmean(MAT_Month_2017.Oct.NH4) nanmean(MAT_Month_2017.Oct.SO4) nanmean(MAT_Month_2017.Oct.BC6);...
    nanmean(MAT_Month_2017.Nov.Org) nanmean(MAT_Month_2017.Nov.NO3) nanmean(MAT_Month_2017.Nov.NH4) nanmean(MAT_Month_2017.Nov.SO4) nanmean(MAT_Month_2017.Nov.BC6);...
    nanmean(MAT_Month_2017.Dec.Org) nanmean(MAT_Month_2017.Dec.NO3) nanmean(MAT_Month_2017.Dec.NH4) nanmean(MAT_Month_2017.Dec.SO4) nanmean(MAT_Month_2017.Dec.BC6);...
    nanmean(MAT_Month_2018.Jan.Org) nanmean(MAT_Month_2018.Jan.NO3) nanmean(MAT_Month_2018.Jan.NH4) nanmean(MAT_Month_2018.Jan.SO4) nanmean(MAT_Month_2018.Jan.BC6);...
    nanmean(MAT_Month_2018.Feb.Org) nanmean(MAT_Month_2018.Feb.NO3) nanmean(MAT_Month_2018.Feb.NH4) nanmean(MAT_Month_2018.Feb.SO4) nanmean(MAT_Month_2018.Feb.BC6);...
    nanmean(MAT_Month_2018.Mar.Org) nanmean(MAT_Month_2018.Mar.NO3) nanmean(MAT_Month_2018.Mar.NH4) nanmean(MAT_Month_2018.Mar.SO4) nanmean(MAT_Month_2018.Mar.BC6);...
    nanmean(MAT_Month_2018.Apr.Org) nanmean(MAT_Month_2018.Apr.NO3) nanmean(MAT_Month_2018.Apr.NH4) nanmean(MAT_Month_2018.Apr.SO4) nanmean(MAT_Month_2018.Apr.BC6);...
    nanmean(MAT_Month_2018.May.Org) nanmean(MAT_Month_2018.May.NO3) nanmean(MAT_Month_2018.May.NH4) nanmean(MAT_Month_2018.May.SO4) nanmean(MAT_Month_2018.May.BC6);...
    nanmean(MAT_Month_2018.Jun.Org) nanmean(MAT_Month_2018.Jun.NO3) nanmean(MAT_Month_2018.Jun.NH4) nanmean(MAT_Month_2018.Jun.SO4) nanmean(MAT_Month_2018.Jun.BC6);...
    nanmean(MAT_Month_2018.Jul.Org) nanmean(MAT_Month_2018.Jul.NO3) nanmean(MAT_Month_2018.Jul.NH4) nanmean(MAT_Month_2018.Jul.SO4) nanmean(MAT_Month_2018.Jul.BC6);...
    nanmean(MAT_Month_2018.Aug.Org) nanmean(MAT_Month_2018.Aug.NO3) nanmean(MAT_Month_2018.Aug.NH4) nanmean(MAT_Month_2018.Aug.SO4) nanmean(MAT_Month_2018.Aug.BC6);...
    nanmean(MAT_Month_2018.Sep.Org) nanmean(MAT_Month_2018.Sep.NO3) nanmean(MAT_Month_2018.Sep.NH4) nanmean(MAT_Month_2018.Sep.SO4) nanmean(MAT_Month_2018.Sep.BC6);...
    nanmean(MAT_Month_2018.Oct.Org) nanmean(MAT_Month_2018.Oct.NO3) nanmean(MAT_Month_2018.Oct.NH4) nanmean(MAT_Month_2018.Oct.SO4) nanmean(MAT_Month_2018.Oct.BC6);...
    nanmean(MAT_Month_2018.Nov.Org) nanmean(MAT_Month_2018.Nov.NO3) nanmean(MAT_Month_2018.Nov.NH4) nanmean(MAT_Month_2018.Nov.SO4) nanmean(MAT_Month_2018.Nov.BC6);...
    nanmean(MAT_Month_2018.Dec.Org) nanmean(MAT_Month_2018.Dec.NO3) nanmean(MAT_Month_2018.Dec.NH4) nanmean(MAT_Month_2018.Dec.SO4) nanmean(MAT_Month_2018.Dec.BC6);...
    nanmean(MAT_Month_2019.Jan.Org) nanmean(MAT_Month_2019.Jan.NO3) nanmean(MAT_Month_2019.Jan.NH4) nanmean(MAT_Month_2019.Jan.SO4) nanmean(MAT_Month_2019.Jan.BC6);...
    nanmean(MAT_Month_2019.Feb.Org) nanmean(MAT_Month_2019.Feb.NO3) nanmean(MAT_Month_2019.Feb.NH4) nanmean(MAT_Month_2019.Feb.SO4) nanmean(MAT_Month_2019.Feb.BC6);...
    nanmean(MAT_Month_2019.Mar.Org) nanmean(MAT_Month_2019.Mar.NO3) nanmean(MAT_Month_2019.Mar.NH4) nanmean(MAT_Month_2019.Mar.SO4) nanmean(MAT_Month_2019.Mar.BC6);...
    nanmean(MAT_Month_2019.Apr.Org) nanmean(MAT_Month_2019.Apr.NO3) nanmean(MAT_Month_2019.Apr.NH4) nanmean(MAT_Month_2019.Apr.SO4) nanmean(MAT_Month_2019.Apr.BC6);...
    nanmean(MAT_Month_2019.May.Org) nanmean(MAT_Month_2019.May.NO3) nanmean(MAT_Month_2019.May.NH4) nanmean(MAT_Month_2019.May.SO4) nanmean(MAT_Month_2019.May.BC6);...
    nanmean(MAT_Month_2019.Jun.Org) nanmean(MAT_Month_2019.Jun.NO3) nanmean(MAT_Month_2019.Jun.NH4) nanmean(MAT_Month_2019.Jun.SO4) nanmean(MAT_Month_2019.Jun.BC6);...
    nanmean(MAT_Month_2019.Jul.Org) nanmean(MAT_Month_2019.Jul.NO3) nanmean(MAT_Month_2019.Jul.NH4) nanmean(MAT_Month_2019.Jul.SO4) nanmean(MAT_Month_2019.Jul.BC6);...
    nanmean(MAT_Month_2019.Aug.Org) nanmean(MAT_Month_2019.Aug.NO3) nanmean(MAT_Month_2019.Aug.NH4) nanmean(MAT_Month_2019.Aug.SO4) nanmean(MAT_Month_2019.Aug.BC6);...
    nanmean(MAT_Month_2019.Sep.Org) nanmean(MAT_Month_2019.Sep.NO3) nanmean(MAT_Month_2019.Sep.NH4) nanmean(MAT_Month_2019.Sep.SO4) nanmean(MAT_Month_2019.Sep.BC6);...
    nanmean(MAT_Month_2019.Oct.Org) nanmean(MAT_Month_2019.Oct.NO3) nanmean(MAT_Month_2019.Oct.NH4) nanmean(MAT_Month_2019.Oct.SO4) nanmean(MAT_Month_2019.Oct.BC6);...
    nanmean(MAT_Month_2019.Nov.Org) nanmean(MAT_Month_2019.Nov.NO3) nanmean(MAT_Month_2019.Nov.NH4) nanmean(MAT_Month_2019.Nov.SO4) nanmean(MAT_Month_2019.Nov.BC6);...
    nanmean(MAT_Month_2019.Dec.Org) nanmean(MAT_Month_2019.Dec.NO3) nanmean(MAT_Month_2019.Dec.NH4) nanmean(MAT_Month_2019.Dec.SO4) nanmean(MAT_Month_2019.Dec.BC6)];
%x = [1:29]; 
PM1_montlhy = [nanmean(MAT_Month_2017.Aug.PM1);nanmean(MAT_Month_2017.Sep.PM1); nanmean(MAT_Month_2017.Oct.PM1);...
    nanmean(MAT_Month_2017.Nov.PM1); nanmean(MAT_Month_2017.Dec.PM1); nanmean(MAT_Month_2018.Jan.PM1);...
    nanmean(MAT_Month_2018.Feb.PM1); nanmean(MAT_Month_2018.Mar.PM1); nanmean(MAT_Month_2018.Apr.PM1);... 
    nanmean(MAT_Month_2018.May.PM1); nanmean(MAT_Month_2018.Jun.PM1); nanmean(MAT_Month_2018.Jul.PM1);... 
    nanmean(MAT_Month_2018.Aug.PM1); nanmean(MAT_Month_2018.Sep.PM1); nanmean(MAT_Month_2018.Oct.PM1);... 
    nanmean(MAT_Month_2018.Nov.PM1); nanmean(MAT_Month_2018.Dec.PM1); nanmean(MAT_Month_2019.Jan.PM1);...
    nanmean(MAT_Month_2019.Feb.PM1); nanmean(MAT_Month_2019.Mar.PM1); nanmean(MAT_Month_2019.Apr.PM1);...
    nanmean(MAT_Month_2019.May.PM1); nanmean(MAT_Month_2019.Jun.PM1); nanmean(MAT_Month_2019.Jul.PM1);... 
    nanmean(MAT_Month_2019.Aug.PM1); nanmean(MAT_Month_2019.Sep.PM1); nanmean(MAT_Month_2019.Oct.PM1);...
    nanmean(MAT_Month_2019.Nov.PM1); nanmean(MAT_Month_2019.Dec.PM1)];

x = categorical({'Aug/17','Sep/17','Oct/17','Nov/17','Dec/17','Jan/18','Feb/18','Mar/18',...
'Apr/18','May/18','Jun/18','Jul/18','Aug/18','Sep/18','Oct/18','Nov/18','Dec/18',...
'Jan/19','Feb/19','Mar/19','Apr/19','May/19','Jun/19','Jul/19','Aug/19','Sep/19',...
'Oct/19','Nov/19','Dec/19'});
x = reordercats(x,{'Aug/17','Sep/17','Oct/17','Nov/17','Dec/17','Jan/18','Feb/18','Mar/18',...
'Apr/18','May/18','Jun/18','Jul/18','Aug/18','Sep/18','Oct/18','Nov/18','Dec/18',...
'Jan/19','Feb/19','Mar/19','Apr/19','May/19','Jun/19','Jul/19','Aug/19','Sep/19',...
'Oct/19','Nov/19','Dec/19'});

clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,150,100])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'FontSize',18);
box on;
bar1 = bar(x,Data,'stacked');
lgd=legend({'Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',28)
%lgd.NumColumns = 3
set(bar1(1), 'FaceColor',[0.200000002980232 0.800000011920929 0]);
set(bar1(2), 'FaceColor',[0 0.600000023841858 1]);
set(bar1(4), 'FaceColor',[0.800000011920929 0 0]);
set(bar1(3), 'FaceColor',[0.929411768913269 0.694117665290833 0.125490203499794]);
set(bar1(5), 'FaceColor',[0 0 0]);
%xticklabels({'A/2017' 'S' 'O' 'N' 'J/2018' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D' 'J/2019' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
%xticklabels({'Aug/2017' 'Dec' 'Jun/2018' 'Oct' 'Mar/2019' 'Aug' 'Jan/2020'});box on;
ylabel('Concentration (µg m^3)','FontSize',35, 'FontWeight','Demi','FontName','times')
set(gca,'FontSize',20, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% hold on
% plot(x,PM1_montlhy)
text(x(1),11,{'910'})%Aug
text(x(2),11.5,{'1028'})%Sep
text(x(3),9.2,{'1024'})%Oct
text(x(4),18,{'1063'})%NOv
text(x(12),14,{'371'})%Jul/2018
text(x(13),6.4,{'801'})%Aug
text(x(14),10.2,{'1089'})%SEp
text(x(15),11.4,{'566'})%Oct
text(x(16),18.5,{'798'})%Nov
text(x(17),13,{'1073'})%Dec
text(x(18),12.6,{'248'})%Jan/2019
text(x(19),15.3,{'920'})%FEb
text(x(20),10,{'815'})%Mar
text(x(21),17.7,{'919'})%APr
text(x(22),10.5,{'1129'})%May
text(x(23),8.6,{'1034'})%Jun
text(x(24),6.3,{'1160'})%Jul
text(x(25),7.45,{'949'})%Aug
text(x(26),4.6,{'668'})%sEP
text(x(27),6.45,{'732'})%Oct
text(x(28),8.6,{'928'})%Nov
text(x(29),10.5,{'768'})
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\Time_series_bar3');
% export_fig Time_series_bar3 -png -r600 -transparent
% savefig(namesauve);

%Xdata = find(month(Time_ACSM)>=1&month(Time_ACSM)<12);
%%%%%%%%%%% 2. Moving average %%%%%%%%%%%%%%%%%%
% Smooth_PM1 = [Org NO3 NH4 SO4 Chl BC6];%PM1
% Smooth_AOPs = [EXT_525 Scat_G Abs_BC3];
% %%AOPs
% for I = 1:3
%     %AOPs(:,I) = smooth(Smooth_AOPs(:,I),'rlowess');
%     %AOPs(:,I) = smooth(Smooth_AOPs(:,I),15,'sgolay',4);%(y,span,'sgolay',degree)
%     AOPs(:,I) = moving_average(Smooth_AOPs(:,I),150);
% end
% %%PM1
% for I = 1:6
%     %CComp(:,I) = smooth(Smooth_PM1(:,I),'rlowess');
%     %CComp(:,I) = smooth(Smooth_PM1(:,I),10,'sgolay',4);
%     CComp(:,I) = moving_average(Smooth_PM1(:,I),150);
%     CCComp(:,I) = moving_average(Smooth_PM1(:,I),300);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Chemical Composition %%%%%%%%%%%%%%%%%%%%%%%
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,50,15])%0,0,150,150
% % Create axes
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'Fontname','Garamond',...
%     'FontSize',32);
% hold(axes1,'on');box on;
% %colororder({'k','k'})
% %%%%Subplot 1
% %s1=subplot(2,1,1)
% %yyaxis right
% %plot(Time_ACSM,AOPs(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
% %ylabel('Mm^{-1}')
% %plot(Time_ACSM,CComp(:,1),'--','Color',[0.5 0.6 0.6],'LineWidth',1.5);%1 0.5 0.5
% %hold on
% %yyaxis left
% plot(Time_ACSM,sum(CComp(:,2:end),2),'-','Color',[0.200000002980232 0.800000011920929 0],'LineWidth',1.5);
% hold on
% plot(Time_ACSM,sum(CComp(:,3:end),2),'-','Color',[0 0.600000023841858 1],'LineWidth',1.5);
% plot(Time_ACSM,sum(CComp(:,4:end),2),'-','Color',[0.929411768913269 0.694117665290833 0.125490203499794],'LineWidth',1.5);
% plot(Time_ACSM,sum(CComp(:,5:end),2),'-','Color',[0.9 0 0],'LineWidth',1.5);
% %plot(Time_ACSM,CComp(:,6),'-','Color',[1 0 1],'LineWidth',2.5);
% plot(Time_ACSM,sum(CComp(:,7:end),2),'-','Color',[0 0 0],'LineWidth',1.5);
% ylabel('µg m^3')
% %lgd=legend({'PM_1','Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','eBC'},'Location','best','Orientation','horizontal')
% lgd=legend({'Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','eBC'},'Location','best','Orientation','horizontal')
% %lgd=legend({'PM_1','Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','Cl^{-}','eBC'},'Location','best','Orientation','horizontal')
% set(gca,'FontSize',21, 'FontWeight','Demi','FontName','times','LineWidth',2.0);
% datetick('x','mmm/yyyy','keepticks')
% ylim([0 45])
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\Time_seriesCC');
% export_fig Time_seriesCC -png -r600 -transparent
% savefig(namesauve);
% %%% Try 1 %%%%%
% figure
% area(Time_ACSM,CComp(:,2:end))
% shading flat
% hold on
% %%%% Try 2 %%%%%
% Aver_Smooth = [];
% for i=1:5
% [Z_downsamp,t_downsamp] = downsample_ts(Smooth_PM1(:,i),Time_ACSM,'weekly','nanmean')
% Aver_Smooth = [Aver_Smooth Z_downsamp' ];
% end
% figure;
% hold on;
% plot(t_downsamp,Aver_Smooth(:,1),'g')
% plot(t_downsamp,Aver_Smooth(:,2),'b')
% plot(Time_ACSM,sum(CComp(:,2:end),2),'-','Color',[0.200000002980232 0.800000011920929 0],'LineWidth',1.5);
% plot(Time_ACSM,sum(CComp(:,3:end),2),'-','Color',[0 0.600000023841858 1],'LineWidth',1.5);
% figure
% plot(Time_ACSM,CCComp(:,1),'-','Color',[0.200000002980232 0.800000011920929 0],'LineWidth',1.5);
% hold on
% plot(Time_ACSM,CCComp(:,2),'-','Color',[0 0.600000023841858 1],'LineWidth',1.5);
% plot(Time_ACSM,CCComp(:,3),'-','Color',[0.929411768913269 0.694117665290833 0.125490203499794],'LineWidth',1.5);
% plot(Time_ACSM,CCComp(:,4),'-','Color',[0.9 0 0],'LineWidth',1.5);
% %plot(Time_ACSM,CComp(:,6),'-','Color',[1 0 1],'LineWidth',2.5);
% plot(Time_ACSM,CCComp(:,6),'-','Color',[0 0 0],'LineWidth',1.5);
% ylabel('µg m^3')
% %lgd=legend({'PM_1','Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','eBC'},'Location','best','Orientation','horizontal')
% lgd=legend({'Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','eBC'},'Location','best','Orientation','horizontal')
% %lgd=legend({'PM_1','Org','NO_{3}^{-}','NH_{4}^{+}','SO_{4}^{2-}','Cl^{-}','eBC'},'Location','best','Orientation','horizontal')
% set(gca,'FontSize',21, 'FontWeight','Demi','FontName','times','LineWidth',2.0);
% datetick('x','mmm/yyyy','keepticks')
% ylim([0 45])
% figure
% area(Time_ACSM,CCComp)
%%%%%%%%%%% AOPs %%%%%%%%%%%%%%%%%%%%%%%
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,50,15])%0,0,150,150
% % Create axes
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'Fontname','Garamond',...
%     'FontSize',32);
% hold(axes1,'on');box on;
% % %%%%Subplot 2 
% % s1=subplot(2,1,2)
% plot(Time_ACSM,AOPs(:,1),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);%EXT
% hold on
% plot(Time_ACSM,AOPs(:,2),'-','Color',[0 1 0],'LineWidth',1.5);
% plot(Time_ACSM,AOPs(:,3),'-','Color',[0 0 0],'LineWidth',1.5);
% ylabel('Mm^{-1}')
% lgd=legend({'Coeff Ext_{525nm}','Coeff Scat_{525nm}','Coeff Abs_{520nm}'},'Location','best','Orientation','vertical')
% dynamicDateTicks([],[],'mm/yy')
% lgd.NumColumns = 4;
% set(gca,'FontSize',25, 'FontWeight','Demi','FontName','times','LineWidth',2.0);
% % namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\Time_seriesAOPs');
% % export_fig Time_seriesAOPs -png -r600 -transparent
% % savefig(namesauve);

return
%% %%%%%%%%%%%%%%%%%%%% Another useful plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SSA vs Ratio
% idx_MAM=find(month(Time_ACSM)>=3&month(Time_ACSM)<6);%%spring
% Ratio_MLR(idx_MAM) = NaN;
% SSA_G(idx_MAM) = NaN;
% Ratio_IMPROVE(idx_MAM) = NaN;
% figure;
% plot(SSA_G,Ratio_MLR,'.k')
% ylabel('Ratio σ_{ext}Measured/_{ext}Retrieved [525nm]');
% xlabel('SSA_{PM1} [525nm]');
% set(gca,'FontSize',25, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% 
% 
% figure;
% plot(SSA_G,Ratio_IMPROVE,'.k')
% ylabel('Ratio σ_{ext}Measured/_{ext}Retrieved [525nm]');
% xlabel('SSA_{PM1} [525nm]');
% set(gca,'FontSize',25, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% 
% figure;
% X_axis = [SSA_G(:)];
% Y_axis = [Ratio_IMPROVE(:)];
% Z =[Time_ACSM-Time_ACSM(1)];
% %Z =[Temperature(:)];
% C=Z;
% scatter3(X_axis,Y_axis,Z,40,Z,'filled')
% view([0 90]);
% %ylim([0 100])
% % xlim([0 300])
% colorbar
% caxis([2 35])
%% Ext measured, retrived:MLR, IMPROVE
% Var  = {'ExtG','ExtG_MLR','ExtG_IMPROVE'};
% labels = {'Measured','MLR','IMPROVE algorithm'};
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,150,100])
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'Fontname','Garamond',...
%     'FontSize',25);
% box on;
% %axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% plot([0.6:11.6],new_meas.Ext_median,'k');
% p1=errorbarxy([0.6:11.6],new_meas.Ext_median,[0.6:11.6].*0,(new_meas.p75-new_meas.Ext_median),[0.6:11.6].*0,(new_meas.Ext_median-new_meas.p25),...
%     'Color',[0 0 0],'LineStyle','-','Marker','s','MarkerFaceColor',[0 0 0],'LineWidth',1,'MarkerSize',9);
% %% b. Extinction MLR
% for i=1:max(size(monthofyear))
%     new_calc.Ext_median(i) = nanmedian(MAT_Month.(monthofyear{i}).ExtG_MLR);
%     %new_2017.SSA_std(i) = prctile(MAT_Month_2017.(monthofyear{i}).SSA_G,prc);
%     new_calc.p25(i) = prctile(MAT_Month.(monthofyear{i}).ExtG_MLR,25);
%     new_calc.p75(i) = prctile(MAT_Month.(monthofyear{i}).ExtG_MLR,75);
% end
% plot([0.7:11.7],new_calc.Ext_median,'r');
% p2=errorbarxy([0.7:11.7],new_calc.Ext_median,[0.7:11.7].*0,(new_calc.p75-new_calc.Ext_median),[0.7:11.7].*0,(new_calc.Ext_median-new_calc.p25),...
%     'Color',[0.9 0 0],'LineStyle','-','Marker','*','MarkerFaceColor',[0.9 0 0],'LineWidth',1,'MarkerSize',9);
% %% c. EXtinction IMPROVE
% for i=1:max(size(monthofyear))
%     new_calc2.Ext_median(i) = nanmedian(MAT_Month.(monthofyear{i}).ExtG_IMPROVE);
%     %new_2017.SSA_std(i) = prctile(MAT_Month_2017.(monthofyear{i}).SSA_G,prc);
%     new_calc2.p25(i) = prctile(MAT_Month.(monthofyear{i}).ExtG_IMPROVE,25);
%     new_calc2.p75(i) = prctile(MAT_Month.(monthofyear{i}).ExtG_IMPROVE,75);
% end
% plot([0.8:11.8],new_calc2.Ext_median,'b');
% p3=errorbarxy([0.8:11.8],new_calc2.Ext_median,[0.8:11.8].*0,(new_calc2.p75-new_calc2.Ext_median),[0.8:11.8].*0,(new_calc2.Ext_median-new_calc2.p25),...
%     'Color',[0 0 1],'LineStyle','-','Marker','o','MarkerFaceColor',[0 0 1],'LineWidth',1,'MarkerSize',9);
% dynamicDateTicks([],[],'dd');
% xlim([min(0) max(12.5)]);
% ylim([0 160]);
% %title('Coeff Extinction_{520-525 nm}');
% ylabel('Coeff Ext (Mm-1)');
% %xlabel('Months');
% box on
% legend([p1(1) p2(1) p3(1)],{'Measured','MLR','IMPROVE'},'Location','bestoutside','FontSize',45);%northeastoutside
% set(axes1,'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13],'XTickLabel',...
%         {'','J','F','M','A','M','J','J','A','S','O','N','D',''});
% %     {'','Jan','Feb','Mar','Abr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan'});
% % t=annotation(figure1,'textbox',...
% %     [0.806299014162669 0.629529837251361 0.155506004628151 0.0411392405063251],...
% %     'String','Median and percentil (75-25)',...
% %     'FitBoxToText','off',...
% %     'EdgeColor','none');
% set(gca,'Fontname','Garamond','Fontsize',50,'linewidth',3);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\Ext_525');
% export_fig Ext_rangeplots -png -r300 -transparent
% savefig(namesauve);
%% Time series retrival
% yPM1 = smooth(Time_ACSM,PM1,0.1,'rloess');
% yEXT = smooth(Time_ACSM,EXT_525,0.1,'rloess');
% yOrg = smooth(Time_ACSM,Org,0.1,'rloess');
% yNO3 = smooth(Time_ACSM,NO3,0.1,'rloess');
% yNH4 = smooth(Time_ACSM,NH4,0.1,'rloess');
% ySO4 = smooth(Time_ACSM,SO4,0.1,'rloess');
% yChl = smooth(Time_ACSM,Chl,0.1,'rloess');
% yeBC = smooth(Time_ACSM,BC6,0.1,'rloess');
% % Crate plot
% colororder({'k','k'})
% yyaxis left
% plot(Time_ACSM,yPM1,'--','Color',[1 0.5 0.5],'LineWidth',2.5);
% hold on
% plot(Time_ACSM,yOrg,'-','Color',[0.200000002980232 0.800000011920929 0],'LineWidth',2.5);
% plot(Time_ACSM,yNO3,'-','Color',[0 0.600000023841858 1],'LineWidth',2.5);
% plot(Time_ACSM,yNH4,'-','Color',[0.929411768913269 0.694117665290833 0.125490203499794],'LineWidth',2.5);
% plot(Time_ACSM,ySO4,'-','Color',[0.9 0 0],'LineWidth',2.5);
% plot(Time_ACSM,yChl,'-','Color',[1 0 1],'LineWidth',2.5);
% plot(Time_ACSM,yeBC,'-','Color',[0 0 0],'LineWidth',2.5);
% ylabel('µg m^3')
% yyaxis right
% plot(Time_ACSM,yEXT,'--','Color',[0.5 0.5 0.5],'LineWidth',2.5);
% ylabel('Mm^{-1}')
%% bext (Mm-1) by IMPROVE, MLR, Measured
% y1 = smooth(Time_ACSM,IMPROVE_G,0.1,'rloess');
% y2 = smooth(Time_ACSM,Ext_retG,0.1,'rloess');
% y3 = smooth(Time_ACSM,EXT_525,0.1,'rloess');
% Create figure
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,150,150])
% % Create axes
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'Fontname','Garamond',...
%     'FontSize',32);
% hold(axes1,'on');box on;
% % Crate plot
% %plot(Time_ACSM(find(isgood==1)),y1(find(isgood==1)),'-','Color',[0 0 1],'LineWidth',2.5);
% plot(Time_ACSM(find(isgood==1)),IMPROVE_G(find(isgood==1)),'-','Color',[0 0 1],'LineWidth',2.5);
% hold on
% % plot(Time_ACSM(find(isgood==1)),y2(find(isgood==1)),'-','Color',[1 0 0],'LineWidth',1.0);
% % plot(Time_ACSM(find(isgood==1)),y3(find(isgood==1)),'--','Color',[0 0 0],'LineWidth',1.0);
% plot(Time_ACSM(find(isgood==1)),Ext_retG(find(isgood==1)),'-','Color',[1 0 0],'LineWidth',1.0);
% plot(Time_ACSM(find(isgood==1)),EXT_525(find(isgood==1)),'.k');
% ylim([0 800])
% dynamicDateTicks([],[],'mm/dd')
% ylabel('Coeff Ext (Mm-1)');
% xlabel('UTC Time');
% legend({'IMPROVE','MLR','Measured'},'Location','bestoutside','FontSize',45);
% %set(gca,'FontSize',50, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% %set(gca,'Fontname','Garamond','Fontsize',50,'linewidth',3);
% set(gca,'FontSize',50, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% yyaxis right
% plot(Time_ACSM(find(isgood==1)),PM1(find(isgood==1)),'--g');
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\time_series_Ext_525');
% export_fig Ext_Timeserie -png -r300 -transparent
% savefig(namesauve);
%% EXT and RH
%clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,150,100])
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'FontSize',18);
% box on;
% yyaxis left
% area(1:12,Cont_Ext_M);
% ylim([0 95]);
% xlim([1 12]);
% newcolors = [0.9 0 0;0 0.2 0.9;0.200000002980232 0.800000011920929 0;0 0 0];
% colororder(newcolors)
% ylabel('Ext Coeff (Mm-1)');
% xlabel('Months')
% xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
% %xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});box on;
% yyaxis right
% plot(1:12,new_meas.RH,'o','Color',[0.2 0.2 0.2],'LineWidth',8)
% legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',50);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
% %title({'Contribution of the chemical composition to the σ_{ext 525-520nm}'});
% set(gca,'FontSize',50, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
%% MLR and EXT
figure
for i=1:10:200
hold on
circ(i,nanmedian(EXT_525(EXT_525>i)./Ext_retG(EXT_525>i)))
end
figure
for i=1:10:200
hold on
circ(i,size(find(EXT_525>i),1))
end
% nanmedian(EXT_525(EXT_525>100)./ExtG_Wang(EXT_525>100))
% nanmedian(EXT_525(EXT_525>100)./ExtG_Yao(EXT_525>100))
% nanmedian(EXT_525(EXT_525>100)./ExtG_Gro(EXT_525>100))
% nanmedian(EXT_525(EXT_525>100)./Ext_retG(EXT_525>100))
%% BOX PLOTS Literature
% Wang = nanmedian(EXT_525(EXT_525>100)./ExtG_Wang(EXT_525>100));%100
% Yao = nanmedian(EXT_525(EXT_525>100)./ExtG_Yao(EXT_525>100));
% Gro = nanmedian(EXT_525(EXT_525>100)./ExtG_Gro(EXT_525>100));
% MLR = nanmedian(EXT_525(EXT_525>100)./Ext_retG(EXT_525>100));
% IMPROVEG = nanmedian(EXT_525(EXT_525>100)./IMPROVE_G(EXT_525>100));
% Ratios = [MLR IMPROVEG Gro Yao Wang];
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,150,150])
% % Create axes
% axes1 = axes('Parent',fig1,...
%     'LineWidth',3,...
%     'FontWeight','Demi','FontName','times',...
%     'FontSize',32);%'Fontname','Garamond',...
% Authors = categorical({'MLR','Pitchoford et al. 2007','Groblicki et al. 1981','Yao et al. 2010','Wang et al. 2015'});
% Authors = reordercats(Authors,{'MLR','Pitchoford et al. 2007','Groblicki et al. 1981','Yao et al. 2010','Wang et al. 2015'});
% b=bar(Authors,Ratios,'FaceColor','black')
% %b.EdgeColor = [0 0 0];
% % set(bar1(1),'DisplayName','MLR',...
% %     'FaceColor',[0 0 1]);
% % set(bar1(2),'DisplayName','Pitchoford et al. 2007','FaceColor',[1 0 0]);
% % set(bar1(3),'DisplayName','Groblicki et al. 1981',...
% %     'FaceColor',[0 1 0]);
% % set(bar1(4),'DisplayName','Yao et al. 2010',...
% %     'FaceColor',[1 1 0]);
% % set(bar1(5),'DisplayName','Wang et al. 2015',...
% %     'FaceColor',[1 0 1]);
% %set(bar1(4),'DisplayName','eBC','FaceColor',[0 0 0]);
% %legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northeastoutside','FontSize',15);
% %legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org'},'Location','northoutside','Orientation','horizontal','FontSize',50);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','FontSize',35
% %title({'Contribution of the chemical composition to the σ_{scatt}'});
% xticklabels({'MLR','Pitchoford et al. 2007','Groblicki et al. 1981','Yao et al. 2010','Wang et al. 2015'});box on;
% %ylabel('m^{2}g^{-1}');ylabel('µg^-3');
% ylabel('Ratio Observed / Retrieved');%xlabel('Wavelength');
% set(gca,'FontSize',35, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% %boxplot(Ratios,{'MLR','Yao et al. 2010','Wang et al. 2015','Groblicki et al. 1981'})
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\EBAS_database\Boxplot_Ext');
% export_fig OBoxplot_Ext -png -r300 -transparent
% savefig(namesauve);
% print( '-djpeg', '-r650', namesauve);
%% Categorias de ext
Var={'EXT_525','Ext_retG','IMPROVE_G','PM1','Time'};
MAT_Cat.PM1 = PM1;
MAT_Cat.EXT_525 = EXT_525;
MAT_Cat.Ext_retG = Ext_retG;
MAT_Cat.IMPROVE_G = IMPROVE_G;
MAT_Cat.Time = Time_ACSM;
idx_P = prctile(EXT_525,75);%63.7 Mm-1
idx_C = prctile(EXT_525,25);%20.3 Mm-1
%% I. Empty matrix 
for j=1:max(size(Var))
    MAT_clean.(Var{j})=[MAT_Cat.(Var{j})];%% PM1 > 25percetile size:2375
    MAT_polluted.(Var{j})=[MAT_Cat.(Var{j})];%% PM1 < 75pertile size:2375
    MAT_middle.(Var{j})=[MAT_Cat.(Var{j})];
end
%% II. Sorting out =usign ext
for j=1:max(size(Var))
    MAT_clean.(Var{j})(find(EXT_525>idx_C))=NaN;%%
    size(find(isnan(MAT_clean.EXT_525)==0));%4267
    %     MAT_clean.(Var{j})(find(PM1>Clean_limit_PM1))=NaN;%%
    %     size(find(isnan(MAT_clean.PM1)==0));%%Means the numbers where there is data
    MAT_polluted.(Var{j})(find(EXT_525<idx_P))=NaN;%%
    size(find(isnan(MAT_polluted.EXT_525)==0));%4267
    %     MAT_polluted.(Var{j})(find(PM1<Polluted_limit_PM1))=NaN;%%20<=
    %     size(find(isnan(MAT_polluted.PM1)==0));
    MAT_middle.(Var{j})(find( EXT_525<idx_C | EXT_525>idx_P))=NaN;%% 50percentile
    size(find(isnan(MAT_middle.EXT_525)==0));%8535
    %     MAT_middle.(Var{j})(find( PM1<Clean_limit_PM1 | PM1>Polluted_limit_PM1))=NaN;%% 50percentile
    %     size(find(isnan(MAT_middle.PM1)==0));
end
%% Removing NaNs in the concentration
isgoodA2=  ~(isnan(MAT_clean.EXT_525));
isgoodB2 =  ~(isnan(MAT_polluted.EXT_525));
isgoodC2 =  ~(isnan(MAT_middle.EXT_525));
idxA=find(isgoodA2==0);%0=>NAN 1=>Numero
idxB=find(isgoodB2==0);
idxC=find(isgoodC2==0);
for j=1:max(size(Var))
    MAT_clean.(Var{j})(idxA)=NaN; %PM1 = 179
    MAT_polluted.(Var{j})(idxB)=NaN;% PM1 = 3316
    MAT_middle.(Var{j})(idxC)=NaN;% PM1= 2889
end
figure;
yyaxis left
plot(MAT_clean.Time,MAT_clean.PM1,'*','Color',[0.6 1 0.6]);
hold on
plot(MAT_polluted.Time,MAT_polluted.PM1,'*','Color',[1, 0.75, 0.75]);
plot(MAT_middle.Time,MAT_middle.PM1,'*','Color',[0.7 0.7 0]);
ylabel('µg m^{-1}')
yyaxis right
plot(MAT_clean.Time,MAT_clean.EXT_525,'.g');
hold on
plot(MAT_polluted.Time,MAT_polluted.EXT_525,'.r');
plot(MAT_middle.Time,MAT_middle.EXT_525,'.','Color',[0.7 0.7 0]);
ylabel('Mm^{-1}')
legend('Clean','Polluted','Middle');
dynamicDateTicks([],[],'mm/dd');


x=MAT_polluted.EXT_525;
y=MAT_polluted.Ext_retG;
figure;
plot(x,y,'.k')
hold on
plot([0 700],[0 700],'--r')

x=MAT_polluted.EXT_525;
y=MAT_polluted.Ext_retG;

x=MAT_clean.EXT_525;
y=MAT_clean.Ext_retG;
%% Percentages 
%% SCATTERPLOT antiguo 
% Create figure
% % clear fig1
% % fig1=figure;
% % set(gcf,'units','centimeters','position',[0,0,150,150])
% % % Create axes
% % axes1 = axes('Parent',fig1,...
% %     'LineWidth',3,...
% %     'FontWeight','Demi','FontName','times',...
% %     'FontSize',32);%'Fontname','Garamond',...
% % hold(axes1,'on');box on;
% % %set(gca,'FontSize',50, 'FontWeight','Demi','FontName','times','LineWidth',3.0);
% % % Crate plot
% % plot([0 500],[0 500],'--k','LineWidth',2.5);
% % xlim([0 500])
% % ylim([0 500])
% % hold on
% % %%%%%%%%%%%%%%%%%%%%%%%%%% MLR %%%%%%%%%%%%%%%%
% % % Crate plot
% % p2=polyfit(EXT_525(find(isgood==1)),Ext_retG(find(isgood==1)),1);
% % plot2=plot(EXT_525(find(isgood==1)),Ext_retG(find(isgood==1)),'.','Color',[1 0.7 0.7]);
% % % Get xdata from plot
% % xdata1 = get(plot2, 'xdata');
% % % Get ydata from plot
% % ydata1 = get(plot2, 'ydata');
% % % Make sure data are column vectors
% % xdata1 = xdata1(:);
% % ydata1 = ydata1(:);
% % % Find x values for plotting the fit based on xlim
% % axesLimits1 = xlim(axes1);
% % xplot1 = linspace(axesLimits1(1), axesLimits1(2));
% % % Preallocate for "Show equations" coefficients
% % coeffs1 = cell(1,1);
% % % Find coefficients for polynomial (order = 1)
% % fitResults1 = polyfit(xdata1,ydata1,1);
% % % Evaluate polynomial
% % yplot1 = polyval(fitResults1,xplot1);
% % % Save type of fit for "Show equations"
% % fittypesArray1(1) = 2;
% % % Save coefficients for "Show Equation"
% % coeffs1{1} = fitResults1;
% % % Plot the fit
% % fitLine1 = plot(xplot1,yplot1,'DisplayName','linear','Tag','linear',...
% %     'Parent',axes1,...
% %     'LineWidth',2,...
% %     'Color',[1 0 0]);
% % % Set new line in proper position
% % setLineOrder(axes1,fitLine1,plot2);
% % ylabel('Ext coeff Retrieved (Mm^{-1})','FontWeight','Demi','FontName','times','Fontsize',48)
% % xlabel('Ext coeff Observed (Mm^{-1})','FontWeight','Demi','FontName','times','Fontsize',48)
% % mdl_MLR=fitlm(EXT_525,Ext_retG);%,'Intercept',false);%'RobustOpts','on');
% % mdl_IMPROVE=fitlm(EXT_525,IMPROVE_G);%,'Intercept',false);%'RobustOpts','on');
% % text(15,380,{'MLR:';'y = 0.91x + 21.61';'R^{2} = 0.78';'RMSE:17.40'},'Color',[0 0 0],'FontWeight','Demi','FontName','times','FontSize',45);%%before 0.97
% % %text(15,380,{'IMPROVE:';'y = 0.8x + 21.0';'R^{2} = 0.76';'RMSE:16.3'},'Color',[0 0 1],'FontWeight','Demi','FontName','times','FontSize',35);
% % % text(175,380,{'MLR:';'y = 0.92x + 22';'R^{2} = 0.80';'RMSE:21.32'},'Color','black','Fontname','Garamond','FontSize',35);%%before 0.97
% % % text(15,380,{'IMPROVE:';'y = 0.64x + 21';'R^{2} = 0.79';'RMSE:15.68'},'Color','black','Fontname','Garamond','FontSize',35);
% % text(250,80,{'# of data = 18 174'},'Color','black','FontWeight','Demi','FontName','times','FontSize',40);
% % 
%% Pollution roses
%% 3D plots
%% WD, FIno & SAE
clear fig1
fig1=figure;
set(gcf,'units','centimeters','position',[0,0,100,100])
axes1 = axes('Parent',fig1,...
    'LineWidth',2,...
    'FontSize',6);box on
%X_axis = [Time_ACSM];
X_axis = [SSA_G];
Y_axis = [EXT_525./Ext_retG];
%Y_axis = [PM1];
% %Z =[MAT.MAM.Time-MAT.MAM.Time(1)];
%Z =[mean_DPg];
Z =[Eff_radius];
C=Z;
scatter3(X_axis(find(isgood==1)),Y_axis(find(isgood==1)),Z(find(isgood==1)),40,Z(find(isgood==1)),'filled')
view([0 90]);grid on;box on
xlabel('Fraction of Inorganics','FontSize',18, 'FontWeight','Demi')
ylabel('Wind direction','FontSize',18, 'FontWeight','Demi');
legend('SAE 450-550nm');
title('Spring period')
set(gca,'Fontname','Garamond','Fontsize',15,'linewidth',2);
colorbar
% namesauve =strcat(['C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Link chemical-Optical properties\MLR_4species\Abs_Scat_Ext_MLR\Seasons\Ext\Spring_analysis\WDvsFInnovsSAE(dots_size)']);
% savefig(namesauve)
% print( '-djpeg', '-r650', namesauve)

%% SIA and Ext
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,150,100])
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'FontSize',18);
% box on;
% area(1:12,Cont_Ext_M);
% hold on;
% plot(1:12,new_meas.Ext_median,'o','Color',[0.2 0.2 0.2],'LineWidth',8)
% ylim([0 95]);
% xlim([1 12]);
% newcolors = [0.9 0 0;0 0.2 0.9;0.200000002980232 0.800000011920929 0;0 0 0];
% colororder(newcolors)
% ylabel('Ext Coeff (Mm-1)');
% xlabel('Months')
% xticklabels({'J' 'F' 'M' 'A' 'M' 'J' 'J' 'A' 'S' 'O' 'N' 'D'});box on;
% %xticklabels({'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});box on;
% legend({'(NH_{4})_{2} SO_{4}','NH_{4}NO_{3}','Org','eBC'},'Location','northoutside','Orientation','horizontal','FontSize',50);%Location','southoutside','Orientation','horizontal'%,'Location','northwestoutside','Fo
% %title({'Contribution of the chemical composition to the σ_{ext 525-520nm}'});
% set(gca,'FontSize',50, 'FontWeight','Demi','FontName','times','LineWidth',3.0);





% %% %%%%%%%%%%%%%%%%%%%%%%% 5.1 PM1 species fraction %%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%% 5.2 Pollution roses %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%% 5.3 Diurnal plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%% 5.4 Ext measured vs retrieved %%%%%%%%%%%%%%%%%%
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%% 5.5 Contribution of the chemical %%%%%%%%%%%%%%%%%
% %% %%%%%%%%%%%%%%%%%%%%%%% composition to the AOPs %%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%% 5.6 3D plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 5.6.1 Ext measured vs retrieved vs: SSA, FIno, FOrg, AAE370-520
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set(gca,'Fontname','Garamond','Fontsize',42,'linewidth',3);
% subplot1=subplot(2,1,1)
% xdata = EXT_525;
% ydata = Ext_retG;
% mdl_S525 = fitlm(xdata,ydata,'Intercept',false,'RobustOpts','on')
% R1=polyfitn(xdata(find(isgood==1)),ydata(find(isgood==1)),1)
% p2=polyfit(xdata(find(isgood==1)),ydata(find(isgood==1)),1);
% plot1=plot(xdata(find(isgood==1)),ydata(find(isgood==1)),'.k');
% hold on
% plot([0 500],[0 500],'--k','LineWidth',2.5);
% xlim([0 500])
% ylim([0 500])
% % Find x values for plotting the fit based on xlim
% axesLimits1 = xlim(subplot1);
% xplot1 = linspace(axesLimits1(1), axesLimits1(2));
% % Get xdata from plot
% xdata1 = get(plot1, 'xdata');
% % Get ydata from plot
% ydata1 = get(plot1, 'ydata');
% % Make sure data are column vectors
% xdata1 = xdata1(:);
% ydata1 = ydata1(:);
% % Preallocate for "Show equations" coefficients
% coeffs1 = cell(1,1);
% % Find coefficients for polynomial (order = 1)
% fitResults1 = polyfit(xdata1,ydata1,1);
% % Evaluate polynomial
% yplot1 = polyval(fitResults1,xplot1);
% % Save type of fit for "Show equations"
% fittypesArray1(1) = 2;
% % Save coefficients for "Show Equation"
% coeffs1{1} = fitResults1;
% % Plot the fit
% fitLine1 = plot(xplot1,yplot1,'DisplayName','   linear','Tag','linear',...
%     'Parent',subplot1,...
%     'LineWidth',2,...
%     'Color',[0 0 0]);
% % Set new line in proper position
% setLineOrder(subplot1,fitLine1,plot1);
% % % "Show equations" was selected
% % showEquations(fittypesArray1,coeffs1,2,subplot1);
% %plot(EXT_525(find(isgood_Ext==1)),polyval(p2,EXT_PM1_IMPROVE(find(isgood_Ext==1))),'LineWidth',2.5);
% legend({'Extinction_{520-525nm}','Linear fit','1:1'},'Location','northeastoutside','FontSize',25);% RMSE:0.83%'y=0.87x+0.08 r^{2}=0.72'
% %ylabel('σ_{ext} [Mm^{-1}] Retrieved','Fontname','Garamond','Fontsize',30)
% %xlabel('σ_{scatt} [Mm^{-1}] Measured','Fontname','Garamond','Fontsize',25)
% title({'MLR method'});
% % text(25,360,{'y = 0.97x + 22';'R^{2} = 0.80';'RMSE:22.52';'n=15686';'Bias: 20.3'},'Color','black','FontSize',14);
% text(15,370,{'y = 0.97x + 22';'R^{2} = 0.80'},'Color','black','FontSize',16);
% text(100,370,{'RMSE:22.52';'Bias: 20.3'},'Color','black','FontSize',16);
% text(180,370,{'n=15686'},'Color','black','FontSize',16);
% % t=annotation(fig1,'textbox',...
% %     [0.806299014162669 0.629529837251361 0.155506004628151 0.0411392405063251],...
% %     'String',{'R^{2} = 0.80';'RMSE:22.52';'n=15686';'Bias: 20.3'},...%FIno 'R^{2} = 0.80';'RMSE:16.7'
% %     'FitBoxToText','off',...
% %     'EdgeColor','none');%%0.71+20.2
% set(gca,'Fontname','Garamond','Fontsize',35,'linewidth',2);
% subplot2=subplot(2,1,2)
% xdata = EXT_525;
% ydata = IMPROVE_G;
% mdl_S525 = fitlm(xdata,ydata,'Intercept',false,'RobustOpts','on')
% R2=polyfitn(xdata(find(isgood==1)),ydata(find(isgood==1)),1)
% p2=polyfit(xdata(find(isgood==1)),ydata(find(isgood==1)),1);
% plot2=plot(xdata(find(isgood==1)),ydata(find(isgood==1)),'.k');
% hold on
% plot([0 500],[0 500],'--k','LineWidth',2.5);
% xlim([0 500])
% ylim([0 500])
% % Find x values for plotting the fit based on xlim
% axesLimits1 = xlim(subplot2);
% xplot1 = linspace(axesLimits1(1), axesLimits1(2));
% % Get xdata from plot
% xdata1 = get(plot2, 'xdata');
% % Get ydata from plot
% ydata1 = get(plot2, 'ydata');
% % Make sure data are column vectors
% xdata1 = xdata1(:);
% ydata1 = ydata1(:);
% % Preallocate for "Show equations" coefficients
% coeffs1 = cell(1,1);
% % Find coefficients for polynomial (order = 1)
% fitResults1 = polyfit(xdata1,ydata1,1);
% % Evaluate polynomial
% yplot1 = polyval(fitResults1,xplot1);
% % Save type of fit for "Show equations"
% fittypesArray1(1) = 2;
% % Save coefficients for "Show Equation"
% coeffs1{1} = fitResults1;
% % Plot the fit
% fitLine1 = plot(xplot1,yplot1,'DisplayName','   linear','Tag','linear',...
%     'Parent',subplot2,...
%     'LineWidth',2,...
%     'Color',[0 0 0]);
% %fitLine1.Fontsize=12;
% % Set new line in proper position
% setLineOrder(subplot2,fitLine1,plot2);
% % % "Show equations" was selected
% % showEquations(fittypesArray1,coeffs1,2,subplot2);
% legend({'Extinction_{520-525nm}','Linear fit','1:1'},'Location','northeastoutside','FontSize',25);% RMSE:0.83%'y=0.87x+0.08 r^{2}=0.72'
% ylabel('σ_{ext} [Mm^{-1}] Retrieved','Fontname','Garamond','Fontsize',80)
% xlabel('σ_{scatt} [Mm^{-1}] Measured','Fontname','Garamond','Fontsize',80)
% title({'IMPROVE algorithm'});
% %text(25,360,{'y = 0.65x + 21';'R^{2} = 0.79';'RMSE:15.68';'n=15686';'Bias: 3.13'},'Color','black','FontSize',14);
% text(15,370,{'y = 0.65x + 21';'R^{2} = 0.79'},'Color','black','FontSize',16);
% text(100,370,{'RMSE:15.68';'Bias: 3.13'},'Color','black','FontSize',16);
% text(180,370,{'n=15686'},'Color','black','FontSize',16);
% % t=annotation(fig1,'textbox',...
% %     [0.806299014162669 0.629529837251361 0.155506004628151 0.0411392405063251],...
% %     'String',{'R^{2} = 0.79';'RMSE:15.68';'n=15686';'Bias: 3.13'},...%FIno 'R^{2} = 0.80';'RMSE:16.7'
% %     'FitBoxToText','off',...
% %     'EdgeColor','none');%%0.71+20.2; Bias:nanmean(Y-X)
% set(gca,'Fontname','Garamond','Fontsize',35,'linewidth',2);
% namesauve =strcat('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\Plots\Papaer\300dp_resolution\IMPROVEvsMLR_525');
% export_fig IMPROVEvsMLR_G -png -r300 -transparent
% savefig(namesauve);
%print( '-djpeg', '-r650', namesauve);
%% a. Extinction measured
% for i=1:max(size(monthofyear))
%     new_meas.median(i) = nanmedian(MAT_Month.(monthofyear{i}).Org);
%     %new_2017.SSA_std(i) = prctile(MAT_Month_2017.(monthofyear{i}).SSA_G,prc);
%     new_meas.p25(i) = prctile(MAT_Month.(monthofyear{i}).Org,25);
%     new_meas.p75(i) = prctile(MAT_Month.(monthofyear{i}).Org,75);
% end
% clear fig1
% fig1=figure;
% set(gcf,'units','centimeters','position',[0,0,150,100])
% axes1 = axes('Parent',fig1,...
%     'LineWidth',2,...
%     'Fontname','Garamond',...
%     'FontSize',25);
% box on;
% %axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% plot([0.6:11.6],new_meas.median,'k');
% p1=errorbarxy([0.6:11.6],new_meas.median,[0.6:11.6].*0,(new_meas.p75-new_meas.median),[0.6:11.6].*0,(new_meas.median-new_meas.p25),...
%     'Color',[0 0 0],'LineStyle','-','Marker','s','MarkerFaceColor',[0 0 0],'LineWidth',1,'MarkerSize',9);