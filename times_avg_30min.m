%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTE CODIGO ES EL MAS ACTUALIZADO 01 JUNE 2022%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

date_Neph_PM1 = datenum(2017,7,28,0,0,0);

%% Download
%% ACSM 2017 - 2019 %% → Latest dataset VR

% load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\ACSM\ACSM_2017_2019');
% 
% Time_ACSM=ACSM_Time;

ACSM=load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\ACSM\ACSM_2017_2019');

idx_PM1 = min(find(ACSM.ACSM_Time>datenum(2017,7,28,23,59,0)));%(idx_PM1:end)%-->2019
size(find(isnan(ACSM.Org)==0))
Time_ACSM = ACSM.ACSM_Time(idx_PM1:end);
NH4 = ACSM.NH4(idx_PM1:end);
NO3 = ACSM.NO3(idx_PM1:end);
Org = ACSM.Org(idx_PM1:end);
SO4 = ACSM.SO4(idx_PM1:end);
Chl = ACSM.Chl(idx_PM1:end);
clear ACSM

%% AE33 2016 - 2021 %% → Latest dataset JB

load mat_files/BC_AE33_ATOLL_flag.mat

%EBC(idx_flag==1,:) = NaN;
idx=find(time_AE33>datenum(2019,1,1)&time_AE33<datenum(2019,2,1));
EBC(idx,:) = NaN;

%% Neph %% → Latest dataset SC
%%%% PM1 → 28-Jul-2017 13:28:00 %%%
Neph = load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Neph\Neph_EBAS_2017_2019_corrected');%%RH bellow 45%
lim_RH = 45; %f(RH) Pitchford = 1.53
Neph.Scat_G(Neph.RH>lim_RH)=NaN;
Neph.Scat_B(Neph.RH>lim_RH)=NaN;
Neph.Scat_R(Neph.RH>lim_RH)=NaN;
%% MTO

MTO=load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\MTO\MTO_2017_2019');

%% ATMO
ATMO = load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\ATMO\ATMO_2017_2019');

%% SMPS_EBAS ->Agregado el 16/02/2021
SMPS = load('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\SMPS\SMPS_2017-2019_EBAS');

%% Checking size
size(find(isnan(Neph.Scat_G)==0))%% 816 697; EBAS:895 757
size(find(isnan(EBC(:,6))==0))%% 887 511; EBAS:870507 / 2246190 since 2016
size(find(isnan(Org)==0))%% 26 549 / 107417 since 2016
size(find(isnan(MTO.Temp)==0))%% 1 199 690
size(find(isnan(ATMO.NO2)==0))%% 1 215 404
size(find(isnan(SMPS.SD)==0))%% 185 108
% size(find(isnan(RI_Real_440)==0))%% 58
%% Loops_time
for i=1:max(size(Time_ACSM))
    %Neph
    idx=find(Neph.Time_Neph>(Time_ACSM(i)-29./(24.*60))&Neph.Time_Neph<=Time_ACSM(i)&Neph.flag_Neph==0);
        
    if size(idx,1)>0 %threshold for data coverage before averaging
        %if size(idx,1)>0.75
        Scat_G(i,1)=nanmean(Neph.Scat_G(idx));
        Scat_B(i,1)=nanmean(Neph.Scat_B(idx)); 
        Scat_R(i,1)=nanmean(Neph.Scat_R(idx)); 
        BkScat_R(i,1) = nanmean(Neph.BkScat_B(idx));
        BkScat_G(i,1) = nanmean(Neph.BkScat_G(idx));
        BkScat_B(i,1) = nanmean(Neph.BkScat_R(idx));
        RH(i,1) = nanmean(Neph.RH(idx));
        Press(i,1) = nanmean(Neph.Press(idx));
    else
        Scat_G(i,1)=NaN;
        Scat_B(i,1)=NaN;
        Scat_R(i,1)=NaN;
        BkScat_R(i,1) = NaN;
        BkScat_G(i,1) = NaN;
        BkScat_B(i,1) = NaN;
        RH(i,1) = NaN;
        Press(i,1) = NaN;
    end
    
    %%% AE33    %%% before→idx=find(Time_Aeth>(Time_ACSM(i)-29./(24.*60))& Time_Aeth<=(Time_ACSM(i)));
    
    idx=find(time_AE33>(Time_ACSM(i)-29./(24.*60))&time_AE33<=Time_ACSM(i)&idx_flag==0);
    
    if size(idx,1)>0 %threshold for data coverage before averaging
        EBC_corr(i,:)=nanmean(EBC(idx,:));        
    else
        EBC_corr(i,:)=ones(7,1).*NaN;
    end
    
    %%% MTO
    idx=find(MTO.Time_MTO>(Time_ACSM(i)-29./(24.*60))&MTO.Time_MTO<=Time_ACSM(i));
    
    if size(idx,1)>0 %threshold for data coverage before averaging
        %if size(idx,1)>0.75
        Temperature(i,1)=nanmean(MTO.Temp(idx));
        Press_abs(i,1)=nanmean(MTO.P_abs(idx));
        Wind_sp(i,1)=nanmean(MTO.W_speed(idx));
        humidity(i,1)=nanmean(MTO.Humidity(idx));
        rain(i,1)=nanmean(MTO.Rain(idx));
        orientacion(i,1)=nanmean(MTO.Orientation(idx));
        Temp_int(i,1)=nanmean(MTO.Temp_INT(idx));
        
    else
        Temperature(i,1)=NaN;
        Press_abs(i,1)=NaN;
        Wind_sp(i,1)=NaN;
        humidity(i,1)=NaN;
        rain(i,1)=NaN;
        orientacion(i,1)=NaN;
        Temp_int(i,1)=NaN;
    end
    
    %%ATMO
    idx=find(ATMO.Time_ATMO>(Time_ACSM(i)-29./(24.*60))&ATMO.Time_ATMO<=Time_ACSM(i));
    
    if size(idx,1)>0 %threshold for data coverage before averaging
        %if size(idx,1)>0.75
        date_init(i,1) = nanmean(ATMO.date_init(idx));
        PM10(i,1) = nanmean(ATMO.PM10(idx));
        PM25(i,1) = nanmean(ATMO.PM25(idx));
        NO(i,1) = nanmean(ATMO.NO(idx));
        N2(i,1) = nanmean(ATMO.N2(idx));
        SO(i,1) = nanmean(ATMO.SO(idx));
        NO2(i,1) = nanmean(ATMO.NO2(idx));
    else
        date_init(i,1) = NaN;
        PM10(i,1) = NaN;
        PM25(i,1) = NaN;
        NO(i,1) = NaN;
        N2(i,1) = NaN;
        SO(i,1) = NaN;
        NO2(i,1) = NaN;
     end
%     %%SMPS_EBAS 12/02/2021
%     
%     %idx=find(Time_SMPS>(Time_ACSM(i)-29./(24.*60))& Time_SMPS<=(Time_ACSM(i)));
%     
%     idx=find(SMPS.Time_SMPS>(Time_ACSM(i)-29./(24.*60))&SMPS.Time_SMPS<=Time_ACSM(i));
%     
%     if size(idx,1)>0 %threshold for data coverage before averaging
%         %if size(idx,1)>0.75
%         %Dia(i,:,1) = nanmean(Dia_A(idx));%%MOROE THAN 1 COLUMN
%         Eff_radius_real(i,1) = nanmean(SMPS.Eff_radius_real(idx));
%         Eff_radius(i,1) = nanmean(SMPS.Eff_radius(idx));
%         PM_SMPS(i,1) = nanmean(SMPS.PM_SMPS(idx));
%         PV_SMPS(i,1) = nanmean(SMPS.PV_SMPS(idx));
%         dS(i,:,1) = nanmean(SMPS.dS(idx,:,1));%%MOROE THAN 1 COLUMN
%         dV(i,:,1) = nanmean(SMPS.dV(idx,:,1));%%MOROE THAN 1 COLUMN
%         P_SMPS(i,1) = nanmean(SMPS.P_SMPS(idx));
%         T_SMPS(i,1) = nanmean(SMPS.T_SMPS(idx));
%         RH_SMPS(i,1) = nanmean(SMPS.RH_SMPS(idx));
%         SD(i,:,1) = nanmean(SMPS.SD(idx,:,1));%%MOROE THAN 1 COLUMN
%         
%     else
%         %Dia(i,1) = NaN;
%         Eff_radius_real(i,1) = NaN;
%         Eff_radius(i,1) = NaN;
%         PM_SMPS(i,1) = NaN;
%         PV_SMPS(i,1) = NaN;
%         dS(i,1) = NaN;
%         dV(i,1) = NaN;
%         P_SMPS(i,1) = NaN;
%         T_SMPS(i,1) = NaN;
%         RH_SMPS(i,1) = NaN;
%         SD(i,1) = NaN;
%         
%     end
 
    
    
end

AE33_Abs_calc

[BC_ff,BC_wb,alpha] = sandradewi(BC1,BC2,BC3,BC4,BC5,BC6,BC7);

idx3=find(BC_ff<0|BC_wb<0);%% size:30 874 
BC_ff(idx3)=NaN;
BC_wb(idx3)=NaN;

%% Using the latest version of AE33, Neph and SMPS
save mat_files/In_situ_obs(2).mat 'BC_ff' 'BC_wb' 'BC1' 'BC2' 'BC3' 'BC4' 'BC5' 'BC6' 'BC7' 'Abs_BC1' 'Abs_BC2' 'Abs_BC3' 'Abs_BC4' 'Abs_BC5' 'Abs_BC6' 'Abs_BC7' ...
    'Scat_G' 'Scat_B' 'Scat_R' 'BkScat_R' 'BkScat_G' 'BkScat_B' 'RH' 'Press' ...
    'Time_ACSM' 'Org' 'NH4' 'NO3' 'Chl' 'SO4' ...
    'Temperature' 'Press_abs' 'Wind_sp' 'rain' 'orientacion' 'humidity' 'Temp_int' ...
    'date_init' 'PM10' 'PM25' 'NO' 'N2' 'SO' 'NO2'% ...
    %'Dia' 'Eff_radius_real' 'Eff_radius' 'PM_SMPS' 'PV_SMPS' 'P_SMPS' 'T_SMPS' 'RH_SMPS' 'dS' 'dV' 'SD'

% save('C:\Users\alejandra.garcia\Documents\MATLAB\LOA\My matrix\Matrix_all_period 2017-2020\Matrix_avareged_29min\LOA_EBASdata_avg_2017_2019_insitu_14MARCH2022',...
%     'BC_ff','BC_wb','BC1','BC2','BC3','BC4','BC5','BC6','BC7','Abs_BC1','Abs_BC2','Abs_BC3','Abs_BC4','Abs_BC5','Abs_BC6','Abs_BC7','idx_flag_AE33','idx_flag_JB_PM1',...
%     'Scat_G','Scat_B','Scat_R','BkScat_R','BkScat_G','BkScat_B','RH','Press',...
%     'Time_ACSM','Org','NH4','NO3','Chl','SO4',...
%     'Temperature','Press_abs','Wind_sp','rain','orientacion','humidity','Temp_int',...
%     'date_init','PM10','PM25','NO','N2','SO','NO2',...
%     'Dia','Eff_radius_real','Eff_radius','PM_SMPS','PV_SMPS','P_SMPS','T_SMPS','RH_SMPS','dS','dV','SD');
