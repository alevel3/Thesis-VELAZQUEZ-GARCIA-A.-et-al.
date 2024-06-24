%% %%%%%%%%%%%%%%%%%%%%%%% Percentiles avg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%%--- best estimations ---%%
%%% 5th perctile (In situ) and mean (BC GP) %%%
%%% R2 → 0.35 %%%%
%%--------------------------%%

load 1_mat_files\Hysplit.mat %%%% Hysplit

lat=interp1(1:73,lat(:,1:end),1:1./6:73);
lon=interp1(1:73,lon(:,1:end),1:1./6:73);
alt=interp1(1:73,alt(:,1:end),1:1./6:73);
PBL=interp1(1:73,PBL(:,1:end),1:1./6:73);
rain=interp1(1:73,rain(:,1:end),1:1./6:73);

load 1_mat_files\INTERPLAY_BC.mat %%%% INTERPLAY BC
%--- test the R2 before --%
% idx=find(nanmax(rain)>=1);
% BC_GP(:,idx)=NaN;
% size(find(isnan(sum(BC_GP))==0))%%14191
% fitlm(EBC_corr(:,6)',sum(BC_GP))%0.165
%--- → R2:0.165 %


load 1_mat_files\BC_AE33_ATOLL_flag.mat %%% AE33

sector={'AWB';'ENE';'IND';'RCO';'REF_TRF';'TNR_Ship';'TRO_noRES';'TRO_RES'};
 

%% Avg calculation %%

val_time=6;%6hrs,%9hrs & 12hrs

for j=1:max(size(val_time))
    
time_BT_avg=time_BT(1):val_time(j)./24:time_BT(end);

for i=1:max(size(time_BT_avg))
    
idx=find(time_AE33>=time_BT_avg(i)&time_AE33<time_BT_avg(i)+val_time(j)./24&idx_flag==0);
idx2=find(time_BT>=time_BT_avg(i)&time_BT<time_BT_avg(i)+val_time(j)./24&nanmax(rain)<1);



if min(size(idx))>0 && min(size(idx2))>0
    
EBC_avg(i,:)=prctile(EBC(idx,:),5);
BC_GP_avg(1:size(lat,1),i)=mean(BC_GP(1:size(lat,1),idx2),2);

%%%--- sectors ---%%%
for k=1:max(size(sector))
    eval(['BC_GP_',sector{k},'_avg','(1:size(lat,1),i)=mean(BC_GP','_',sector{k},'(1:size(lat,1),idx2),2);'])
end
%------------------%

 %--- dis point ---%
 for lol=1:size(lat,1)
     dist_point_avg(lol,i)=pos2dist(lat(1,1),lon(1,1),lat(lol,i),lon(lol,i),2);
 end
% %--------------------%


else
    
EBC_avg(i,:)=ones(7,1).*NaN;
BC_GP_avg(1:size(lat,1),i)=NaN;
count_filter(1:size(lat,1),i)=NaN;

%%%--- sectors ---%%%
for k=1:max(size(sector))
    eval(['BC_GP_',sector{k},'_avg','(1:size(lat,1),i)=NaN;'])
end
%------------------%

end


end

LM=fitlm(EBC_avg(:,6),sum(BC_GP_avg))
% R2(j)=LM.Rsquared.Adjusted;
%%12hrs → 0.358
%%9hrs → 0.322 / yo 0.335
%%6hrs  → 0.302 / yo 0.288
end

save 1_mat_files/INTERPLAY_baseline_HRS6.mat EBC_avg BC_GP_avg time_BT_avg...
    BC_GP_AWB_avg BC_GP_ENE_avg BC_GP_IND_avg BC_GP_RCO_avg BC_GP_REF_TRF_avg BC_GP_TNR_Ship_avg...
    BC_GP_TRO_noRES_avg BC_GP_TRO_RES_avg dist_point_avg

%'AWB';'ENE';'IND';'RCO';'REF_TRF';'TNR_Ship';'TRO_noRES';'TRO_RES'
%EBC_corr