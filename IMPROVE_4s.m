function [Ext_PM1_calc,RC] = IMPROVE_4s(NH42SO4,NH4NO3,OM,BC)%,NO2,Rayleigh);

% NH42SO4= Data.(Season{i}).NH42SO4;
% NH4NO3=Data.(Season{i}).NH4NO3;
% NH4Cl=Data.(Season{i}).NH4Cl;
% OM=Data.(Season{i}).Org;
% BC=Data.(Season{i}).BC_ff;
% BrC=Data.(Season{i}).BC_wb;
% NO2=Data.(Season{i}).NO2;


MEE_NH42SO4 = 2.2; %m2/g
MEE_NH4NO3 = 2.4;%m2/g
%MEE_NH4Cl = 4.3;%m2/g
MEE_OM = 2.8;%m2/g
MEE_BC = 10;%m2/g--Black carbon
%MEE_BrC = 17.5;%m2/g---<brown carbon    no  bc      wb
MEE_NO2 = 0.33;%m2/g


Ext_NH42SO4 = MEE_NH42SO4.*NH42SO4;%*1.53;
Ext_NH4NO3 = MEE_NH4NO3.*NH4NO3;%*1.53;%slope decrease
Ext_OM = MEE_OM.*OM;
Ext_BC = MEE_BC.*BC;

%Ext_NO2 = MEE_NO2.*NO2;


Ext_PM1_calc = sum([Ext_NH42SO4 Ext_NH4NO3 Ext_OM Ext_BC],2);%Ext_NO2  Rayleigh

RC = [MEE_NH42SO4*nanmean(NH42SO4);...
    MEE_NH4NO3*nanmean(NH4NO3);...
    MEE_OM*nanmean(OM)
    MEE_BC*nanmean(BC)];

end