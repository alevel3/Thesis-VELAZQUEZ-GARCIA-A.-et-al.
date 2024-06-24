
%load mat_files/BC_AE33_ATOLL_flag.mat


SG = [18.47 ;14.54 ;13.14 ; 11.58 ;10.35;7.77 ; 7.19];%%
%%%%%% ACTRIS harmonization factor %%%%%%%%%%%%%%%%%%%
h = 2.21;%% For filters M8060 (ATOLL: 2016 - 2017)
H = 1.73;%% For filters M8060 (ATOLL: 2018 - ongoing)

BC1=EBC_corr(:,1);
BC2=EBC_corr(:,2);
BC3=EBC_corr(:,3);
BC4=EBC_corr(:,4);
BC5=EBC_corr(:,5);
BC6=EBC_corr(:,6);
BC7=EBC_corr(:,7);

%%%%%%% 2016-2017 %%%%%%
idx=find(Time_ACSM>datenum(2016,12,01) & Time_ACSM<=datenum(2018,01,01));
BC1_2016_2017 = (BC1(idx).*SG(1))./h;
BC2_2016_2017 = (BC2(idx).*SG(2))./h;
BC3_2016_2017 = (BC3(idx).*SG(3))./h;
BC4_2016_2017 = (BC4(idx).*SG(4))./h;
BC5_2016_2017 = (BC5(idx).*SG(5))./h;
BC6_2016_2017 = (BC6(idx).*SG(6))./h;
BC7_2016_2017 = (BC7(idx).*SG(7))./h;

%%%%%%% 2018 - 2021 %%%%%%
idx1=find(Time_ACSM>datenum(2018,01,01) & Time_ACSM<=datenum(2022,01,01));
BC1_2018_2021 = (BC1(idx1).*SG(1))./H;
BC2_2018_2021 = (BC2(idx1).*SG(2))./H;
BC3_2018_2021 = (BC3(idx1).*SG(3))./H;
BC4_2018_2021 = (BC4(idx1).*SG(4))./H;
BC5_2018_2021 = (BC5(idx1).*SG(5))./H;
BC6_2018_2021 = (BC6(idx1).*SG(6))./H;
BC7_2018_2021 = (BC7(idx1).*SG(7))./H;

Abs_370nm = [];
Abs_470nm = [];
Abs_520nm = [];
Abs_590nm = [];
Abs_660nm = [];
Abs_880nm = [];
Abs_940nm = [];
%%%%%%% 2016-2017 %%%%%%
Abs_370nm = [Abs_370nm;BC1_2016_2017];
Abs_470nm = [Abs_470nm;BC2_2016_2017];
Abs_520nm = [Abs_520nm;BC3_2016_2017];
Abs_590nm = [Abs_590nm;BC4_2016_2017];
Abs_660nm = [Abs_660nm;BC5_2016_2017];
Abs_880nm = [Abs_880nm;BC6_2016_2017];
Abs_940nm = [Abs_940nm;BC7_2016_2017];
%%%%%%% 2018 - 2021 %%%%%%
Abs_370nm = [Abs_370nm;BC1_2018_2021];
Abs_470nm = [Abs_470nm;BC2_2018_2021];
Abs_520nm = [Abs_520nm;BC3_2018_2021];
Abs_590nm = [Abs_590nm;BC4_2018_2021];
Abs_660nm = [Abs_660nm;BC5_2018_2021];
Abs_880nm = [Abs_880nm;BC6_2018_2021];
Abs_940nm = [Abs_940nm;BC7_2018_2021];

Abs_BC1=Abs_370nm;
Abs_BC2=Abs_470nm;
Abs_BC3=Abs_520nm;
Abs_BC4=Abs_590nm;
Abs_BC5=Abs_660nm;
Abs_BC6=Abs_880nm;
Abs_BC7=Abs_940nm;

clearvars SG H h idx idx1 BC1_2016_2017 BC2_2016_2017 BC3_2016_2017 BC4_2016_2017 BC5_2016_2017 BC6_2016_2017 BC7_2016_2017 BC1_2018_2021 BC2_2018_2021 BC3_2018_2021 BC4_2018_2021 BC5_2018_2021 BC6_2018_2021 BC7_2018_2021 Abs_370nm Abs_470nm Abs_520nm Abs_590nm Abs_660nm Abs_880nm Abs_940nm
%BC1 BC2 BC3 BC4 BC5 BC6 BC7