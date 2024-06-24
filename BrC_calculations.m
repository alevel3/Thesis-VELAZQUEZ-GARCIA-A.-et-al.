function [Contrib_BrC] = BrC_calculations(BC1,BC2,BC3,BC4,BC5,BC6,BC7)


lambda = [370	470	525	590	660	880	940];
SG_Ch1 = 18.47 ;
SG_Ch2 = 14.54 ;
SG_Ch3 = 13.14 ;
SG_Ch4 = 11.58 ;
SG_Ch5 = 10.35;
SG_Ch6 = 7.77 ;
SG_Ch7 = 7.19;

BABS1=BC1%.*(SG_Ch1/1000);
BABS2=BC2%.*(SG_Ch2/1000);
BABS3=BC3%.*(SG_Ch3/1000);
BABS4=BC4%.*(SG_Ch4/1000);
BABS5=BC5%.*(SG_Ch5/1000);
BABS6=BC6%.*(SG_Ch6/1000);
BABS7=BC7%.*(SG_Ch7/1000);

x = log(lambda./1e9);
LOG_BABS = log([BABS1 BABS2 BABS3 BABS4 BABS5 BABS6 BABS7]);
for i=1:max(size(BABS1))
    ans = polyfit(x,LOG_BABS(i,:),1);
    alpha(i,1) = real(ans(1));
    clear ans
end


AAE_880_940 = -log10(BABS7./BABS6)./log10(lambda(7)./lambda(6));%%%BC
AAE_cste = 1 ;%%AAE for BC


Babs_BC_370_AAE880_940(:,1) = BABS6(:).*(lambda(1)./lambda(6)).^(-AAE_880_940);
Babs_BrC_370_AAE880_940 = BABS1 - Babs_BC_370_AAE880_940;

Babs_BC_370_AAE1(:,1) = BABS6(:).*(lambda(1)./lambda(6)).^(-AAE_cste);
Babs_BrC_370_AAE1 = BABS1 - Babs_BC_370_AAE1;%% using AAE=1


Babs_BrC_370_Zhang  = BABS6(:).*(lambda(6)./lambda(1)).^(-AAE_cste);
Babs_BrC_370_Zhang = BABS1 - Babs_BrC_370_Zhang;%% using AAE=1



Contrib_BrC.AAE880_940  = Babs_BrC_370_AAE880_940./BABS1;
Contrib_BrC.AAE_1  = Babs_BrC_370_AAE1./BABS1;%% contribution to the ABS
Contrib_BrC.AAE_Zhang  = Babs_BrC_370_Zhang./BABS1;%% contribution to the ABS

% BrC = Babs_BrC_370./(SG_Ch1./1000);
% BC = BC6;

end








