function [BC_ff,BC_wb,alpha] = sandradewi(BC1,BC2,BC3,BC4,BC5,BC6,BC7)

lambda = [370	470	525	590	660	880	940];
SG = [18.47 ;14.54 ;13.14 ; 11.58 ;10.35;7.77 ; 7.19];
% % BABS1=BC1;%.*SG(1)./1000;%% already converted in my data (Ale)!!!!Attenuation to mass
% % BABS2=BC2;%.*SG(2)./1000;
% % BABS3=BC3;%.*SG(3)./1000;
% % BABS4=BC4;%.*SG(4)./1000;
% % BABS5=BC5;%.*SG(5)./1000;
% % BABS6=BC6;%.*SG(6)./1000;
% % BABS7=BC7;%.*SG(7)./1000;
SG_Ch1 = 18.47 ;
SG_Ch2 = 14.54 ;
SG_Ch3 = 13.14 ;
SG_Ch4 = 11.58 ;
SG_Ch5 = 10.35;
SG_Ch6 = 7.77 ;
SG_Ch7 = 7.19;
BABS1=BC1*(SG_Ch1/1000);%%absportion
BABS2=BC2*(SG_Ch2/1000);
BABS3=BC3*(SG_Ch3/1000);
BABS4=BC4*(SG_Ch4/1000);
BABS5=BC5*(SG_Ch5/1000);
BABS6=BC6*(SG_Ch6/1000);
BABS7=BC7*(SG_Ch7/1000);

x = log(lambda./1e9);%%log10
LOG_BABS = log([BABS1 BABS2 BABS3 BABS4 BABS5 BABS6 BABS7]);



for i=1:max(size(BABS1))
    ans = polyfit(x,LOG_BABS(i,:),1);
    alpha(i,1) = real(ans(1));
    clear ans
end

Alpha_ff = 1;
Alpha_wb = 2;




BABS_WB_470nm = (1./(1-(lambda(2)/lambda(6))^Alpha_wb*(lambda(2)/lambda(6))^(-Alpha_ff)))*(BABS2 - (lambda(2)/lambda(6))^(-Alpha_ff)*BABS6);
BABS_WB_880nm = BABS_WB_470nm*(lambda(2)/lambda(6))^Alpha_wb;
BABS_FF_470nm = BABS2 - BABS_WB_470nm;
BABS_FF_880nm = BABS6 - BABS_WB_880nm;



% %%attenuation
% BC66 = BC6*1000./SG(6);%%BC6*(SG_Ch6/1000);

BC_ff = BC6.*  BABS_FF_880nm./BABS6;%%BC6->attenuation*abs/conc
BC_wb = BC6 - BC_ff;

end












