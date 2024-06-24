%% AAE in time %%
for i = 1:max(size(Time_ACSM))
    
    x(i,:)=log10(Abs_BC1(i)./Abs_BC2(i));
    y(i,:)=log10(370/470);
    AAE_370_470(i,:)=-(x(i)/y(i));
    
    xx(i,:)=log10(Abs_BC4(i)./Abs_BC4(i));
    yy(i,:)=log10(590/660);
    AAE_590_660(i,:)=-(xx(i)/yy(i));
        
end

lambda1 = 370;% La longitud de onda del primer punto (Prop. Opt. Medida) 470;%AE16
lambda2 = 470;% La longitud de onda del punto final (Prop. Opt. Medida)880;%AE33
lambda_x = 450;%% La onda en la que quieres tener la propiedad optica (ej.lambda Neph)
x=log10(Abs_BC1./Abs_BC2);
y=log10(lambda1/lambda2);
AAE=-(x/y);
[Abs_450] = change_wavelength(Abs_BC2,AAE,lambda2,lambda_x);
clearvars lambda1 lambda2 lambda_x x y AAE 

EXT_450 = Scat_B+abs(Abs_450);

EXT_525 = Scat_G+Abs_BC3;

lambda4 = 590;% La longitud de onda del primer punto (Prop. Opt. Medida) 470;%AE16
lambda5 = 660;% La longitud de onda del punto final (Prop. Opt. Medida)880;%AE33
lambda_x = 635;%% La onda en la que quieres tener la propiedad optica (ej.lambda Neph)
x=log10(Abs_BC4./Abs_BC5);
y=log10(lambda4/lambda5);
AAE=-(x/y);
[Abs_635] = change_wavelength(Abs_BC4,AAE,lambda5,lambda_x);
clearvars lambda4 lambda5 lambda_x x y AAE 

EXT_635 = Scat_R+abs(Abs_635);

    Abs_BC(i,:) = Abs_880(i).*((lambda_x/lambda_880).^(-AAE));
    AbsBrC(i,:) = Abs_x(i) - Abs_BC(i,:);


