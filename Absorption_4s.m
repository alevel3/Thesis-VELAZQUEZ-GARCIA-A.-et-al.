%%% FUNCTION TO CALCULATE MAE OF 4 species: NH42SO4,NHNO3,Org,BC %%%%%
%%% Rayilegt => scattering of gases
%%fitlm => 
function[MAEs,Abs_370_calc,Abs_470_calc,Abs_520_calc,Abs_590_calc,Abs_660_calc,Abs_880_calc,Abs_950_calc,RC,...
    res1,res2,res3,res4,res5,res6,res7,hi] = Absorption_4s(NH42SO4,NH4NO3,Org,BC,Abs_BC1,Abs_BC2,Abs_BC3,Abs_BC4,Abs_BC5,Abs_BC6,Abs_BC7)%,NO2)
%% Xdata -> Input
%Org = ACSM.Org;
xdata=[NH42SO4 NH4NO3 Org BC];
%% Absorption in PM1 (Bap -> ydata)
ydata = [Abs_BC1 Abs_BC2 Abs_BC3 Abs_BC4 Abs_BC5 Abs_BC6 Abs_BC7];
%% Removing NaN data points in xdata and ydata
%     isnan(xdata);
idx=max(isnan(xdata),[],2);
xdata(idx,:)=[];
ydata(idx,:)=[];

%     isnan(ydata);
%     max(isnan(ydata),[],2);
idx=max(isnan(ydata),[],2);
ydata(idx,:)=[];
xdata(idx,:)=[];

% % Leverage is a measure of the effect of a particular observation on 
% % the regression predictions due to the position of that observation in 
% % the space of the inputs. In general, the farther a point is from the 
% % center of the input space, the more leverage it has. Because the sum of 
% % the leverage values is p, an observation i can be considered as an outlier 
% % if its leverage substantially exceeds the mean leverage value, p/n, 
% % for example, a value larger than 2*p/n.
sum(leverage(xdata))
hi=leverage(xdata);
%hi_AOPs=leverage(ydata);

%% %decay model:y=b(1)*xdata(:,1)+b(2)*xdata(:,2)+b(3)*xdata(:,3)+b(4)*xdata(:,4)+b(5)*xdata(:,5)+b(6)*xdata(:,6)
%% %Parameters bounded%
%ub=[0.397,0.035,x,0.112,4.604,0.246];%% Maximum values of each specie: NH4SO4, NH4NO3, NH4Cl, Org, Bc, BrC%

%%%--- Intercept ---%%
% X = [ones(size(xdata(:,1))) xdata(:,1) xdata(:,2) xdata(:,3) xdata(:,4)];
% lb=[0,0,0,0,0];
% ub=[];
% x0=[0,0,0,0,0];

%%%--- No intercepti ---%%
lb=[0,0,0,0];
ub=[];
x0=[0,0,0,0];
%%------------------%%%
%% %%%CReate the model
%Absportion -> MAE %%%replace b by  'MSE'|'MEE'|'MAE'
%fun = @(MAE,xdata)MAE(1)*xdata(:,1)+MAE(2)*xdata(:,2)+MAE(3)*xdata(:,3)+MAE(4)*xdata(:,4)+MAE(5)*xdata(:,5)+MAE(6)*xdata(:,6);
%No Scattering
% fun = @(b,xdata)b(1)*xdata(:,1)+b(2)*xdata(:,2)+b(3)*xdata(:,3)+b(4)*xdata(:,4)+b(5)*xdata(:,5)+b(6)*xdata(:,6)
%Extinction -> MEE %%%replace b by  'MSE'|'MEE'|'MAE'
%fun = @(MEE,xdata)MEE(1)*xdata(:,1)+MEE(2)*xdata(:,2)+MEE(3)*xdata(:,3)+MEE(4)*xdata(:,4)+MEE(5)*xdata(:,5)+MEE(6)*xdata(:,6);
%% %%%Initian guess
%x0=[0,0,0,0];%%%No particles% No exctinction
% x0=[0,0,0,0,0,0,0]%% b(0)

%% %% Solve the bounded fitting problem %WITH Bap
% %%%Calculating MAE%%
%period_5min stand for Sca. 450, 525, 635
% fun=@(coe,X)(coe(1)*X(:,1)+coe(2)*X(:,2)+coe(3)*X(:,3)+coe(4)*X(:,4)+coe(5)*X(:,5)+coe(6)*X(:,6));
% MAE_470=lsqcurvefit(fun,x0,xdata,ydata(:,1),lb,ub);
% MAE_525=lsqcurvefit(fun,x0,xdata,ydata(:,2),lb,ub);
% MAE_880=lsqcurvefit(fun,x0,xdata,ydata(:,3),lb,ub);
% MAE_370=lsqcurvefit(fun,x0,xdata,ydata(:,4),lb,ub);
% B = regress(ydata(:,2),xdata);
% scatter(MAE_525,B,20,1:6,'filled')
% colorbar
%% Model using "fitlm"
%mdl = fitlm(xdata,ydata,'RobustOpts','on')
%dlm = fitlm(X,y,'Intercept',false);
%dlm = fitlm(X,y,'y~x1+x2+x3+x4+x5+x6');
%% Abs_370
%%FITLM function
%mdl_370 = fitlm(xdata,ydata(:,1),'RobustOpts','on')
% mdl_370 = fitlm(xdata,ydata(:,1),'Intercept',false)
% %mdl_370 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_370 = fitlm(xdata,ydata(:,1))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC1 = table2struct(mdl_370.Coefficients(:,1),'ToScalar',true);
% MAE_370 = MAE_BC1.Estimate
% CI_370 = coefCI(mdl_370)%%(,aplha)
%%lsqlin

%%%--- Intercept ---%%%
% [x,resnorm,res1,exitflag,output,lambda] = lsqlin(X,ydata(:,1),[],[],[],[],lb,ub,x0);
% MAE_370 = x
% yfit = X*MAE_370;
% resid = ydata(:,1) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_370 = ci';
% P=polyfitn(ydata(:,1),X*x,1);
% P.R2
% P.Coefficients
%%%%%----------------%%%%%%%%

%%%--- No intercept ---%%%
[x,resnorm,res1,exitflag,output,lambda] = lsqlin(xdata,ydata(:,1),[],[],[],[],lb,ub,x0);
MAE_370 = x
yfit = xdata*MAE_370;
resid = ydata(:,1) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_370 = ci'
P=polyfitn(ydata(:,1),xdata*x,1)
P.R2
P.Coefficients
%%%%%----------------%%%%%%%%
%% Abs_470
%%fitlm
% %mdl_470 = fitlm(xdata,ydata(:,2),'RobustOpts','on')
% mdl_470 = fitlm(xdata,ydata(:,2),'Intercept',false)
% %mdl_470 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_470 = fitlm(xdata,ydata(:,2))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC2 = table2struct(mdl_470.Coefficients(:,1),'ToScalar',true);
% MAE_470 = MAE_BC2.Estimate
% CI_470 = coefCI(mdl_470)%%(,aplha)
%%lsqlin
%%%%-------Intercept ------%%%%%
% [x,resnorm,res2,exitflag,output,lambda] = lsqlin(X,ydata(:,2),[],[],[],[],lb,ub,x0);
% MAE_470 = x
% yfit = X*MAE_470;
% resid = ydata(:,2) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_470 = ci';
% P=polyfitn(ydata(:,2),X*x,1);
% P.R2
% P.Coefficients

%%%%------No intercept --------%%%
[x,resnorm,res2,exitflag,output,lambda] = lsqlin(xdata,ydata(:,2),[],[],[],[],lb,ub,x0);
MAE_470 = x
yfit = xdata*MAE_470;
resid = ydata(:,2) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_470 = ci'
P=polyfitn(ydata(:,2),xdata*x,1)
P.R2
P.Coefficients
%% Abs_520
%%fitlm
% %mdl_520 = fitlm(xdata,ydata(:,3),'RobustOpts','on')
% mdl_520 = fitlm(xdata,ydata(:,3),'Intercept',false)
% %mdl_520 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_520 = fitlm(xdata,ydata(:,3))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC3 = table2struct(mdl_520.Coefficients(:,1),'ToScalar',true);
% MAE_520 = MAE_BC3.Estimate
% CI_520 = coefCI(mdl_520)%%(,aplha)
%%lsqlin

%%%%------Intercept-------%%%%
% [x,resnorm,res3,exitflag,output,lambda] = lsqlin(X,ydata(:,3),[],[],[],[],lb,ub,x0);
% MAE_520 = x
% yfit = X*MAE_520;
% resid = ydata(:,3) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_520 = ci';
% P=polyfitn(ydata(:,3),X*x,1);
% P.R2
% P.Coefficients


%%%%------- No intercept ------%%%%
[x,resnorm,res3,exitflag,output,lambda] = lsqlin(xdata,ydata(:,3),[],[],[],[],lb,ub,x0);
MAE_520 = x
yfit = xdata*MAE_520;
resid = ydata(:,3) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_520 = ci'
P=polyfitn(ydata(:,3),xdata*x,1)
P.R2
P.Coefficients
%% Abs_590
%%fitlm
% %mdl_590 = fitlm(xdata,ydata(:,4),'RobustOpts','on')
% mdl_590 = fitlm(xdata,ydata(:,4),'Intercept',false)
% %mdl_590 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_590 = fitlm(xdata,ydata(:,4))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC4 = table2struct(mdl_590.Coefficients(:,1),'ToScalar',true);
% MAE_590 = MAE_BC4.Estimate
% CI_590 = coefCI(mdl_590)%%(,aplha)
%%lsqlin

%%%%-----Intercept-------%%%
% [x,resnorm,res4,exitflag,output,lambda] = lsqlin(X,ydata(:,4),[],[],[],[],lb,ub,x0);
% MAE_590 = x
% yfit = X*MAE_590;
% resid = ydata(:,4) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_590 = ci';
% P=polyfitn(ydata(:,4),X*x,1);
% P.R2
% P.Coefficients


%%%%------No intercept------%%%%
[x,resnorm,res4,exitflag,output,lambda] = lsqlin(xdata,ydata(:,4),[],[],[],[],lb,ub,x0);
MAE_590 = x
yfit = xdata*MAE_590;
resid = ydata(:,4) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_590 = ci'
P=polyfitn(ydata(:,4),xdata*x,1)
P.R2
P.Coefficients
%% Abs_660
%%%fitlm
% %mdl_660 = fitlm(xdata,ydata(:,5),'RobustOpts','on')
% mdl_660 = fitlm(xdata,ydata(:,5),'Intercept',false)
% %mdl_660 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_660 = fitlm(xdata,ydata(:,5))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC5 = table2struct(mdl_660.Coefficients(:,1),'ToScalar',true);
% MAE_660 = MAE_BC5.Estimate
% CI_660 = coefCI(mdl_660)%%(,aplha)
%%lsqlin

%%%----Intercerpt---%%%
% [x,resnorm,res5,exitflag,output,lambda] = lsqlin(X,ydata(:,5),[],[],[],[],lb,ub,x0);
% MAE_660 = x
% yfit = X*MAE_660;
% resid = ydata(:,5) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_660 = ci';
% P=polyfitn(ydata(:,5),X*x,1);
% P.R2
% P.Coefficients

%%%-----NO intercep----%%%
[x,resnorm,res5,exitflag,output,lambda] = lsqlin(xdata,ydata(:,5),[],[],[],[],lb,ub,x0);
MAE_660 = x
yfit = xdata*MAE_660;
resid = ydata(:,5) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_660 = ci'
P=polyfitn(ydata(:,5),xdata*x,1)
P.R2
P.Coefficients
%% ABs_880
%%fitlm
% %mdl_880 = fitlm(xdata,ydata(:,6),'RobustOpts','on')
% mdl_880 = fitlm(xdata,ydata(:,6),'Intercept',false)
% %mdl_880 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_880 = fitlm(xdata,ydata(:,6))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC6 = table2struct(mdl_880.Coefficients(:,1),'ToScalar',true);
% MAE_880 = MAE_BC6.Estimate
% CI_880 = coefCI(mdl_880)%%(,aplha)
%%%lsqlin
%%----Interept-----%%
% [x,resnorm,res6,exitflag,output,lambda] = lsqlin(X,ydata(:,6),[],[],[],[],lb,ub,x0);
% MAE_880 = x
% yfit = X*MAE_880;
% resid = ydata(:,6) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_880 = ci';
% P=polyfitn(ydata(:,6),X*x,1);
% P.R2
% P.Coefficients

%%%%--- No intercept ----%%%
[x,resnorm,res6,exitflag,output,lambda] = lsqlin(xdata,ydata(:,6),[],[],[],[],lb,ub,x0);
MAE_880 = x
yfit = xdata*MAE_880;
resid = ydata(:,6) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_880 = ci'
P=polyfitn(ydata(:,6),xdata*x,1)
P.R2
P.Coefficients
%% Abs_940
% %mdl_940 = fitlm(xdata,ydata(:,7),'RobustOpts','on')
% mdl_940 = fitlm(xdata,ydata(:,7),'Intercept',false)
% %mdl_940 = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
% %mdl_940 = fitlm(xdata,ydata(:,7))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MAE_BC7 = table2struct(mdl_940.Coefficients(:,1),'ToScalar',true);
% MAE_940 = MAE_BC7.Estimate
% CI_940 = coefCI(mdl_940)%%(,aplha)
%%%lsqlin
%%%----Interept----%%
% [x,resnorm,res7,exitflag,output,lambda] = lsqlin(X,ydata(:,7),[],[],[],[],lb,ub,x0);
% MAE_950 = x
% yfit = X*MAE_950;
% resid = ydata(:,7) - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_950 = ci';
% P=polyfitn(ydata(:,7),X*x,1);
% P.R2
% P.Coefficients

%%%----No interept----%%%
[x,resnorm,res7,exitflag,output,lambda] = lsqlin(xdata,ydata(:,7),[],[],[],[],lb,ub,x0);
MAE_950 = x
yfit = xdata*MAE_950;
resid = ydata(:,7) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_950 = ci'
P=polyfitn(ydata(:,7),xdata*x,1)
P.R2
P.Coefficients
%% %SAVING FILES%%%
% MAE_470=[MAE_470]';
% MAE_525=[MAE_525]';
% MAE_880=[MAE_880]';
% MAE_370=[MAE_370]';
% Species={'NH42SO2' 'NH42NO3' 'NH42Chl' 'Org' 'BC' 'BrC'}
% Table_MEE=table(Species,MEE_450,MEE_525,MEE_635)

% MAEs.W370 = MAE_370;
% MAEs.W470 = MAE_470;
% MAEs.W520 = MAE_520;
% MAEs.W590 = MAE_590;
% MAEs.W660 = MAE_660;
% MAEs.W880 = MAE_880;
% MAEs.W940 = MAE_940;

% MAEs.W370 = round(MAE_370,2);
% MAEs.W470 = round(MAE_470,2);
% MAEs.W520 = round(MAE_520,2);
% MAEs.W590 = round(MAE_590,2);
% MAEs.W660 = round(MAE_660,2);
% MAEs.W880 = round(MAE_880,2);
% MAEs.W950 = round(MAE_950,2);

MAEs.W370 = round(MAE_370,1);
MAEs.W470 = round(MAE_470,1);
MAEs.W520 = round(MAE_520,1);
MAEs.W590 = round(MAE_590,1);
MAEs.W660 = round(MAE_660,1);
MAEs.W880 = round(MAE_880,1);
MAEs.W950 = round(MAE_950,1);

%% Retriving Abs
%% 370
Abs_370_NH42SO4 = MAE_370(1).*NH42SO4;
Abs_370_NH4NO3 = MAE_370(2).*NH4NO3;
Abs_370_Org = MAE_370(3).*Org;
Abs_370_BC = MAE_370(4).*BC;

%Abs_370_NO2 = 0.33.*NO2;

Abs_370_calc = sum([Abs_370_NH42SO4 Abs_370_NH4NO3 Abs_370_Org Abs_370_BC],2);%Abs_370_NO2

% RC.W370 = [MAE_370(1)*nanmean(NH42SO4);...
%     MAE_370(2)*nanmean(NH4NO3);...
%     MAE_370(3)*nanmean(Org);...
%     MAE_370(4)*nanmean(BC)];

RC.W370 = [MAEs.W370(1)*nanmean(NH42SO4);...
    MAEs.W370(2)*nanmean(NH4NO3);...
    MAEs.W370(3)*nanmean(Org);...
    MAEs.W370(4)*nanmean(BC)];

%% 470
Abs_470_NH42SO4 = MAE_470(1).*NH42SO4;
Abs_470_NH4NO3 = MAE_470(2).*NH4NO3;
Abs_470_Org = MAE_470(3).*Org;
Abs_470_BC = MAE_470(4).*BC;

%Abs_470_NO2 = 0.33.*NO2;

Abs_470_calc = sum([Abs_470_NH42SO4 Abs_470_NH4NO3 Abs_470_Org Abs_470_BC],2);%Abs_470_NO2

% RC.W470 = [MAE_470(1)*nanmean(NH42SO4);...
%     MAE_470(2)*nanmean(NH4NO3);...
%     MAE_470(3)*nanmean(Org);...
%     MAE_470(4)*nanmean(BC)];

RC.W470 = [MAEs.W470(1)*nanmean(NH42SO4);...
    MAEs.W470(2)*nanmean(NH4NO3);...
    MAEs.W470(3)*nanmean(Org);...
    MAEs.W470(4)*nanmean(BC)];
%% 520
Abs_520_NH42SO4 = MAE_520(1).*NH42SO4;
Abs_520_NH4NO3 = MAE_520(2).*NH4NO3;
Abs_520_Org = MAE_520(3).*Org;
Abs_520_BC = MAE_520(4).*BC;

%Abs_520_NO2 = 0.33.*NO2;

Abs_520_calc = sum([Abs_520_NH42SO4 Abs_520_NH4NO3 Abs_520_Org Abs_520_BC],2);%Abs_520_NO2

% RC.W520 = [MAE_520(1)*nanmean(NH42SO4);...
%     MAE_520(2)*nanmean(NH4NO3);...
%     MAE_520(3)*nanmean(Org);...
%     MAE_520(4)*nanmean(BC)];

RC.W520 = [MAEs.W520(1)*nanmean(NH42SO4);...
    MAEs.W520(2)*nanmean(NH4NO3);...
    MAEs.W520(3)*nanmean(Org);...
    MAEs.W520(4)*nanmean(BC)];
%% 590
Abs_590_NH42SO4 = MAE_590(1).*NH42SO4;
Abs_590_NH4NO3 = MAE_590(2).*NH4NO3;
Abs_590_Org = MAE_590(3).*Org;
Abs_590_BC = MAE_590(4).*BC;

%Abs_590_NO2 = 0.33.*NO2;

Abs_590_calc = sum([Abs_590_NH42SO4 Abs_590_NH4NO3 Abs_590_Org Abs_590_BC],2);%Abs_590_NO2

% RC.W590 = [MAE_590(1)*nanmean(NH42SO4);...
%     MAE_590(2)*nanmean(NH4NO3);...
%     MAE_590(3)*nanmean(Org);...
%     MAE_590(4)*nanmean(BC)];

RC.W590 = [MAEs.W590(1)*nanmean(NH42SO4);...
    MAEs.W590(2)*nanmean(NH4NO3);...
    MAEs.W590(3)*nanmean(Org);...
    MAEs.W590(4)*nanmean(BC)];

%% 660
Abs_660_NH42SO4 = MAE_660(1).*NH42SO4;
Abs_660_NH4NO3 = MAE_660(2).*NH4NO3;
Abs_660_Org = MAE_660(3).*Org;
Abs_660_BC = MAE_660(4).*BC;

%Abs_660_NO2 = 0.33.*NO2;

Abs_660_calc = sum([Abs_660_NH42SO4 Abs_660_NH4NO3 Abs_660_Org Abs_660_BC],2);%Abs_660_NO2

% RC.W660 = [MAE_660(1)*nanmean(NH42SO4);...
%     MAE_660(2)*nanmean(NH4NO3);...
%     MAE_660(3)*nanmean(Org);...
%     MAE_660(4)*nanmean(BC)];

RC.W660 = [MAEs.W660(1)*nanmean(NH42SO4);...
    MAEs.W660(2)*nanmean(NH4NO3);...
    MAEs.W660(3)*nanmean(Org);...
    MAEs.W660(4)*nanmean(BC)];

%% 880
Abs_880_NH42SO4 = MAE_880(1).*NH42SO4;
Abs_880_NH4NO3 = MAE_880(2).*NH4NO3;
Abs_880_Org = MAE_880(3).*Org;
Abs_880_BC = MAE_880(4).*BC;

%Abs_880_NO2 = 0.33.*NO2;

Abs_880_calc = sum([Abs_880_NH42SO4 Abs_880_NH4NO3 Abs_880_Org Abs_880_BC],2);%Abs_880_NO2

% RC.W880 = [MAE_880(1)*nanmean(NH42SO4);...
%     MAE_880(2)*nanmean(NH4NO3);...
%     MAE_880(3)*nanmean(Org);...
%     MAE_880(4)*nanmean(BC)];

RC.W880 = [MAEs.W880(1)*nanmean(NH42SO4);...
    MAEs.W880(2)*nanmean(NH4NO3);...
    MAEs.W880(3)*nanmean(Org);...
    MAEs.W880(4)*nanmean(BC)];

%% 950
Abs_950_NH42SO4 = MAE_950(1).*NH42SO4;
Abs_950_NH4NO3 = MAE_950(2).*NH4NO3;
Abs_950_Org = MAE_950(3).*Org;
Abs_950_BC = MAE_950(4).*BC;

%Abs_940_NO2 = 0.33.*NO2;

Abs_950_calc = sum([Abs_950_NH42SO4 Abs_950_NH4NO3 Abs_950_Org Abs_950_BC],2);%Abs_940_NO2

% RC.W940 = [MAE_940(1)*nanmean(NH42SO4);...
%     MAE_940(2)*nanmean(NH4NO3);...
%     MAE_940(3)*nanmean(Org);...
%     MAE_940(4)*nanmean(BC)];

RC.W940 = [MAEs.W950(1)*nanmean(NH42SO4);...
    MAEs.W950(2)*nanmean(NH4NO3);...
    MAEs.W950(3)*nanmean(Org);...
    MAEs.W950(4)*nanmean(BC)];

end
