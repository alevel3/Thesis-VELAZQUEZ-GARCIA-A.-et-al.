%%% FUNCTION TO CALCULATE MEE OF 4 species: NH42SO4,NHNO3,Org,BC %%%
%%fitlm mdl

%%% USING MEEs => Cheng et al 2008

function[MEE,Bext_PM1_calc_MLR,RC] = Extinction_Cheng(NH42SO4,NH4NO3,Org,BC,EXT_525)%,NO2,Rayleigh);
%% Xdata -> Input
% isgood = ~(isnan(SSA_P5)|isnan(SSA_C11));
%isgood= ~(isnan(NH42SO4)|isnan(NH4NO3)|isnan(NH4Cl)|isnan(Org)|isnan(BC_ff)|isnan(BC_wb)|isnan(EXT_525)|isnan(NO2)|isnan(Rayleigh));
%Org = ACSM.Org;
xdata=[NH42SO4 NH4NO3 Org BC];
%% Absorption in PM1 (Bap -> ydata)
ydata = [EXT_525];
%ydata = [EXT_450,EXT_525,EXT_635];
%% MEEs Cheng et al 2008
MEE = [7.5 8.0 7.1 11.0];
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
%% Model using "regress function"
%X = [ones(size(xdata(:,1))) xdata(:,1) xdata(:,2) xdata(:,2) xdata(:,4) xdata(:,5) xdata(:,6) ];

% [b,bint,r,rint,stats] = regress(ydata,X)
% clc
% %[___] = regress(ydata,X,alpha)
% R2=r
% MEE_525 = b
% CI = bint
% STATS =stats

%% Model using "fitlm"
% mdl = fitlm(xdata,ydata,'Intercept',false)%% same coefficients than using lqlin
% %mdl = fitlm(xdata,ydata,'RobustOpts','on')
% %options = fitoptions('Lower', [0 0 0 0]);
% mdl = fitlm(xdata,ydata,'Intercept',false,'RobustOpts','on')
% %mdl = fitlm(xdata,ydata)%,'y~x1+x2+x3+x4+x5+x6')%% robust fit: mdl = fitlm(ingredients,heat,'RobustOpts','on')
% MEE_525 = table2struct(mdl.Coefficients(:,1),'ToScalar',true);
% MEE = MEE_525.Estimate
% ci = coefCI(mdl)%%(,aplha)
% %[p,F,r] = coefTest(mdl) %% Linear hypothesis test on linear regression model coefficients
% %S = table2struct(T,'ToScalar',true)
% %mdl = fitlm(xdata,ydata,'RobustOpts','on')
%% Model using lsqlin ydata(:,1) ydata(:,2) ydata(:,3) => 3wv
% lb=[0,0,0,0,0,0];
% ub=[];
% x0=[0,0,0,0,0,0];
% [x,resnorm,residual,exitflag,output,lambda] = lsqlin(xdata,ydata,[],[],[],[],lb,ub,x0);
% MEE = x
% yfit = xdata*MEE;
% resid = ydata - yfit;
% ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CInt = ci'
%% Compute Mean Absolute Error Using Cross-Validation
%rng('default') % For reproducibility
% function errors = regf(X1train,X2train,ytrain,X1test,X2test,ytest)
% tbltrain = table(X1train,X2train,ytrain, ...
%     'VariableNames',{'Acceleration','Displacement','Weight'});
% tbltest = table(X1test,X2test,ytest, ...
%     'VariableNames',{'Acceleration','Displacement','Weight'});
% mdl = fitlm(tbltrain,'Weight ~ Acceleration + Displacement');
% yfit = predict(mdl,tbltest);
% MAE = mean(abs(yfit-tbltest.Weight));
% adjMAE = MAE/range(tbltest.Weight);
% errors = [MAE adjMAE];
% end
% values = crossval(@regf,xdata(:,1),xdata(:,2),ydata)
% mean(values)

%% Specify Response and Predictor Variables for Linear Model
% mdl = fitlm(hospital,'interactions','ResponseVar','Weight',...
%     'PredictorVars',{'Sex','Age','Smoker'},...
%     'CategoricalVar',{'Sex','Smoker'}) weight


%%METLABT EXAMPLE
% load carsmall
% X = [Weight,Horsepower,Acceleration];%%%THere is not "ordenada al origen (1)"
% mdl = fitlm(X,MPG)


%% Model using "lsqcurvefit"
% %% %decay model:y=b(1)*xdata(:,1)+b(2)*xdata(:,2)+b(3)*xdata(:,3)+b(4)*xdata(:,4)+b(5)*xdata(:,5)+b(6)*xdata(:,6)
%     %% %Parameters bounded%
%     %ub=[3.8,5.1,4.3,4.5,11,6.3];%% Maximum values of each specie: NH4SO4, NH4NO3, NH4Cl, Org, Bc, BrC%
%     lb=[0,0,0,0,0,0]; %%NO PARTICLES NO Bap%%
%     %lb=[0.08,0.11,0.1,0.1,0.2,0.12];
%     ub=[];
%     %% %%%CReate the model
%     %Absportion -> MAE %%%replace b by  'MSE'|'MEE'|'MAE'
%     fun = @(MEE,xdata)MEE(1)*xdata(:,1)+MEE(2)*xdata(:,2)+MEE(3)*xdata(:,3)+MEE(4)*xdata(:,4)+MEE(5)*xdata(:,5)+MEE(6)*xdata(:,6);
%     %No Scattering
%     % fun = @(b,xdata)b(1)*xdata(:,1)+b(2)*xdata(:,2)+b(3)*xdata(:,3)+b(4)*xdata(:,4)+b(5)*xdata(:,5)+b(6)*xdata(:,6)
%     %Extinction -> MEE %%%replace b by  'MSE'|'MEE'|'MAE'
%     %fun = @(MEE,xdata)MEE(1)*xdata(:,1)+MEE(2)*xdata(:,2)+MEE(3)*xdata(:,3)+MEE(4)*xdata(:,4)+MEE(5)*xdata(:,5)+MEE(6)*xdata(:,6);
%     %% %%%Initian guess
%     x0=[0,0,0,0,0,0];%%%No particles% No exctinction
%     % x0=[0,0,0,0,0,0,0]%% b(0)
%      
% %% %% Solve the bounded fitting problem %WITH SCAT
% % %%%Calculating MSE%%
% %period_5min stand for Sca. 450, 525, 635
% % fun=@(coe,X)(coe(1)*X(:,1)+coe(2)*X(:,2)+coe(3)*X(:,3)+coe(4)*X(:,4)+coe(5)*X(:,5)+coe(6)*X(:,6));
% %MEE_450=lsqcurvefit(fun,x0,xdata,ydata(:,1),lb,ub);
% [x,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
% conf = nlparci(x,residual,'jacobian',jacobian)
% MEE_525 = [x]

%% Retriving Ext
%%%fitlm
% Ext_NH42SO4 = MEE_525.Estimate(1).*NH42SO4;
% Ext_NH4NO3 = MEE_525.Estimate(2).*NH4NO3;
% %Ext_NH4Cl = MEE_525.Estimate(4).*NH4Cl;
% Ext_Org = MEE_525.Estimate(3).*Org;
% Ext_BC = MEE_525.Estimate(4).*BC;
% %Ext_BrC = MEE_525.Estimate(7).*BC_wb;
% Ext_NO2 = 0.33.*NO2;

%%lsqlin 
Ext_NH42SO4 = MEE(1).*NH42SO4;
Ext_NH4NO3 = MEE(2).*NH4NO3;
Ext_Org = MEE(3).*Org;
Ext_BC = MEE(4).*BC;

%Ext_NO2 = 0.33.*NO2;

%%calcule
Bext_PM1_calc_MLR = nansum([Ext_NH42SO4 Ext_NH4NO3 Ext_Org Ext_BC],2);%Ext_NO2  Rayleigh

RC = [MEE(1)*nanmean(NH42SO4);...
    MEE(2)*nanmean(NH4NO3);...
    MEE(3)*nanmean(Org)
    MEE(4)*nanmean(BC)];
end