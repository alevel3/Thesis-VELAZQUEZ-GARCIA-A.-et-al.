%%% FUNCTION TO CALCULATE MSE OF 4 species: NH42SO4,NHNO3,Org,BC %%%
%%% NO2 -> absorband gass
%%fitlm mdl_B mdl_G mdl_R
function[MSE_B,CI_B,MSE_G,CI_G,MSE_R,CI_R,Scat_CB,Scat_CG,Scat_CR,RC_B,RC_G,RC_R,res_B,res_G,res_R,hi] = Scattering_3s(NH42SO4,NH4NO3,Org,Scat_B,Scat_G,Scat_R)%,Rayleigh_B,Rayleigh_G,Rayleigh_R)
%Only BC6 without NH4CHL
%% Xdata -> Input
%Org = ACSM.Org
%xdata=[NH42SO4 NH4NO3 Org BC];% BC_ff BC_wb];
xdata=[NH42SO4 NH4NO3 Org];
%% Absorption in PM1 (Bap -> ydata)
ydata = [Scat_B,Scat_G,Scat_R];
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

sum(leverage(xdata))
hi=leverage(xdata);

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
%mdl = fitlm(xdata,ydata,'RobustOpts','on')
%dlm = fitlm(xdata,ydata,'Intercept',false);
%dlm = fitlm(xdata,ydata,'y~x1+x2+x3+x4+x5+x6');
%% Model using "lsqlin"
% x = lsqlin(C,d,A,b)
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
% x = lsqlin(problem)
% [x,resnorm,residual,exitflag,output,lambda] = lsqlin(___)
% [wsout,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b,Aeq,beq,lb,ub,ws)
%fun = @(MSE,xdata)MSE(1)*xdata(:,1)+MSE(2)*xdata(:,2)+MSE(3)*xdata(:,3);
% %fx=@(x)lsqlin(xdata,ydata(:,1),[],[],[],[],lb,ub,x0);

%%---intercept ---
% X = [ones(size(xdata(:,1))) xdata(:,1) xdata(:,2) xdata(:,3)];
% lb=[0,0,0,0];
% ub=[];
% x0=[0,0,0,0];

%%% --- No intercept ---%
lb=[0,0,0];
ub=[];
x0=[0,0,0];
%% Scat_CB,
%%fitlm (it works)
% mdl_B = fitlm(xdata,ydata(:,1),'Intercept',false)
% %mdl_B = fitlm(xdata,ydata(:,1),'RobustOpts','on')
% %mdl_B = fitlm(xdata,ydata(:,1))%% mdl = fitlm(xdata,ydata,'RobustOpts','on')
% MSE_450 = table2struct(mdl_B.Coefficients(:,1),'ToScalar',true);
% MSE_B = MSE_450.Estimate
% CI_B = coefCI(mdl_B)%%(,aplha)

%%% ---lsqlin intercept --- %
% [x,resnorm,res_B,exitflag,output,lambda] = lsqlin(X,ydata(:,1),[],[],[],[],lb,ub,x0);
% MSE_B = x
% yfit = X*MSE_B;
% resid = ydata(:,1) - yfit;
% ci = bootci(100,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_B = ci';
% P=polyfitn(ydata(:,1),X*MSE_B,1);
% P.R2
% P.Coefficients

%%%---lslin NO intercept --%
[x,resnorm,res_B,exitflag,output,lambda] = lsqlin(xdata,ydata(:,1),[],[],[],[],lb,ub,x0);
MSE_B = x
yfit = xdata*MSE_B;
resid = ydata(:,1) - yfit;
ci = bootci(100,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_B = ci';
P=polyfitn(ydata(:,1),xdata*x,1);
P.R2
P.Coefficients
%% Scat_CG

%%% ---lsqlin intercept --- %
% [x,resnorm,res_G,exitflag,output,lambda] = lsqlin(X,ydata(:,2),[],[],[],[],lb,ub,x0);
% MSE_G = x
% yfit = X*MSE_G;
% resid = ydata(:,2) - yfit;
% ci = bootci(100,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_G = ci';
% P=polyfitn(ydata(:,2),X*MSE_G,1);
% P.R2
% P.Coefficients

%%%---lslin NO intercept --%
[x,resnorm,res_G,exitflag,output,lambda] = lsqlin(xdata,ydata(:,2),[],[],[],[],lb,ub,x0);
MSE_G = x
yfit = xdata*MSE_G;
resid = ydata(:,2) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_G = ci'
P=polyfitn(ydata(:,2),xdata*x,1)
P.R2
P.Coefficients
%% Scat_CR
%%fitlm(it works)
% mdl_R = fitlm(xdata,ydata(:,3),'Intercept',false)
% %mdl_R = fitlm(xdata,ydata(:,3),'RobustOpts','on')
% %mdl_R = fitlm(xdata,ydata(:,3))
% MSE_730 = table2struct(mdl_R.Coefficients(:,1),'ToScalar',true);
% MSE_R = MSE_730.Estimate
% CI_R = coefCI(mdl_R)%%(,aplha)

%%% ---lsqlin intercept --- %
% [x,resnorm,res_R,exitflag,output,lambda] = lsqlin(X,ydata(:,3),[],[],[],[],lb,ub,x0);
% MSE_R = x
% yfit = X*MSE_R;
% resid = ydata(:,3) - yfit;
% ci = bootci(100,{@(bootr)lsqlin(X,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
% CI_R = ci';
% P=polyfitn(ydata(:,3),X*MSE_R,1);
% P.R2
% P.Coefficients

%%%---lslin NO intercept --%
[x,resnorm,res_R,exitflag,output,lambda] = lsqlin(xdata,ydata(:,3),[],[],[],[],lb,ub,x0);
MSE_R = x
yfit = xdata*MSE_R;
resid = ydata(:,3) - yfit;
ci = bootci(1000,{@(bootr)lsqlin(xdata,yfit+bootr,[],[],[],[],lb,ub,x0),resid},'Type','normal');
CI_R = ci'
P=polyfitn(ydata(:,3),xdata*x,1)
P.R2
P.Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S = table2struct(T,'ToScalar',true)
%mdl = fitlm(xdata,ydata,'RobustOpts','on')

%%METLABT EXAMPLE
% load carsmall
% X = [Weight,Horsepower,Acceleration];%%%THere is not "ordenada al origen (1)"
% mdl = fitlm(X,MPG)

%% Model using "lsqcurvefit"
 %decay model:y=b(1)*xdata(:,1)+b(2)*xdata(:,2)+b(3)*xdata(:,3)+b(4)*xdata(:,4)+b(5)*xdata(:,5)+b(6)*xdata(:,6)
%% %Parameters bounded%
%ub=[0.397,0.035,x,0.112,4.604,0.246];%% Maximum values of each specie: NH4SO4, NH4NO3, NH4Cl, Org, Bc, BrC%
% lb=[0,0,0,0,0,0]; %%NO PARTICLES NO SCATTERING%%
% %lb=[0.036,0.047,x,0.0,0.099,0.0];
% ub=[];
% %% %%%CReate the model
% %Absportion -> MAE %%%replace b by  'MSE'|'MEE'|'MAE'
% fun = @(MAE,xdata)MAE(1)*xdata(:,1)+MAE(2)*xdata(:,2)+MAE(3)*xdata(:,3)+MAE(4)*xdata(:,4)+MAE(5)*xdata(:,5)+MAE(6)*xdata(:,6);
% %No Scattering
% % fun = @(b,xdata)b(1)*xdata(:,1)+b(2)*xdata(:,2)+b(3)*xdata(:,3)+b(4)*xdata(:,4)+b(5)*xdata(:,5)+b(6)*xdata(:,6)
% %Extinction -> MEE %%%replace b by  'MSE'|'MEE'|'MAE'
% %fun = @(MEE,xdata)MEE(1)*xdata(:,1)+MEE(2)*xdata(:,2)+MEE(3)*xdata(:,3)+MEE(4)*xdata(:,4)+MEE(5)*xdata(:,5)+MEE(6)*xdata(:,6);
% %% %%%Initian guess
% x0=[0,0,0,0,0,0];%%%No particles% No exctinction
% % x0=[0,0,0,0,0,0,0]%% b(0)
% %% %% Solve the bounded fitting problem %WITH SCAT+ABS
% % %%%Calculating MEE%%
% %period_5min stand for Sca. 450, 525, 635
% % fun=@(coe,X)(coe(1)*X(:,1)+coe(2)*X(:,2)+coe(3)*X(:,3)+coe(4)*X(:,4)+coe(5)*X(:,5)+coe(6)*X(:,6));
% %MEE_450=lsqcurvefit(fun,x0,xdata,ydata(:,1),lb,ub);
% 
% [x,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(fun,x0,xdata,ydata(:,1),lb,ub);
% conf = nlparci(x,residual,'jacobian',jacobian)
% MSE_450 = [x]'
% 
% [x,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(fun,x0,xdata,ydata(:,2),lb,ub);
% conf = nlparci(x,residual,'jacobian',jacobian)
% MSE_525 = [x]'
% 
% [x,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(fun,x0,xdata,ydata(:,3),lb,ub);
% conf = nlparci(x,residual,'jacobian',jacobian)
% MSE_635 = [x]'

% MSE_450=lsqcurvefit(fun,x0,xdata,ydata(:,1),lb,ub);
% MSE_525=lsqcurvefit(fun,x0,xdata,ydata(:,2),lb,ub);
% MSE_635=lsqcurvefit(fun,x0,xdata,ydata(:,3),lb,ub);
% MSE_450=[MSE_450]'
% MSE_525=[MSE_525]'
% MSE_635=[MSE_635]'
% Species={'NH42SO2' 'NH42NO3' 'NH42Chl' 'Org' 'BC' 'BrC'}
% Table_MEE=table(Species,MEE_450,MEE_525,MEE_635)
%MSEs=[MSE_450 MSE_525 MSE_635];
%% Retrieving Scatt
%% Scat_CB
SB_NH42SO4 = MSE_B(1).*NH42SO4;
SB_NH4NO3 = MSE_B(2).*NH4NO3;
%SB_NH4Cl = MSE_B(4).*NH4Cl;
SB_Org = MSE_B(3).*Org;
%SB_BC = MSE_B(4).*BC;
%SB_BC = MSE_B(6).*BC_ff;
%SB_BrC = MSE_B(6).*BC_wb;

Scat_CB = sum([SB_NH42SO4 SB_NH4NO3 SB_Org],2);%Rayleigh_B / before 27 June nansum
RC_B = [MSE_B(1)*nanmean(NH42SO4);...
    MSE_B(2)*nanmean(NH4NO3);...
    MSE_B(3)*nanmean(Org)];

%% Scat_CG
SG_NH42SO4 = MSE_G(1).*NH42SO4;
SG_NH4NO3 = MSE_G(2).*NH4NO3;
%SG_NH4Cl = MSE_G(4).*NH4Cl;
SG_Org = MSE_G(3).*Org;
%SG_BC = MSE_G(4).*BC;
%SG_BC = MSE_G(6).*BC_ff;
%SG_BrC = MSE_G(6).*BC_wb;

Scat_CG = sum([SG_NH42SO4 SG_NH4NO3 SG_Org],2);%Rayleigh_G
RC_G = [MSE_G(1)*nanmean(NH42SO4);...
    MSE_G(2)*nanmean(NH4NO3);...
    MSE_G(3)*nanmean(Org)];

%% Scat_CR
SR_NH42SO4 = MSE_R(1).*NH42SO4;
SR_NH4NO3 = MSE_R(2).*NH4NO3;
%SR_NH4Cl = MSE_R(4).*NH4Cl;
SR_Org = MSE_R(3).*Org;
%SR_BC = MSE_R(4).*BC;
%SR_BC = MSE_R(6).*BC_ff;
%SR_BrC = MSE_R(6).*BC_wb;

Scat_CR = sum([SR_NH42SO4 SR_NH4NO3 SR_Org],2);%Rayleigh_R
RC_R = [MSE_R(1)*nanmean(NH42SO4);...
    MSE_R(2)*nanmean(NH4NO3);...
    MSE_R(3)*nanmean(Org)];

end