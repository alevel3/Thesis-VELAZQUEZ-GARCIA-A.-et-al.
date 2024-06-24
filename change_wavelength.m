function [data_lambda3] = change_wavelength(data_lambda1,AE1,lambda1,lambda3)
%AE1=-log10(data_lambda1/data_lambda2)/(log10(lambda1/lambda2));
data_lambda3 = data_lambda1.*(lambda3/lambda1).^(-AE1);
end

