%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   M-function for fitting probabilistic models to Malawi hazard values   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Probability model fitting is based on the least squares method

% Candidate models
% 1) Lognormal
% 2) Gumbel
% 3) Frechet
% 4) Weibull

function [para,coef,corr] = ModelFit_MalawiSeismicHazard(F,H,model_option)

H = H(:);

if model_option == 1 % Lognormal distribution
    
    % z = norminv(F) and x = log(y)
    % z = -para(1)/para(2) + 1/para(2)*x
    % coef(1) = -para(1)/para(2)
    % coef(2) = 1/para(2)
    [coef,~,~,~,corr] = regress(norminv(F),[ones(length(H),1) log(H)]); 
    para(1) = -coef(1)/coef(2);
    para(2) = 1/coef(2);

elseif model_option == 2 % Gumbel distribution
    
    % z = -log(-log(F)) and x = x
    % z = -para(1)*para(2) + para(2)*x
    % coef(1) = -para(1)*para(2)
    % coef(2) = para(2)
    [coef,~,~,~,corr] = regress(-log(-log(F)),[ones(length(H),1) H]); 
    para(1) = -coef(1)/coef(2);
    para(2) = coef(2);

elseif model_option == 3 % Frechet distribution
    
    % z = -log(-log(F)) and x = log(x)
    % z = -log(para(1))*para(2) + para(2)*x
    % coef(1) = -log(para(1))*para(2)
    % coef(2) = para(2)
    [coef,~,~,~,corr] = regress(-log(-log(F)),[ones(length(H),1) log(H)]); 
    para(1) = exp(-coef(1)/coef(2));
    para(2) = coef(2);

elseif model_option == 4 % Weibull distribution
    
    % z = log(-log(1-F)) and x = log(x)
    % z = -log(para(1))*para(2) + para(2)*x
    % coef(1) = -log(para(1))*para(2)
    % coef(2) = para(2)
    [coef,~,~,~,corr] = regress(log(-log(1-F)),[ones(length(H),1) log(H)]); 
    para(1) = exp(-coef(1)/coef(2));
    para(2) = coef(2);

end















