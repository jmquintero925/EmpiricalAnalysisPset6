%% Import data
clear all;
close all;
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["period", "logC_diff", "logElogC", "logDlogC"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("Problem_16.csv", opts);

period      = flip(tbl.period);
logC_diff   = flip(tbl.logC_diff);
logElogC    = flip(tbl.logElogC);
logDlogC    = flip(tbl.logDlogC);

clear opts tbl
clc;

%% Create matrix of controls

% Create index for c
idxC = repmat((2:279)',1,4)+(0:3);
% Create index for the rest
idx   = repmat((2:279)',1,5)+(0:4);
% Create matrix of controls
Xt = [logC_diff(idxC),logElogC(idx),logDlogC(idx)];

clear idxC idx

%% Linear models

% Index for controls 
idx = (1:size(Xt,1))';

Z1 = fitlm(Xt,logC_diff(idx));
Z2 = fitlm([logC_diff(idx),Xt],logElogC(idx));
Z3 = fitlm([logC_diff(idx),logElogC(idx),Xt],logDlogC(idx));

% Beta
struct.beta1 = table2array(Z1.Coefficients(:,1));
struct.beta2 = table2array(Z2.Coefficients(:,1));
struct.beta3 = table2array(Z3.Coefficients(:,1));
% Sigma
struct.sigma1   = (1/Z1.NumObservations)*Z1.SSE;
struct.sigma2   = (1/Z2.NumObservations)*Z2.SSE;
struct.sigma3   = (1/Z3.NumObservations)*Z3.SSE;
% Number of observations
struct.n    = Z1.NumObservations;

%% Calculate the long term IRF from the estimated values

% Calculate moments
moments = auxFunctions.calculateMatrices(struct);
% Add the var structure
moments = auxFunctions.calculateVAR(moments);
% Calculate long term shock
irf = auxFunctions.calculateIRF(moments);
% Full impulse response 
irfTrans = NaN(1000,1);
for i = 1:length(irfTrans)
    % fill the vector
    irfTrans(i) = auxFunctions.calculateIRF(moments,i);
end


%% Plot Estimated IRF 
save = false;
try
   global c
   myColors();
   myCols = true;
catch
    warning("Using standard Matlab palette for plot")
end

% Plot the IRF
if(myCols)
    figure;
    hold on
    plot(irfTrans,'Color',c.maroon)
    plot([1,length(irfTrans)],[irf,irf],'--k')
    xlabel('Period')
    ylabel('IRF')
    legend('Impulse Response','Long Term Effect','Location','southeast','box','off')
else
    figure;
    hold on
    plot(irfTrans,'Color')
    plot([1,length(irfTrans)],[irf,irf],'--k')
    xlabel('Period')
    ylabel('IRF')
    legend('Impulse Response','Long Term Effect','Location','southeast','box','off')
end














