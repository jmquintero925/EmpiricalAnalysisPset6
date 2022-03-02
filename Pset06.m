%% Import data
clear all;
close all; clc;
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

% 
save = true;
try
   global c
   myColors();
   myCols = true;
catch
    warning("Using standard Matlab palette for plot")
end

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
if(save)
    export_fig('Pset06Q1','-pdf','-transparent'); 
end

%% Simulate

rng(1)
%---Store simulation parameters----
% Mean
simpar.b1 = struct.beta1;
simpar.b2 = struct.beta2;
simpar.b3 = struct.beta3;
% Matrix Lambda
cons      = ones(size(Xt,1),1);
simpar.L1 = auxFunctions.calculateGamma([cons,Xt],0);
simpar.L2 = auxFunctions.calculateGamma([cons,logC_diff(idx),Xt],0);
simpar.L3 = auxFunctions.calculateGamma([cons,logC_diff(idx),logElogC(idx),Xt],0);
% C
simpar.c  = Z1.NumObservations;
% d
simpar.d1 = Z1.SSE;
simpar.d2 = Z2.SSE;
simpar.d3 = Z3.SSE;

% Simulate
simStruct = struct;
nSim      = 1e02;

% Vector of IRF
IRF     = NaN(nSim,1);
minEig  = NaN(nSim,1);

tic
for s = 1:nSim
    % Simulate the beta and sigma for each regressions
    for reg = 1:3
        % Which regression
        r = string(reg);
        % Simulate sigma
        simStruct.(strcat('sigma',r))   = 1/gamrnd(simpar.c/2,2/simpar.(strcat('d',r)));
        % Get the inverse of Lambda
        L = inv(simpar.(strcat('L',r)));
        % Account for numerical errors
        L = (L+L')/2;
        % Simulate Beta
        simStruct.(strcat('beta',r))    = mvnrnd(simpar.(strcat('b',r)),L*simStruct.(strcat('sigma',r)))';
    end
    % Calculate the model primitives for simulated values
    simMoments = auxFunctions.calculateMatrices(simStruct);
    % Add the var structure
    simMoments = auxFunctions.calculateVAR(simMoments);
    % Check if the matrix is stable
    if(all(abs(eig(simMoments.A))<1))
        minEig(s) = min(abs(eig(simMoments.A)));
        IRF(s)    = auxFunctions.calculateIRF(simMoments);
        
    end
end
toc








