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
simpar.c  = Z1.NumObservations-2;
% d
simpar.d1 = Z1.SSE;
simpar.d2 = Z2.SSE;
simpar.d3 = Z3.SSE;

% Simulate
simStruct = struct;
nSim      = 1e04;

% Vector of IRF
IRF     = NaN(nSim,1);
% Warning: Next line makes the code slow. 
IRh_h   = NaN(nSim,500);

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
        % Long term IRF
        IRF(s)    = auxFunctions.calculateIRF(simMoments);
        % Warning: Loop makes code very slow
        % Transition
        for i = 1:500
            % fill the vector
            IRh_h(s,i) = auxFunctions.calculateIRF(simMoments,i);
        end

        
    end
end
toc
% Clear auxiliary variables
N      = s;
clear s reg i r

% Count surviving draws
idx    = ~isnan(IRF);
numObs = sum(idx);


% Remove NaN observations
IRF     = IRF(idx);
IRh_h   = IRh_h(idx,:);

% Sort observations
[IRF,idx]   = sort(IRF);
IRh_h       = IRh_h(idx,:);

%% Plots of IRF  

% Trim outliers (keep distribution between 0.1 and 99.9 percent)
p1      = round(length(IRF)*1e-03);
IRF     = IRF((p1+1):end-p1);
IRh_h   = IRh_h((p1+1):end-p1,:);

% Calculate percentiles
p10 = round(length(IRF)*1e-01);
p50 = round(length(IRF)*5e-01);
p90 = round(length(IRF)*9e-01);

% Plot impulse response for those percentiles
figure;
hold on
plot(IRh_h(p10,:),'Color',c.maroon)
plot(IRh_h(p50,:),'Color',c.nvyBlue)
plot(IRh_h(p90,:),'Color',c.dkGreen)
plot([1,500],[IRF(p10),IRF(p10)],'--k','LineWidth',1.2)
plot([1,500],[IRF(p50),IRF(p50)],'--k','LineWidth',1.2)
plot([1,500],[IRF(p90),IRF(p90)],'--k','LineWidth',1.2)
xlabel('Period')
ylabel('IRF')
legend('Percentile 10','Median','Percentile 90','Location','northwest','box','off')
ylim(max(ylim(),0))
export_fig('Pset06Q1d1','-pdf','-transparent'); 

%% Plot distribution of IRF
% Fit Kernel to have a smooth distribution
pd = fitdist(IRF,'Lognormal');
x  = linspace(0.01,3,1e03);
y  = pdf(pd,x);

% Plot distribution 
dist = figure;
hold on
patch([[0,5],flip(x)],[[0,0],flip(y)],[0.12 0.18 0.48],'FaceAlpha',0.27)
pl10 = plot([IRF(p10),IRF(p10)],ylim(),'Color',c.maroon,'LineWidth',1.2);
pl50 = plot([IRF(p50),IRF(p50)],ylim(),'Color',c.nvyBlue,'LineWidth',1.2);
pl90 = plot([IRF(p90),IRF(p90)],ylim(),'Color',c.dkGreen,'LineWidth',1.2);
xlim([0,0.5])
legend([pl10,pl50,pl90],'Percentile 10','Median','Percentile 90','Location','northeast','box','off')
ylabel('Density')
xlabel('Long Term Impulse Response')
saveas(dist,'Pset06Q1d2.eps','epsc')


