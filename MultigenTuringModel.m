% Written by E. Simmons, 2023
% Note: You must upload a correctly formatted parameter set before running
% the code!

clearvars -except param_F1
close all
clc

%call configuration function
configVars=config();

loadedParams=exist('param_F1','var');
if loadedParams==0
    error(['Looks like you havent correctly loaded a parameter set. ' ...
        'Make sure you have all the parameters you need, then try again.'])
end

%call function to run all simulations in F2 generation
[configVars,U_data]=runF2Simulations(configVars,param_F1);

%optially, call function to output and save a figure
if configs.makeF2Fig==1
    [F2Figure]=buildF2Fig(configs,U_data);
    savefig(['F2SimulationsSz',num2str(configs.N),'.fig'])
end

%save simulation results
save(['UdataSz',num2str(configs.N),'.dat'],'U_data','-ascii')
save(['UdataSz',num2str(configs.N),'.mat'],'U_data')