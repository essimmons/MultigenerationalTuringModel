% Written by E. Simmons, 2023
% Note: You must upload a correctly formatted parameter set before running
% the code!

clearvars -except param_F1
close all
clc

configs.N=25; %choose the number of mesh points in x and y
configs.axisSize=5; %how many wavelengths of the F1 pattern you will see
configs.makeF2Fig=1; %turn on or off to build the figure with all simulations.
configs.showProgress=1; %turn on or off progress photos

loadedParams=exist('param_F1','var');

if loadedParams==0
    error(['Looks like you havent correctly loaded a parameter set. ' ...
        'Make sure you have all the parameters you need, then try again.'])
end

[configs,U_data]=runF2Simulations(configs,param_F1);

if configs.makeF2Fig==1
    [F2Figure]=buildF2Fig(configs,U_data);
    savefig(['F2SimulationsSz',num2str(configs.N),'.fig'])
end
save(['UdataSz',num2str(configs.N),'.dat'],'U_data','-ascii')
save(['UdataSz',num2str(configs.N),'.mat'],'U_data')