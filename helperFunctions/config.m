function [configs]=config()

configs.N=25; %choose the number of mesh points in x and y
configs.axisSize=5; %how many wavelengths of the F1 pattern you will see
configs.makeF2Fig=1; %turn on or off to build the figure with all simulations.
configs.showProgress=1; %turn on or off progress photos

loadedParams=exist('param_F1','var');

if loadedParams==0
    error(['Looks like you havent correctly loaded a parameter set. ' ...
        'Make sure you have all the parameters you need, then try again.'])
end

end