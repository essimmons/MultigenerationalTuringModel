function [configs,U_data]=runF2Simulations(configs,param_F1)

if configs.showProgress==1
    f20=figure('Name','Progress Figure');
end

%parameter setup
[UA,UB,VA,VB]=revertFromF1ToGenes(param_F1);
[param_PA]=interpretZygosity(UA,UA,VA,VA);
[param_PB]=interpretZygosity(UB,UB,VB,VB);
param_PA.D_u=param_F1.D_u; param_PA.D_v=param_F1.D_v;
param_PB.D_u=param_F1.D_u; param_PB.D_v=param_F1.D_v;

%initial calculations
J_F1=findJacobian(param_F1);
[~,~,~,~,configs.w] = findCriticalk2(J_F1,param_F1);

%Parent and F1 simulations
[U_PA,~,~,~]=ImplicitRxnDiff2D(param_PA,configs);
[U_PB,~,x,y]=ImplicitRxnDiff2D(param_PB,configs);
[U_F1,~,~,~]=ImplicitRxnDiff2D(param_F1,configs);

%% F2 simulations

%gather data for spatial grid, parent and F1 simulations
[X,Y]=meshgrid(x,y);
U_data=zeros((length(x))^2,14);
U_data(:,1)=reshape(X,length(U_data),1);
U_data(:,2)=reshape(Y,length(U_data),1);
U_data(:,3)=reshape(U_PA,length(U_data),1);
U_data(:,4)=reshape(U_PB,length(U_data),1);
U_data(:,5)=reshape(U_F1,length(U_data),1);

%prepare for F2 simulations
U=[UA UB];
V=[VA VB];
F2Count=1;
paramF2=struct('alpha_u1',{},'beta_u1',{},'gamma_u1',{}, ...
    'alpha_u2',{},'beta_u2',{},'gamma_u2',{},...
    'alpha_v1',{},'beta_v1',{},'gamma_v1',{},...
    'alpha_v2',{},'beta_v2',{},'gamma_v2',{},...
    'ku1u1',{},'ku1u2',{},'ku1v1',{},'ku1v2',{},...
    'ku2u1',{},'ku2u2',{},'ku2v1',{},'ku2v2',{},...
    'kv1u1',{},'kv1u2',{},'kv1v1',{},'kv1v2',{},...
    'kv2u1',{},'kv2u2',{},'kv2v1',{},'kv2v2',{},...
    'D_u',0.1,'D_v',{});

%iterate through all allele combos, run simulations on each
for i=1:2
    for j=i:2
        for k=1:2
            for l=k:2
                [paramTmp]=interpretZygosity(U(i),U(j),V(k),V(l));
                paramTmp.D_u=param_F1.D_u; paramTmp.D_v=param_F1.D_v;
                paramF2(F2Count)=paramTmp;
                if configs.showProgress==1
                    figure(f20)
                end
                [U_F2,~,~,~]=ImplicitRxnDiff2D(paramTmp,configs);
                U_data(:,F2Count+5)=reshape(U_F2,length(U_data),1);
                F2Count=F2Count+1;
            end
        end
    end
end