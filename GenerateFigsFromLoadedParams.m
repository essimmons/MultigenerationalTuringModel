% Written by E. Simmons, 2023
% Note: You must upload a correctly formatted parameter set before running
% the code!


loadedParams=exist('param_F1')&&exist('param_PA')&&exist('param_PB')...
    &&exist('UA')&&exist('UB')&&exist('VA')&&exist('VB');

if loadedParams==0
    error(['Looks like you havent correctly loaded a parameter set. ' ...
        'Make sure you have all the parameters you need, then try again.'])
end

N=25; %choose the number of mesh points in x and y
axisSize=5; %how many wavelengths of the F1 pattern you will see
makeF2Fig=1; %turn on or off to build the figure with all simulations.


D_u=param_F1.D_u;
D_v=param_F1.D_v;
J_F1=findJacobian(param_F1);
[k2,realEigsAp,k2crit,realEigCrit,w] = findCriticalk2(J_F1,param_F1,D_v);
        
             
        [U_PA,V_PA,~,~]=ImplicitRxnDiff2D(param_PA,D_u,D_v,w,N,axisSize);
        
        [U_PB,V_PB,x,y]=ImplicitRxnDiff2D(param_PB,D_u,D_v,w,N,axisSize);
        
        [U_F1,V_F1,~,~]=ImplicitRxnDiff2D(param_F1,D_u,D_v,w,N,axisSize);
        figParams.x=x; figParams.y=y; figParams.w=w;
        %}
        [X,Y]=meshgrid(x,y);
        U_data=zeros((length(x))^2,14);
        U_data(:,1)=reshape(X,length(U_data),1);
        U_data(:,2)=reshape(Y,length(U_data),1);
        U_data(:,3)=reshape(U_PA,length(U_data),1);
        U_data(:,4)=reshape(U_PB,length(U_data),1);
        U_data(:,5)=reshape(U_F1,length(U_data),1);
        %}
        %%
    
ucaxisMax=ceil(max(U_F1,[],'all')); 
%vcaxisMax=ceil(max(V_F1,[],'all'));
myAxes=[0 5*w 0 5*w];

    if makeF2Fig==true
        U=[UA UB];
        V=[VA VB];


        F2labels={'F2aaaa', 'F2aaab', 'F2aabb', 'F2abaa', 'F2abab', 'F2abbb', 'F2bbaa', 'F2bbab', 'F2bbbb'};
        F2Count=1;
        paramF2=struct('alpha_u1',{},'beta_u1',{},'gamma_u1',{}, ...
            'alpha_u2',{},'beta_u2',{},'gamma_u2',{},...
            'alpha_v1',{},'beta_v1',{},'gamma_v1',{},...
            'alpha_v2',{},'beta_v2',{},'gamma_v2',{},...
            'ku1u1',{},'ku1u2',{},'ku1v1',{},'ku1v2',{},...
            'ku2u1',{},'ku2u2',{},'ku2v1',{},'ku2v2',{},...
            'kv1u1',{},'kv1u2',{},'kv1v1',{},'kv1v2',{},...
            'kv2u1',{},'kv2u2',{},'kv2v1',{},'kv2v2',{},...
            'D_u',0.1,'F2_label',{});

        f3=figure('Name','F2 Simulation');
        f3.Position=[1201 925 1000 1000];
        pos=zeros(9,4);
        pos(:,3)=0.27;
        pos(:,4)=0.25;

        for i=1:3
            pos(3*i,1)=0.69;
            pos(3*i-1,1)=0.37;
            pos(3*i-2,1)=0.05;

            pos(i,2)=0.69;
            pos(i+3,2)=0.37;
            pos(i+6,2)=0.05;
        end

        
        figure(f3)
       


        for i=1:2
            for j=i:2
                for k=1:2
                    for l=k:2
                        %if F2Count ~= 5 && F2Count ~= 1 && F2Count ~=9
                            
                            [paramTmp,M_eff]=interpretZygosity(U(i),U(j),V(k),V(l),D_u);
                            paramTmp.F2_label=F2labels(F2Count);
                            paramF2(F2Count)=paramTmp;

                            figure(20)
                            [U_F2,V_F2,x,y]=ImplicitRxnDiff2D(paramTmp,D_u,D_v,w,N,axisSize);
                            U_data(:,F2Count+5)=reshape(U_F2,length(U_data),1);


                            figure(f3)
                            subplot('Position',pos(F2Count,:))
                            surf(x,y,U_F2,"EdgeColor","none")
                            view(2)
                            colorbar
                            %colormap bone
                            clim([0 ucaxisMax])
                            axis([0 5*w 0 5*w])
                            title(paramF2(F2Count).F2_label)
                            %}
                        %end

                        F2Count=F2Count+1;

                    end
                end
            end
        end
        
        save(['UdataSz',num2str(N),'.dat'],'U_data','-ascii')
        save(['UdataSz',num2str(N),'.mat'],'U_data')

    end
