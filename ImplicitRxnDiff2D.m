function [U_tot,V_tot,x,y] = ImplicitRxnDiff2D(param,configs)
D_u=param.D_u; D_v=param.D_v;
au1 = param.alpha_u1; bu1 = param.beta_u1; gu1 = param.gamma_u1;
au2 = param.alpha_u2; bu2 = param.beta_u2; gu2 = param.gamma_u2;
av1 = param.alpha_v1; bv1 = param.beta_v1; gv1 = param.gamma_v1;
av2 = param.alpha_v2; bv2 = param.beta_v2; gv2 = param.gamma_v2;
ku1u1 = param.ku1u1; ku1u2 = param.ku1u2; ku1v1 = param.ku1v1; ku1v2 = param.ku1v2;
ku2u1 = param.ku2u1; ku2u2 = param.ku2u2; ku2v1 = param.ku2v1; ku2v2 = param.ku2v2;
kv1u1 = param.kv1u1; kv1u2 = param.kv1u2; kv1v1 = param.kv1v1; kv1v2 = param.kv1v2;
kv2u1 = param.kv2u1; kv2u2 = param.kv2u2; kv2v1 = param.kv2v1; kv2v2 = param.kv2v2;

N=configs.N;
w=configs.w;
axisSize=configs.axisSize;
prog=configs.showProgress;

%%Define spatial grid
xmax=axisSize*w;
dx=xmax/(N-1);
x=0:dx:xmax;
y=0:dx:xmax;

%%Define time grid
dt=20;
tmax=2e6;
t=0:dt:tmax;

epsilon=1e-10;
maxDev=0.045;

[u1_ss,u2_ss,v1_ss,v2_ss] = findSteadyState(param);
U1=u1_ss+maxDev*u1_ss*randn(N);
k=U1<0;
U1(k)=0;
U2=u2_ss+maxDev*u2_ss*rand(N);
k=U2<0;
U2(k)=0;
V1=v1_ss+maxDev*v1_ss*rand(N);
k=V1<0;
V1(k)=0;
V2=v2_ss+maxDev*v2_ss*rand(N);
k=V2<0;
V2(k)=0;
U_tot=U1+U2;
V_tot=V1+V2;

%Generate tridiagonal matrices A and A_rhs
r_u=D_u*dt/(2*dx^2);
r_v=D_v*dt/(2*dx^2);
while r_v>0.5
    dt=0.5*dt;
    t=0:dt:tmax;
    r_u=D_u*dt/(2*dx^2);
    r_v=D_v*dt/(2*dx^2);
end

e1=ones(N,1);
Bu_a=[-r_u*e1 (1+2*r_u)*e1 -r_u*e1];
Bv_a=[-r_v*e1 (1+2*r_v)*e1 -r_v*e1];
Lux=spdiags(Bu_a, [-1 0 1], N, N);
Lvx=spdiags(Bv_a, [-1 0 1], N, N);
Ix=speye(N);
A_u=kron(Ix,Lux);
A_v=kron(Ix,Lvx);

pos1=1:N:N^2;
for i = 1:length(pos1)
    A_u(pos1(i), pos1(i))=1+r_u;
    A_v(pos1(i), pos1(i))=1+r_v;
end

pos2=N:N:N^2;
for i=1:length(pos2)
    A_u(pos2(i),pos2(i))=1+r_u;
    A_v(pos2(i),pos2(i))=1+r_v;
end

e2=ones(N^2,1);
Bu_b=[r_u*e2 (1-r_u)*e2 r_u*e2];
Bv_b=[r_v*e2 (1-r_v)*e2 r_v*e2];
A_uRhsDiff=spdiags(Bu_b,[-N 0 N], N^2, N^2);
A_vRhsDiff=spdiags(Bv_b,[-N 0 N], N^2, N^2);

pos3 = N+1:N^2-N;
for i = 1:length(pos3)
    A_uRhsDiff(pos3(i),pos3(i))=1-2*r_u;
    A_vRhsDiff(pos3(i),pos3(i))=1-2*r_v;
end

for i=2:length(t)
    if prog==1
        if rem(i,1000)==0
            surf(x,y,U_tot,"EdgeColor","none")
            colorbar
            view(2)
            myAxes=[0 xmax 0 xmax];
            axis(myAxes)
            title({['Time is ', num2str(i*dt), ' time units'];'Activator'});
            drawnow();
        end
    end

    tmpU1 = reshape(flipud(U1)',N^2,1);
    tmpU2 = reshape(flipud(U2)',N^2,1);
    tmpV1 = reshape(flipud(V1)',N^2,1);
    tmpV2 = reshape(flipud(V2)',N^2,1);

    A_u1RhsRxn=(-bu1*tmpU1+au1+gu1*((ku1u1*tmpU1+ku1u2*tmpU2)./(1+ku1u1*tmpU1+ku1u2*tmpU2+ku1v1*tmpV1+ku1v2*tmpV2)).^2)*dt/2;
    A_u2RhsRxn=(-bu2*tmpU2+au2+gu2*((ku2u1*tmpU1+ku2u2*tmpU2)./(1+ku2u1*tmpU1+ku2u2*tmpU2+ku2v1*tmpV1+ku2v2*tmpV2)).^2)*dt/2;
    A_v1RhsRxn=(-bv1*tmpV1+av1+gv1*((kv1u1*tmpU1+kv1u2*tmpU2)./(1+kv1u1*tmpU1+kv1u2*tmpU2+kv1v1*tmpV1+kv1v2*tmpV2)).^2)*dt/2;
    A_v2RhsRxn=(-bv2*tmpV2+av2+gv2*((kv2u1*tmpU1+kv2u2*tmpU2)./(1+kv2u1*tmpU1+kv2u2*tmpU2+kv2v1*tmpV1+kv2v2*tmpV2)).^2)*dt/2;

    A_u1Rhs=A_uRhsDiff*tmpU1+A_u1RhsRxn;
    A_u2Rhs=A_uRhsDiff*tmpU2+A_u2RhsRxn;
    A_v1Rhs=A_vRhsDiff*tmpV1+A_v1RhsRxn;
    A_v2Rhs=A_vRhsDiff*tmpV2+A_v2RhsRxn;

    solU1 = A_u\A_u1Rhs;
    solU2 = A_u\A_u2Rhs;
    solV1 = A_v\A_v1Rhs;
    solV2 = A_v\A_v2Rhs;
    U1 = flipud(reshape(solU1,N,N)');
    U2 = flipud(reshape(solU2,N,N)');
    V1 = flipud(reshape(solV1,N,N)');
    V2 = flipud(reshape(solV2,N,N)');

    tmpU1 = reshape(flipud(U1),N^2,1);
    tmpU2 = reshape(flipud(U2),N^2,1);
    tmpV1 = reshape(flipud(V1),N^2,1);
    tmpV2 = reshape(flipud(V2),N^2,1);

    A_u1RhsRxn=(-bu1*tmpU1+au1+gu1*((ku1u1*tmpU1+ku1u2*tmpU2)./(1+ku1u1*tmpU1+ku1u2*tmpU2+ku1v1*tmpV1+ku1v2*tmpV2)).^2)*dt/2;
    A_u2RhsRxn=(-bu2*tmpU2+au2+gu2*((ku2u1*tmpU1+ku2u2*tmpU2)./(1+ku2u1*tmpU1+ku2u2*tmpU2+ku2v1*tmpV1+ku2v2*tmpV2)).^2)*dt/2;
    A_v1RhsRxn=(-bv1*tmpV1+av1+gv1*((kv1u1*tmpU1+kv1u2*tmpU2)./(1+kv1u1*tmpU1+kv1u2*tmpU2+kv1v1*tmpV1+kv1v2*tmpV2)).^2)*dt/2;
    A_v2RhsRxn=(-bv2*tmpV2+av2+gv2*((kv2u1*tmpU1+kv2u2*tmpU2)./(1+kv2u1*tmpU1+kv2u2*tmpU2+kv2v1*tmpV1+kv2v2*tmpV2)).^2)*dt/2;

    A_u1Rhs = A_uRhsDiff*tmpU1+A_u1RhsRxn;
    A_u2Rhs = A_uRhsDiff*tmpU2+A_u2RhsRxn;
    A_v1Rhs = A_vRhsDiff*tmpV1+A_v1RhsRxn;
    A_v2Rhs = A_vRhsDiff*tmpV2+A_v2RhsRxn;
    solU1 = A_u\A_u1Rhs;
    solU2 = A_u\A_u2Rhs;
    solV1 = A_v\A_v1Rhs;
    solV2 = A_v\A_v2Rhs;
    U1 = flipud(reshape(solU1,N,N));
    U2 = flipud(reshape(solU2,N,N));
    V1 = flipud(reshape(solV1,N,N));
    V2 = flipud(reshape(solV2,N,N));

    U_tot_tmp=U_tot;
    V_tot_tmp=V_tot;

    U_tot=U1+U2;
    V_tot=V1+V2;

    dt_changeU=abs(U_tot_tmp-U_tot);
    dt_changeV=abs(V_tot_tmp-V_tot);

    thresholdU=epsilon*U_tot;
    thresholdV=epsilon*V_tot;

    checkU=dt_changeU<thresholdU;
    checkV=dt_changeV<thresholdV;

    if all(checkU(:)==1) && all(checkV(:)==1)
        return
    end
end
end