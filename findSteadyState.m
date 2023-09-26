function [u1_ss,u2_ss,v1_ss,v2_ss] = findSteadyState(param)
    dt=0.1; %step size to be used
    N=1000000; %max steps to run
    u1=zeros(1,N); %create empty vector for u
    u1(1)=0.1; %initial value u
    u2=zeros(1,N); %create empty vector for u
    u2(1)=0.1; %initial value u
    v1=zeros(1,N); %create empty vector for v
    v1(1)=0.1; %initial value v
    v2=zeros(1,N); %create empty vector for v
    v2(1)=0.1; %initial value v
    epsilon = 10e-15; %threshold to determine when steady-state is reached
    
    for n=1:N-1
        [f1,f2,g1,g2] = findUdot(u1(n),u2(n),v1(n),v2(n),param);
        u1(n+1)=u1(n)+dt*f1;
        u2(n+1)=u2(n)+dt*f2;
        v1(n+1)=v1(n)+dt*g1;
        v2(n+1)=v2(n)+dt*g2;
        if abs(f1) < abs(epsilon*u1(n)) && abs(f2) < abs(epsilon*u2(n)) && abs(g1) < abs(epsilon*v1(n)) && abs(g2) < abs(epsilon*v2(n))
            break
        end
    end 
    u1_ss=u1(n); u2_ss=u2(n); v1_ss=v1(n); v2_ss=v2(n);
end

function [f1,f2,g1,g2] = findUdot(u1,u2,v1,v2,param)
    au1 = param.alpha_u1; bu1 = param.beta_u1; gu1 = param.gamma_u1;
    au2 = param.alpha_u2; bu2 = param.beta_u2; gu2 = param.gamma_u2;
    av1 = param.alpha_v1; bv1 = param.beta_v1; gv1 = param.gamma_v1;
    av2 = param.alpha_v2; bv2 = param.beta_v2; gv2 = param.gamma_v2;
    ku1u1 = param.ku1u1; ku1u2 = param.ku1u2; ku1v1 = param.ku1v1; ku1v2 = param.ku1v2;
    ku2u1 = param.ku2u1; ku2u2 = param.ku2u2; ku2v1 = param.ku2v1; ku2v2 = param.ku2v2;
    kv1u1 = param.kv1u1; kv1u2 = param.kv1u2; kv1v1 = param.kv1v1; kv1v2 = param.kv1v2;
    kv2u1 = param.kv2u1; kv2u2 = param.kv2u2; kv2v1 = param.kv2v1; kv2v2 = param.kv2v2;
    
    
    f1 = -bu1*u1+au1+gu1*((ku1u1*u1+ku1u2*u2)/(1+ku1u1*u1+ku1u2*u2+ku1v1*v1+ku1v2*v2))^2;
    f2 = -bu2*u2+au2+gu2*((ku2u1*u1+ku2u2*u2)/(1+ku2u1*u1+ku2u2*u2+ku2v1*v1+ku2v2*v2))^2;
    g1 = -bv1*v1+av1+gv1*((kv1u1*u1+kv1u2*u2)/(1+kv1u1*u1+kv1u2*u2+kv1v1*v1+kv1v2*v2))^2;
    g2 = -bv2*v2+av2+gv2*((kv2u1*u1+kv2u2*u2)/(1+kv2u1*u1+kv2u2*u2+kv2v1*v1+kv2v2*v2))^2;
end