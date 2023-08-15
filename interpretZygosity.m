function [param,M_eff] = interpretZygosity(U1,U2,V1,V2,D_u)
    param.D_u = D_u;
    
    param.alpha_u1=U1.alpha_u;
    param.beta_u1=U1.beta_u;
    param.gamma_u1=U1.gamma_u;
    param.ku1u1=U1.kuu;
    
    param.alpha_u2=U2.alpha_u;
    param.beta_u2=U2.beta_u;
    param.gamma_u2=U2.gamma_u;
    param.ku2u2=U2.kuu;
    
    param.alpha_v1=V1.alpha_v;
    param.beta_v1=V1.beta_v;
    param.gamma_v1=V1.gamma_v;
    param.kv1v1=V1.kvv;
    
    param.alpha_v2=V2.alpha_v;
    param.beta_v2=V2.beta_v;
    param.gamma_v2=V2.gamma_v;
    param.kv2v2=V2.kvv;
    
    M_eff=zeros(4,1);
    
    if U1.type == V1.type
        param.ku1v1=U1.kuv;
        param.kv1u1=V1.kvu;
    else
        param.ku1v1=U1.kuv_off;
        param.kv1u1=V1.kvu_off;
    end
    
    if U1.type == V2.type
        param.ku1v2=U1.kuv;
        param.kv2u1=V2.kvu;
    else
        param.ku1v2=U1.kuv_off;
        param.kv2u1=V2.kvu_off;
    end
    
    if U2.type == V1.type
        param.ku2v1=U2.kuv;
        param.kv1u2=V1.kvu;
    else
        param.ku2v1=U2.kuv_off;
        param.kv1u2=V1.kvu_off;
    end
    
    if U2.type == V2.type
        param.ku2v2=U2.kuv;
        param.kv2u2=V2.kvu;
    else
        param.ku2v2=U2.kuv_off;
        param.kv2u2=V2.kvu_off;
    end
    
    if U1.type == U2.type
        param.ku1u2=U1.kuu;
        param.ku2u1=U2.kuu;
    else
        param.ku1u2=U1.kuu_off;
        param.ku2u1=U2.kuu_off;
    end
    
    if V1.type == V2.type
        param.kv1v2=V1.kvv;
        param.kv2v1=V2.kvv;
    else
        param.kv1v2=V1.kvv_off;
        param.kv2v1=V2.kvv_off;
    end
    
    
    M_eff(1)=(param.ku1u1+param.ku1u2+param.ku2u1+param.ku2u2)/(2*U1.kuu+2*U2.kuu);
    M_eff(2)=(param.ku1v1+param.ku1v2+param.ku2v1+param.ku2v2)/(2*U1.kuv+2*U2.kuv);
    M_eff(3)=(param.kv1u1+param.kv1u2+param.kv2u1+param.kv2u2)/(2*V1.kvu+2*V2.kvu);
    M_eff(4)=(param.kv1v1+param.kv1v2+param.kv2v1+param.kv2v2)/(2*V1.kvv+2*V2.kvv);
    
end