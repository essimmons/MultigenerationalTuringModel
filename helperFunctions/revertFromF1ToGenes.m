function [UA,UB,VA,VB]=revertFromF1ToGenes(param)
    UA.alpha_u=param.alpha_u1;
    UA.beta_u=param.beta_u1;
    UA.gamma_u=param.gamma_u1;
    UA.kuu=param.ku1u1;
    UA.kuv=param.ku1v1;
    UA.kuu_off=param.ku1u2;
    UA.kuv_off=param.ku1v2;
    UA.type="A";

    VA.alpha_v=param.alpha_v1;
    VA.beta_v=param.beta_v1;
    VA.gamma_v=param.gamma_v1;
    VA.kvu=param.kv1u1;
    VA.kvv=param.kv1v1;
    VA.kvu_off=param.kv1u2;
    VA.kvv_off=param.kv1v2;
    VA.type="A";


    UB.alpha_u=param.alpha_u2;
    UB.beta_u=param.beta_u2;
    UB.gamma_u=param.gamma_u2;
    UB.kuu=param.ku2u2;
    UB.kuv=param.ku2v2;
    UB.kuu_off=param.ku2u1;
    UB.kuv_off=param.ku2v1;
    UB.type="B";

    VB.alpha_v=param.alpha_v2;
    VB.beta_v=param.beta_v2;
    VB.gamma_v=param.gamma_v2;
    VB.kvu=param.kv2u2;
    VB.kvv=param.kv2v2;
    VB.kvu_off=param.kv2u1;
    VB.kvv_off=param.kv2v1;
    VB.type="B";
end