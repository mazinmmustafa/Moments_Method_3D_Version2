//
#include "projection.hpp"

projection_1d_para prjection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> p){
    vector_t<real_t> v21=v2-v1;
    const real_t alpha=v21*(p-v1)/(v21*v21);
    vector_t<real_t> p_0=v1+alpha*v21;
    projection_1d_para para;
    para.p_0 = p_0;
    vector_t<real_t> P_0=para.p_0-p;
    para.P_0 = mag(P_0);
    para.P_0_unit = unit(P_0);
    para.l_unit = unit(v2-v1);
    para.l_m = (v1-para.p_0)*para.l_unit;
    para.l_p = (v2-para.p_0)*para.l_unit;
    para.P_m = mag(v1-p);
    para.P_p = mag(v2-p);
    return para;
}