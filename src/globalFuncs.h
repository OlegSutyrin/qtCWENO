#ifndef qtCWENO_globalFuncs_H //include guard
#define qtCWENO_globalFuncs_H

//коэффициенты Рунге-Кутты
const double RKcoeff[RK_ORDER_MAX][RK_ORDER_MAX][RK_ORDER_MAX + 1] = { //[order-1][step][Qn], [..][..][RK_ORDER_MAX] - производные
    {
        {1.0, 0.0, 0.0, 1.0}, //1 порядок: explicit Euler
        {},
        {}
    },
    {
        {1.0, 0.0, 0.0, 1.0}, //2 Heun's method
        {0.5, 0.5, 0.0, 0.5}, //(Butcher tableau: a21 = 1 | b1 = 1/2, b2 = 1/2)
        {}
    },
    {
        {1.0, 0.0, 0.0, 1.0}, //3 Strong stability preserving (TVD) Runge-Kutta (SSPRK3)
        {0.75, 0.25, 0.0, 0.25},
        {1.0 / 3.0, 0.0, 2.0 / 3.0, 2.0 / 3.0} //(Butcher tableau: a21 = 1 | a31 = 1/4, a32 = 1/4 | b1 = 1/6, b2 = 1/6, b3 = 2/3)
    }
};

//глобальные функции (названия с большой буквы)

#endif


