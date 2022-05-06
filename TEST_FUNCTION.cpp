#include <cmath>
#include <iostream>

void ACKLEY_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void RASTRIGIN_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void ROSENBROCK_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (1,....,1)
void SPHERE_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Michalewicz_Function(double *x, double *f, int DIM); //global optima f(x*) = -1.8013 when DIM = 2
void Bent_Cigar_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Schaffer_F7_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Zakharov_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Griewank_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Schwefel_226_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (420.9687,....,420.9687)
void SUM_SQUARES_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void LEVY03_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (1,....,1)
void STYBLINSKI_TANG_Function(double *x, double *f, int DIM); //global optima f(x*) = -39.16599 * DIM at x* = (-2.903534,....,-2.903534)
void Salomon_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Whitley(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (1,....,1)
void Brown_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void Exponential_Function(double *x, double *f, int DIM); //global optima f(x*) = -1 at x* = (0,....,0)
void Rotated_Hyper_Ellipsoid_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at x* = (0,....,0)
void DIXON_PRICE_Function(double *x, double *f, int DIM); //global optima f(x*) = 0 at every x_i = pow(2, -((pow(2, i) - 2) / pow(2, i)))
void Trid_Function(double *x, double *f, int DIM); //global optima f(x*) = -DIM * (DIM + 4) * (DIM - 1) / 6 at every x_i = i * (DIM + 1 - i)
void CosineMixture(double *x, double *f, int DIM); //global optima f(x*) = -0.1 * DIM at x* = (0,....,0)
void Perm_beta_Function(double *x, double *f, int DIM); //global optima f(x*) = 0, at x* = (1, 1/2 , ... , 1/DIM)
void Shubert_N4(double *x, double *f, int DIM); //global optima f(x*) ≈ −25.740858 when DIM = 2
void Powell_Singular(double *x, double *f, int DIM); //global optima f(x*) ≈ 0, at x* = (3,−1,0,1,··· ,3,−1,0,1)
void Powell_Sum(double *x, double *f, int DIM);//global optima f(x*) ≈ 0
void Ripple_1(double *x, double *f, int DIM);
void Ripple_25(double *x, double *f, int DIM);
void Deb_1(double *x, double *f, int DIM);
void Deb_3(double *x, double *f, int DIM);
void Csendes(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Alpine_1(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Alpine_2(double *x, double *f, int DIM);//global optima f(x*) = -pow(2.808, DIM) at x* = (7.917,....,7.917)
void Chung_Reynolds(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Pathological(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Pinter(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Qing(double *x, double *f, int DIM);//global optima f(x*) = 0  at every x_i = +sqrt(i) or -sqrt(i)
void Sargan(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Streched_V_Sine_Wave(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Trigonometric_1(double *x, double *f, int DIM);//global optima f(x*) = 0 at x* = (0,....,0)
void Trigonometric_2(double *x, double *f, int DIM);//global optima f(x*) = 1 at x* = (0.9,....,0.9)

void Unimodal_Range(double *max, double *min, int dim, int function_num);
void Unimodal_Function(double *x, double *f, int dim, int function_num);
void Multimodal_Range(double *max, double *min, int dim, int function_num);
void Multimodal_Function(double *x, double *f, int dim, int function_num);

void Set_Range(double *max, double *min, int dim, int function_num){
    if(function_num >= 1 && function_num <= 14){
        Unimodal_Range(max, min, dim, function_num);
    }
    else if(function_num >= 15 && function_num <= 40){
        Multimodal_Range(max, min, dim, function_num - 14);
    }
    else{
        std::cout << "Invaild Function Number" << std::endl;
        abort();
    }
}

void ALL_TEST_FUNCTION(double *x, double *f, int dim, int function_num){
    if(function_num >= 1 && function_num <= 14){
        Unimodal_Function(x, f, dim, function_num);
    }
    else if(function_num >= 15 && function_num <= 40){
        Multimodal_Function(x, f, dim, function_num - 14);
    }
    else{
        std::cout << "Invaild Function Number" << std::endl;
        abort();
    }
}


/*
*Unimodal Function:
*No.1 ROSENBROCK_Function /done
*No.2 SPHERE_Function /done
*No.3 Bent_Cigar_Function /done
*No.4 Zakharov_Function /done
*No.5 SUM_SQUARES_Function /done
*No.6 Brown_Function /done
*No.7 Exponential_Function /done
*No.8 Rotated_Hyper_Ellipsoid_Function /done
*No.9 DIXON_PRICE_FUNCTION /uncertain
*No.10 Trid_Function /done
*NO.11 Powell_Singular /done
*NO.12 Powell_Sum /done
*NO.13 Chung_Reynolds /done
*NO.14 Streched_V_Sine_Wave /done
*/
void Unimodal_Range(double *max, double *min, int dim, int function_num){
    switch (function_num){
        case 1:
            *max = 10.0;
            *min = -5.0;
            break;
        case 2:
            *max = 5.12;
            *min = -5.12;
            break;
        case 3:
            *max = 100.0;
            *min = -100.0;
            break;
        case 4:
            *max = 10.0;
            *min = -5.0;
            break;
        case 5:
            *max = 10.0;
            *min = -10.0;
            break;
        case 6:
            *max = 4.0;
            *min = -1.0;
            break;
        case 7:
            *max = 1.0;
            *min = -1.0;
            break;
        case 8:
            *max = 65.536;
            *min = -65.536;
            break;
        case 9:
            *max = 10.0;
            *min = -10.0;
            break;
        case 10:
            *max = pow(dim, 2);
            *min = -pow(dim, 2);
            break;
        case 11:
            *max = 5.0;
            *min = -4.0;
            break;
        case 12:
            *max = 1.0;
            *min = -1.0;
            break;
        case 13:
            *max = 100.0;
            *min = -100.0;
            break;
        case 14:
            *max = 10.0;
            *min = -10.0;
            break;
        default:
            break;
    }
}

void Unimodal_Function(double *x, double *f, int dim, int function_num){
    switch(function_num){
        case 1:
            ROSENBROCK_Function(x, f, dim);
            break;
        case 2:
            SPHERE_Function(x, f, dim);
            break;
        case 3:
            Bent_Cigar_Function(x, f, dim);
            break;
        case 4:
            Zakharov_Function(x, f, dim);
            break;
        case 5:
            SUM_SQUARES_Function(x, f, dim);
            break;
        case 6:
            Brown_Function(x, f, dim);
            break;
        case 7:
            Exponential_Function(x, f, dim);
            break;
        case 8:
            Rotated_Hyper_Ellipsoid_Function(x, f, dim);
            break;
        case 9:
            DIXON_PRICE_Function(x, f, dim);
            break;
        case 10:
            Trid_Function(x, f, dim);
            break;
        case 11:
            Powell_Singular(x, f, dim);
            break;
        case 12:
            Powell_Sum(x, f, dim);
            break;
        case 13:
            Chung_Reynolds(x, f, dim);
            break;
        case 14:
            Streched_V_Sine_Wave(x, f, dim);
            break;
        default:
            break;
    }
}

/*
*Multimodal Function:
*No.1 RASTRIGIN_Function /done
*No.2 ACKLEY_Function /done
*No.3 CosineMixture /done
*N0.4 Schaffer_F7 /done
*No.5 Griewank_Function /done
*No.6 Schwefel 2.26 /done
*No.7 LEVY03_Function /done
*No.8 STYBLINSKI_TANG_Function /done
*No.9 Salomon_Function /done
*No.10 Whitley_Function /done
*No.11 Perm 0,d beta Function /done
*No.12 Shubert N.4 /done
*No.13 Michalewicz_Function/ done
*NO.14 Ripple_1 /done
*NO.15 Ripple_25 /done
*NO.16 Deb1 /done
*NO.17 Deb3 /done
*NO.18 Csendes /done
*NO.19 Alpine_1 /done
*NO.20 Alpine_2 /done
*NO.21 Pathological /done
*NO.22 Pinter /done
*NO.23 Qing /done
*No.24 Sargan /done
*NO.25 Trigonometric1 /done
*NO.26 Trigonometric2 /done
*/

void Multimodal_Range(double *max, double *min, int dim, int function_num){
    switch (function_num){
        case 1:
            *max = 5.12;
            *min = -5.12;
            break;
        case 2:
            *max = 32.768;
            *min = -32.768;
            break;
        case 3:
            *max = 1.0;
            *min = -1.0;
            break;
        case 4:
            *max = 100.0;
            *min = -100.0;
            break;
        case 5:
            *max = 600.0;
            *min = -600.0;
            break;
        case 6:
            *max = 500.0;
            *min = -500.0;
            break;
        case 7:
            *max = 5.0;
            *min = -5.0;
            break;
        case 8:
            *max = 5.0;
            *min = -5.0;
            break;
        case 9:
            *max = 100.0;
            *min = -100.0;
            break;
        case 10:
            *max = 10.24;
            *min = -10.24;
            break;
        case 11:
            *max = dim;
            *min = -dim;
            break;
        case 12:
            *max = 10.0;
            *min = -10.0;
            break;
        case 13:
            *max = M_PI;
            *min = 0.0;
            break;
        case 14:
            *max = 1.0;
            *min = 0.0;
            break;
        case 15:
            *max = 1.0;
            *min = 0.0;
            break;
        case 16:
            *max = 1.0;
            *min = -1.0;
            break;
        case 17:
            *max = 1.0;
            *min = -1.0;
            break;
        case 18:
            *max = 1.0;
            *min = -1.0;
            break;
        case 19:
            *max = 10.0;
            *min = -10.0;
            break;
        case 20:
            *max = 10.0;
            *min = 0.0;
            break;
        case 21:
            *max = 100.0;
            *min = -100.0;
            break;
        case 22:
            *max = 10.0;
            *min = -10.0;
            break;
        case 23:
            *max = 500.0;
            *min = -500.0;
            break;
        case 24:
            *max = 100.0;
            *min = -100.0;
            break;
        case 25:
            *max = M_PI;
            *min = 0.0;
            break;
        case 26:
            *max = 500.0;
            *min = -500.0;
            break;
        default:
            break;
    }
}

void Multimodal_Function(double *x, double *f, int dim, int function_num){
    switch(function_num){
        case 1:
            RASTRIGIN_Function(x, f, dim);
            break;
        case 2:
            ACKLEY_Function(x, f, dim);
            break;
        case 3:
            CosineMixture(x, f, dim);
            break;
        case 4:
            Schaffer_F7_Function(x, f, dim);
            break;
        case 5:
            Griewank_Function(x, f, dim);
            break;
        case 6:
            Schwefel_226_Function(x, f, dim);
            break;
        case 7:
            LEVY03_Function(x, f, dim);
            break;
        case 8:
            STYBLINSKI_TANG_Function(x, f, dim);
            break;
        case 9:
            Salomon_Function(x, f, dim);
            break;
        case 10:
            Whitley(x, f, dim);
            break;
        case 11:
            Perm_beta_Function(x, f, dim);
            break;
        case 12:
            Shubert_N4(x, f, dim);
            break;
        case 13:
            Michalewicz_Function(x, f, dim);
            break;
        case 14:
            Ripple_1(x, f, dim);
            break;
        case 15:
            Ripple_25(x, f, dim);
            break;
        case 16:
            Deb_1(x, f, dim);
            break;
        case 17:
            Deb_3(x, f, dim);
            break;
        case 18:
            Csendes(x, f, dim);
            break;
        case 19:
            Alpine_1(x, f, dim);
            break;
        case 20:
            Alpine_2(x, f, dim);
            break;
        case 21:
            Pathological(x, f, dim);
            break;
        case 22:
            Pinter(x, f, dim);
            break;
        case 23:
            Qing(x, f, dim);
            break;
        case 24:
            Sargan(x, f, dim);
            break;
        case 25:
            Trigonometric_1(x, f, dim);
            break;
        case 26:
            Trigonometric_2(x, f, dim);
            break;
        default:
            break;
    }
}


void ACKLEY_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    double t2 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
        t2 += cos(2 * M_PI * x[i]);
    }
    t1 /= DIM;
    t2 /= DIM;

    f[0] = (-20.0) * exp((-0.2) * sqrt(t1)) - exp(t2) + 20.0 + exp(1);
}

void RASTRIGIN_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2) - (10 * cos(2 * M_PI * x[i]));
    }
    f[0] = 10.0 * DIM + t1;
}

void ROSENBROCK_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM - 1; ++i){
        t1 += 100.0 * pow(x[i + 1] - pow(x[i], 2), 2) + pow(x[i] - 1, 2);
    }
    f[0] = t1;
}

void SPHERE_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
    }
    f[0] = t1;
}
//uncertain
void Michalewicz_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double a1 = sin(x[i]);
        double a2 = pow(sin(((i + 1) * x[i] * x[i]) / M_PI), 2 * 10);
        t1 += a1 * a2;
    }
    f[0] = -t1;
}

void Bent_Cigar_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 1; i < DIM; ++i){
        t1 += pow(x[i], 2);
    }
    f[0] = pow(x[0], 2) + pow(10.0, 6) * t1;
}

void Schaffer_F7_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM - 1; ++i){
        double t2 = sqrt(pow(x[i], 2) + pow(x[i + 1], 2));
        t1 += sqrt(t2) + sqrt(t2) * pow(sin(50.0 * pow(t2, 0.2)), 2);
    }
    f[0] = pow(1.0 / (DIM - 1) * t1, 2);
}

void Zakharov_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    double t2 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
        t2 += 0.5 * x[i] * (i + 1);
    }
    f[0] = t1 + pow(t2, 2) + pow(t2, 4);
}

void Griewank_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    double t2 = 1.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
        t2 *= cos((x[i] / sqrt(i + 1)));
    }
    f[0] = t1 / 4000.0 - t2 + 1;
}

void Schwefel_226_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += x[i] * sin(sqrt(fabs(x[i])));
    }
    f[0] = 418.9829 * DIM - t1;
}

void SUM_SQUARES_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += (i + 1) * pow(x[i], 2);
    }
    f[0] = t1;
}

void LEVY03_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    double Wd = 1.0 + ((x[DIM - 1] - 1.0) / 4.0);
    for(int i = 0; i < DIM - 1; ++i){
        double Wi = 1.0 + ((x[i] - 1.0) / 4.0);
        t1 += pow(Wi - 1, 2) * (1.0 + 10 * pow(sin(M_PI * Wi + 1.0), 2)) + pow(Wd - 1, 2) * (1.0 + pow(sin(2 * M_PI * Wd), 2));
    }
    double W1 = 1.0 + ((x[0] - 1.0) / 4.0);
    f[0] = pow(sin(W1 * M_PI), 2) + t1;
}

void STYBLINSKI_TANG_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 4) - 16 * pow(x[i], 2) + 5 * x[i];
    }
    f[0] = 0.5 * t1;
}

void Salomon_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
    }
    f[0] =  -cos(2.0 * M_PI * sqrt(t1)) + 0.1 * sqrt(t1) + 1.0;
}

void Whitley(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 0.0;
        for(int j = 0; j < DIM; ++j){
            double a = pow(100.0 * pow(pow(x[i], 2) - x[j], 2) + pow(1.0 - x[j], 2), 2) / 4000.0;
            double b = cos(100.0 * pow(pow(x[i], 2) - x[j], 2) + pow(1.0 - x[j], 2));
            t2 += a - b + 1.0;
        }
        t1 += t2;
    }
    f[0] = t1;
}

void Brown_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM - 1; ++i){
        t1 += (pow(pow(x[i], 2), pow(x[i + 1], 2) + 1) + pow(pow(x[i + 1], 2), pow(x[i], 2) + 1));
    }
    f[0] = t1;
}

void Exponential_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
    }
    f[0] = -exp(-0.5 * t1);
}

void Rotated_Hyper_Ellipsoid_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 0.0;
        for(int j = 0; j <= i; ++j){
            t2 += x[j];
        }
        t1 += pow(t2, 2);
    }
    f[0] = t1;
}

void DIXON_PRICE_Function(double *x, double *f, int DIM){//*
    double t1 = 0.0;
    for(int i = 1; i < DIM; ++i){
        t1 += (i + 1) * pow(2 * pow(x[i], 2) - x[i - 1], 2);
    }
    f[0] = pow(x[0] - 1, 2) + t1;
}

void Trid_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    double t2 = 0.0;
    for(int i = 0; i < DIM; ++i){
        if(i != 0){
            t2 += x[i] * x[i - 1];
        }
        t1 += pow(x[i] - 1, 2);
    }
    f[0] = t1 - t2;
}

void CosineMixture(double *x, double *f, int DIM){
    double t1 = 0.0;
    double t2 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += cos(5.0 * M_PI * x[i]);
        t2 += pow(x[i], 2);
    }

    f[0] = -(0.1 * t1 - t2);
}

void Perm_beta_Function(double *x, double *f, int DIM){
    // Beta default value : 10
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 0.0;
        for(int j = 0; j < DIM; ++j){
            t2 += ((j + 1) + 10) * (pow(x[j], i + 1) - pow(j + 1, -(i + 1)));
        }
        t1 += pow(t2, 2);
    }
    f[0] = t1;
}

void Shubert_N4(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 0.0;
        for(int j = 1; j < 6; ++j){
            t2 += j * cos((j + 1) * x[i] + j);
        }
        t1 += t2;
    }
    f[0] = t1;
}

void Ripple_1(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
       t1 += -exp(-2.0 * log(2.0 * pow((x[i] - 0.1) / 0.8, 2))) * (pow(sin(5.0 * M_PI * x[i]), 6) + 0.1 * pow(cos(500.0 * M_PI * x[i]), 2));
    }
    f[0] = t1;
}

void Ripple_25(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
       t1 += -exp(-2.0 * log(2.0 * pow((x[i] - 0.1) / 0.8, 2))) * (pow(sin(5.0 * M_PI * x[i]), 6));
    }
    f[0] = t1;
}

void Deb_1(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(sin(5.0 * M_PI * x[i]), 6);
    }
    f[0] = ((-1.0) / DIM) * t1;
}

void Deb_3(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(sin(5.0 * M_PI * (pow(x[i], 3.0 / 4.0) - 0.05)), 6);
    }
    f[0] = ((-1.0) / DIM) * t1;
}

void Csendes(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        if(x[i] == 0){
            t1 += pow(x[i], 6) * (2.0 + sin(0.0));
        }
        else{
            t1 += pow(x[i], 6) * (2.0 + sin(1.0 / x[i]));
        }
        
    }
    f[0] = t1;
}

void Alpine_1(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += fabs(x[i] * sin(x[i]) + 0.1 * x[i]);
    }
    f[0] = t1;
}

void Alpine_2(double *x, double *f, int DIM){
    double t1 = 1.0;
    for(int i = 0; i < DIM; ++i){
        t1 *= sqrt(x[i]) * sin(x[i]);
    }
    f[0] = -t1;
    
}

void Pathological(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM - 1; ++i){
        double temp = (pow(sin(sqrt(100.0 * pow(x[i], 2) + pow(x[i + 1], 2))), 2) - 0.5) / (1.0 + 0.001 * pow(pow(x[i] - 2.0 * x[i] * x[i + 1] + pow(x[i + 1], 2), 2), 2));
        t1 += 0.5 + temp;
    }
    f[0] = t1;
}

void Pinter(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double A;
        double B ;
        if(i == 0){
            A = x[DIM - 1] * sin(x[i] + sin(x[i + 1]));
            B = pow(x[DIM - 1], 2) - 2.0 * x[i] + 3.0 * x[i + 1] - cos(x[i]) + 1.0;
        }
        else if(i == DIM - 1){
            A = x[DIM - 1] * sin(x[i] + sin(x[0]));
            B = pow(x[i - 1], 2) - 2.0 * x[i] + 3.0 * x[0] - cos(x[i]) + 1.0;
        }
        else{
            A = x[i - 1] * sin(x[i] + sin(x[i + 1]));
            B = pow(x[i - 1], 2) - 2.0 * x[i] + 3.0 * x[i + 1] - cos(x[i]) + 1.0;
        }
        t1 += (i + 1) * pow(x[i], 2) + 20.0 * (i + 1) * pow(sin(A), 2) + (i + 1) * log10(1.0 + (i + 1) * pow(B, 2));
    }
    f[0] = t1;
}

void Qing(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(pow(x[i], 2) - (i + 1), 2);
    }
    f[0] = t1;
}

void Sargan(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 0.0;
        for(int j = 1; j < DIM; ++j){
            t2 += x[i] * x[j];
        }
        t1 += DIM * (pow(x[i], 2) + 0.4 * t2);
    }
    f[0] = t1;
}

void Trigonometric_1(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 0.0;
        for(int j = 0; j < DIM; ++j){
            t2 += cos(x[j]);
        }
        t1 += pow((DIM - t2 + (i + 1) * (1.0 - cos(x[i]) - sin(x[i]))), 2);
    }
    f[0] = t1;
}

void Trigonometric_2(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        double t2 = 8.0 * pow(sin(7.0 * pow(x[i] - 0.9, 2)), 2);
        double t3 = 6.0 * pow(sin(14.0 * pow(x[i] - 0.9, 2)), 2);
        t1 += (t2 + t3 + pow(x[i] - 0.9, 2));
    }
    f[0] = 1.0 + t1;
}

void Powell_Singular(double *x, double *f, int DIM){
    if(DIM / 4 > 0){
        double t1 = 0.0;
        for(int i = 0; i < DIM / 4; ++i){
            int c = 4 * (i + 1);
            t1 += (pow(x[c - 4] + 10.0 * x[c - 3], 2) + 5.0 * pow(x[c - 2] - x[c - 1], 2) + pow(x[c - 3] - x[c - 2], 4) + 10.0 * pow(x[c - 4] - x[c - 1], 4));
        }
        f[0] = t1;
    }
    else{
        std::cout << "illegal DIM" << std::endl;
        abort();
    }
}

void Powell_Sum(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(fabs(x[i]), i + 2);
    }
    f[0] = t1;
}

void Chung_Reynolds(double *x, double *f, int DIM){
    double t1 = 0;
    for(int i = 0; i < DIM; ++i){
        t1 += pow(x[i], 2);
    }
    f[0] = pow(t1, 2);
}

void Streched_V_Sine_Wave(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM - 1; ++i){
        double t2 = pow(pow(x[i + 1], 2) + pow(x[i], 2), 0.25);
        double t3 = pow(sin(50.0 * pow(pow(x[i + 1], 2) + pow(x[i], 2), 0.1)), 2) + 0.1;
        t1 += t2 * t3;
    }
    f[0] = t1;
}