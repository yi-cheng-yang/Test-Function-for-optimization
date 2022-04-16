#include <cmath>
#include <iostream>

void ACKLEY_Function(double *x, double *f, int DIM);
void RASTRIGIN_Function(double *x, double *f, int DIM);
void ROSENBROCK_Function(double *x, double *f, int DIM);
void SPHERE_Function(double *x, double *f, int DIM);
void Michalewicz_Function(double *x, double *f, int DIM); 
void Bent_Cigar_Function(double *x, double *f, int DIM);
void Schaffer_F7_Function(double *x, double *f, int DIM); 
void Zakharov_Function(double *x, double *f, int DIM);
void Griewank_Function(double *x, double *f, int DIM);
void Schwefel_226_Function(double *x, double *f, int DIM);
void SUM_SQUARES_Function(double *x, double *f, int DIM);
void LEVY03_Function(double *x, double *f, int DIM);
void LEVY05_Function(double *x, double *f, int DIM);
void STYBLINSKI_TANG_Function(double *x, double *f, int DIM);
void Salomon_Function(double *x, double *f, int DIM);
void Whitley(double *x, double *f, int DIM);
void Brown_Function(double *x, double *f, int DIM);
void Exponential_Function(double *x, double *f, int DIM);
void Rotated_Hyper_Ellipsoid_Function(double *x, double *f, int DIM);
void DIXON_PRICE_Function(double *x, double *f, int DIM);
void Trid_Function(double *x, double *f, int DIM);

void Unimodal_Range(double *max, double *min, int dim, int function_num);
void Unimodal_Function(double *x, double *f, int dim, int function_num);
void Multimodal_Range(double *max, double *min, int dim, int function_num);
void Multimodal_Function(double *x, double *f, int dim, int function_num);

void Set_Range(double *max, double *min, int dim, int Land_Structure, int function_num){
    if(Land_Structure == 0){
        Unimodal_Range(max, min, dim, function_num);
    }
    else{
        Multimodal_Range(max, min, dim, function_num);
    }
}

void ALL_TEST_FUNCTION(double *x, double *f, int dim, int function_num){
    if(function_num >= 1 && function_num <= 10){
        Unimodal_Function(x, f, dim, function_num);
    }
    else if(function_num >= 11 && function_num <= 21){
        Multimodal_Function(x, f, dim, function_num - 10);
    }
    else{
        std::cout << "Invaild Function Number" << std::endl;
        abort();
    }
}


/*
*Unimodal Function:
*No.1 ROSENBROCK_Function
*No.2 SPHERE_Function
*No.3 Bent_Cigar_Function
*No.4 Zakharov_Function
*No.5 SUM_SQUARES_Function
*No.6 Brown_Function
*No.7 Exponential_Function
*No.8 Rotated_Hyper_Ellipsoid_Function
*No.9 DIXON_PRICE_FUNCTION
*No.10 Trid_Function
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
        default:
            break;
    }
}

/*
*Multimodal Function:
*No.1 RASTRIGIN_Function
*No.2 ACKLEY_Function
*No.3 Michalewicz_Function
*N0.4 Schaffer_F7
*No.5 Griewank_Function
*No.6 Schwefel 2.26
*No.7 LEVY03_Function
*No.8 LEVY05_Function
*No.9 STYBLINSKI_TANG_Function
*No.10 Salomon_Function
*No.11 Whitley_Function
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
            *max = 5.12;
            *min = -5.12;
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
            *max = 10.0;
            *min = -10.0;
            break;
        case 8:
            *max = 2.0;
            *min = -2.0;
            break;
        case 9:
            *max = 5.0;
            *min = -5.0;
            break;
        case 10:
            *max = 100.0;
            *min = -100.0;
            break;
        case 11:
            *max = 10.24;
            *min = -10.24;
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
            Michalewicz_Function(x, f, dim);
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
            LEVY05_Function(x, f, dim);
            break;
        case 9:
            STYBLINSKI_TANG_Function(x, f, dim);
            break;
        case 10:
            Salomon_Function(x, f, dim);
            break;
        case 11:
            Whitley(x, f, dim);
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

    f[0] = (-20) * exp((-0.2) * sqrt(t1)) - exp(t2) + 20 + exp(1);
}

void RASTRIGIN_Function(double *x, double *f, int DIM){

}

void ROSENBROCK_Function(double *x, double *f, int DIM){

}

void SPHERE_Function(double *x, double *f, int DIM){

}

void Michalewicz_Function(double *x, double *f, int DIM){

}

void Bent_Cigar_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 1; i < DIM; ++i){
        t1 += x[i] * x[i];
    }
    f[0] = x[0] * x[0] + pow(10.0, 6) * t1;
}

void Schaffer_F7_Function(double *x, double *f, int DIM){

}

void Zakharov_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    double t2 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += x[i] * x[i];
        t2 += 0.5 * x[i] * i;
    }
    f[0] = t1 + pow(t2, 2) + pow(t2, 4);
}

void Griewank_Function(double *x, double *f, int DIM){

}

void Schwefel_226_Function(double *x, double *f, int DIM){

}

void SUM_SQUARES_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += i * x[i] * x[i];
    }
    f[0] = t1;
}

void LEVY03_Function(double *x, double *f, int DIM){

}

void LEVY05_Function(double *x, double *f, int DIM){

}

void STYBLINSKI_TANG_Function(double *x, double *f, int DIM){

}

void Salomon_Function(double *x, double *f, int DIM){

}

void Whitley(double *x, double *f, int DIM){

}

void Brown_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM - 1; ++i){
        t1 += (pow(x[i] * x[i], x[i + 1] * x[i + 1] + 1) + pow(x[i + 1] * x[i + 1], x[i] * x[i] + 1));
    }
    f[0] = t1;
}

void Exponential_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 0; i < DIM; ++i){
        t1 += x[i] * x[i];
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
        t1 += t2 * t2;
    }
    f[0] = t1;
}

void DIXON_PRICE_Function(double *x, double *f, int DIM){
    double t1 = 0.0;
    for(int i = 1; i < DIM; ++i){
        t1 += i * (2 * x[i] * x[i] - x[i - 1]) * (2 * x[i] * x[i] - x[i - 1]);
    }
    f[0] = (x[0] - 1) * (x[0] - 1) + t1;
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

