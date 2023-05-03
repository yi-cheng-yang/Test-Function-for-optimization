/*
* Test function for NSYSU Swarm intelligence
* Editor: yi cheng yang
*/


#include <math.h>
#include <float.h>
#include <stdio.h>

double cal_test_function(const double *x, const int d, const int func_num);
void set_search_range(double *range_max, double *range_min, const int func_num);

double Ackley(const double *x, const int d);
double Griewank(const double *x, const int d);
double BentCigar(const double *x, const int d);
double Michalewicz(const double *x, const int d);
double Rosenbrock(const double *x, const int d);
double Schwefel(const double *x, const int d);
double Zakharov(const double *x, const int d);
double HappyCat(const double *x, const int d);
double Rastrigin(const double *x, const int d);
double HGBat(const double *x, const int d);
/*appendix*/
double SchafferF7(const double *x, const int d);
double ExpandedSchafferF6(const double *x, const int d);
double HighConditionedElliptic(const double *x, const int d);
double Weierstrass(const double *x, const int d);
double DropWave(const double *x, const int d);
double Quintic(const double *x, const int d);
double Salomon(const double *x, const int d);
double DixonPrice(const double *x, const int d);
double Katsuura(const double *x, const int d);
double ExpandedGriewankPlusRosenbrock(const double *x, const int d);
double Pinter(const double *x, const int d);
double StyblinskiTang(const double *x, const int d);
double Step(const double *x, const int d);
double Sphere(const double *x, const int d);
double Whitley(const double *x, const int d);
double Discus(const double *x, const int d);
/*-------*/

void set_search_range(double *range_max, double *range_min, const int func_num){
    switch (func_num){
        case 1:
            *range_max =  32.768;
            *range_min = -32.768;
            printf("Initialize Ackley function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 2:
            *range_max =  600.0;
            *range_min = -600.0;
            printf("Initialize Griewank function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 3:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize BentCigar function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 4:
            *range_max =  M_PI;
            *range_min =  0.0;
            printf("Initialize Michalewicz function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 5:
            *range_max =  10.0;
            *range_min = -10.0;
            printf("Initialize Rosenbrock function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 6:
            *range_max =  500.0;
            *range_min = -500.0;
            printf("Initialize Schwefel function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 7:
            *range_max =  10.0;
            *range_min = -10.0;
            printf("Initialize Zakharov function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 8:
            *range_max =  20.0;
            *range_min = -20.0;
            printf("Initialize HappyCat function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 9:
            *range_max =  5.12;
            *range_min = -5.12;
            printf("Initialize Rastrigin function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 10:
            *range_max =  15.0;
            *range_min = -15.0;
            printf("Initialize HGBat function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 11:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize SchafferF7 function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 12:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize ExpandedSchafferF6 function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 13:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize HighConditionedElliptic function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 14:
            *range_max =  0.5;
            *range_min = -0.5;
            printf("Initialize Weierstrass function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 15:
            *range_max =  5.12;
            *range_min = -5.12;
            printf("Initialize DropWave function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 16:
            *range_max =  20.0;
            *range_min = -20.0;
            printf("Initialize Quintic function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 17:
            *range_max =  20.0;
            *range_min = -20.0;
            printf("Initialize Salomon function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 18:
            *range_max =  10.0;
            *range_min = -10.0;
            printf("Initialize DixonPrice function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 19:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize Katsuura function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 20:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize Expanded Griewank plus Rosenbrock function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 21:
            *range_max =  10.0;
            *range_min = -10.0;
            printf("Initialize Pinter function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 22:
            *range_max =  5.0;
            *range_min = -5.0;
            printf("Initialize StyblinskiTang function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 23:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize Step function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 24:
            *range_max =  1.0;
            *range_min = -1.0;
            printf("Initialize Sphere function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 25:
            *range_max =  10.24;
            *range_min = -10.24;
            printf("Initialize Whitley function with range [%f, %f]\n", *range_min, *range_max);
            break;
        case 26:
            *range_max =  100.0;
            *range_min = -100.0;
            printf("Initialize Discus function with range [%f, %f]\n", *range_min, *range_max);
            break;
        default:
            printf("invalid function.\n");
            abort();
            break;
    }
}
double cal_test_function(const double *x, const int d, const int func_num){
    double f = DBL_MAX;
    switch (func_num){
        case 1:
            f = Ackley(x, d);
            break;
        case 2:
            f = Griewank(x, d);
            break;
        case 3:
            f = BentCigar(x, d);
            break;
        case 4:
            f = Michalewicz(x, d);
            break;
        case 5:
            f = Rosenbrock(x, d);
            break;
        case 6:
            f = Schwefel(x, d);
            break;
        case 7:
            f = Zakharov(x, d);
            break;
        case 8:
            f = HappyCat(x, d);
            break;
        case 9:
            f = Rastrigin(x, d);
            break;
        case 10:
            f = HGBat(x, d);
            break;
        case 11:
            f = SchafferF7(x, d);
            break;
        case 12:
            f = ExpandedSchafferF6(x, d);
            break;
        case 13:
            f = HighConditionedElliptic(x, d);
            break;
        case 14:
            f = Weierstrass(x, d);
            break;
        case 15:
            f = DropWave(x, d);
            break;
        case 16:
            f = Quintic(x, d);
            break;
        case 17:
            f = Salomon(x, d);
            break;
        case 18:
            f = DixonPrice(x, d);
            break;
        case 19:
            f = Katsuura(x, d);
            break;
        case 20:
            f = ExpandedGriewankPlusRosenbrock(x, d);
            break;
        case 21:
            f = Pinter(x, d);
            break;
        case 22:
            f = StyblinskiTang(x, d);
            break;
        case 23:
            f = Step(x, d);
            break;
        case 24:
            f = Sphere(x, d);
            break;
        case 25:
            f = Whitley(x, d);
            break;
        case 26:
            f = Discus(x, d);
            break;
        default:
            printf("invalid function.\n");
            break;
    }
    return f;
}

double Ackley(const double *x, const int d){
    double sum1 = 0.0, sum2 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
        sum2 += cos(2.0 * M_PI * x[i]);
    }
    return (-20.0) * exp(-0.2 * sqrt(sum1 / double(d))) - exp(sum2 / double(d)) + 20.0 + exp(1.0);
}
double Griewank(const double *x, const int d){
    double sum1 = 0.0, product1 = 1.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
        product1 *= cos(x[i] / sqrt(i + 1));
    }
    return (sum1 / 4000.0) - product1 + 1.0;
}
double BentCigar(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 1; i < d; ++i){
        sum1 += x[i] * x[i];
    }
    return x[0] * x[0] + pow(10.0, 6) * sum1;
}
double Michalewicz(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += sin(x[i]) * pow(sin((double(i + 1) * x[i] * x[i]) / M_PI), 20.0);
    }
    return sum1 * (-1);
}
double Rosenbrock(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d - 1; ++i){
        sum1 += 100.0 * (x[i + 1] - (x[i] * x[i])) * (x[i + 1] - (x[i] * x[i])) + ((x[i] - 1.0) * (x[i] - 1.0));
    }
    return sum1;
}
double Schwefel(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * sin(sqrt(fabs(x[i])));
    }
    return 418.9829 * double(d) - sum1;
}
double Zakharov(const double *x, const int d){
    double sum1 = 0.0, sum2 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
        sum2 += 0.5 * (i + 1) * x[i];
    }
    return sum1 + pow(sum2, 2) + pow(sum2, 4);
}
double HappyCat(const double *x, const int d){
    double sum1 = 0.0, sum2 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
        sum2 += x[i];
    }
    return pow(fabs(sum1 - double(d)), 0.25) + (0.5 * sum1 + sum2) / double(d) + 0.5;
}
double Rastrigin(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += (x[i] * x[i]) - (10.0 * cos(2.0 * M_PI * x[i]));
    }
    return sum1 + 10.0 * double(d);
}
double HGBat(const double *x, const int d){
    double sum1 = 0.0, sum2 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
        sum2 += x[i];
    }
    return sqrt(fabs((sum1 * sum1) - (sum2 * sum2))) + (0.5 * sum1 + sum2) / double(d) + 0.5;
}

/*appendix*/
double SchafferF7(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d - 1; ++i){
        double v = sqrt((x[i] * x[i]) + (x[i + 1] * x[i + 1]));
        sum1 += sqrt(v) + sqrt(v) * (50.0 * pow(v, 0.2)) * (50.0 * pow(v, 0.2));
    }
    return (1.0 / double(d) * sum1) * (1.0 / double(d) * sum1);
}
double ExpandedSchafferF6(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d - 1; ++i){
        double v = (x[i] * x[i]) + (x[i + 1] * x[i + 1]);
        sum1 += 0.5 + ((sin(sqrt(v)) * sin(sqrt(v))) - 0.5) / ((1.0 + 0.001 * v) * (1.0 + 0.001 * v));
    }
    double v = (x[d - 1] * x[d - 1]) + (x[0] * x[0]);
    return sum1 + 0.5 + ((sin(sqrt(v)) * sin(sqrt(v))) - 0.5) / ((1.0 + 0.001 * v) * (1.0 + 0.001 * v));
}
double HighConditionedElliptic(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        double p = double(i) / double(d - 1);
        sum1 += pow(10.0, 6.0 * p) * x[i] * x[i];
    }
    return sum1;
}
double Weierstrass(const double *x, const int d){
    double sum1 = 0.0, sum2 = 0.0;
    for(int i = 0; i < d; ++i){
        double temp1 = 0.0;
        for(int k = 0; k < 20; ++k){
            temp1 += pow(0.5, k) * cos(2.0 * M_PI * pow(3.0, k) * (x[i] + 0.5));
        }
        sum1 += temp1;
    }
    for(int k = 0; k < 20; ++k){
        sum2 += pow(0.5, k) * cos(M_PI * pow(3.0, k));
    }
    return sum1 - d * sum2;
}
double DropWave(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
    }
    return 1.0 - ((1.0 + cos(12.0 * sqrt(sum1))) / (0.5 * sum1 + 2.0));
}
double Quintic(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += fabs(pow(x[i], 5) - 3.0 * pow(x[i], 4) + 4.0 * pow(x[i], 3) + 2.0 * pow(x[i], 2) - 10.0 * x[i] - 4.0);
    }
    return sum1;
}
double Salomon(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
    }
    return 1.0 - cos(2.0 * M_PI * sqrt(sum1)) + 0.1 * sqrt(sum1);
}
double DixonPrice(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 1; i < d; ++i){
        sum1 += double(i + 1) * (2.0 * x[i] * x[i] - x[i - 1]) * (2.0 * x[i] * x[i] - x[i - 1]);
    }
    return (x[0] - 1) * (x[0] - 1) + sum1;
}
double Katsuura(const double *x, const int d){
    double product1 = 1.0;
    for(int i = 0; i < d; ++i){
        double temp = 0.0;
        for(int j = 0; j < 32; ++j){
            temp += fabs(pow(2.0, j + 1) * x[i] - round(pow(2.0, j + 1) * x[i])) * pow(2.0, -(j + 1));
        }
        product1 *= pow(1.0 + double(i + 1) * temp, 10.0 * pow(d, -1.2));
    }
    return product1 * 10.0 * pow(d, -2) - 10.0 * pow(d, -2);
}
double ExpandedGriewankPlusRosenbrock(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d - 1; ++i){
        double z1 = 0.05 * x[i] + 1.0;
        double z2 = 0.05 * x[i + 1] + 1.0;
        double f1 = 100.0 * (z1 * z1 - z2) * (z1 * z1 - z2) + (z1 - 1.0) * (z1 - 1.0);
        sum1 += (f1 * f1) / 4000.0 - cos(f1) + 1.0;
    }
    double z1 = 0.05 * x[d - 1] + 1.0;
    double z2 = 0.05 * x[0] + 1.0;
    double f1 = 100.0 * (z1 * z1 - z2) * (z1 * z1 - z2) + (z1 - 1.0) * (z1 - 1.0);
    sum1 += (f1 * f1) / 4000.0 - cos(f1) + 1.0;
    return sum1;
}
double Pinter(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        double A;
        double B;
        if(i == 0){
            A = x[d - 1] * sin(x[i] + sin(x[i + 1]));
            B = pow(x[d - 1], 2) - 2.0 * x[i] + 3.0 * x[i + 1] - cos(x[i]) + 1.0;
        }
        else if(i == d - 1){
            A = x[d - 1] * sin(x[i] + sin(x[0]));
            B = pow(x[i - 1], 2) - 2.0 * x[i] + 3.0 * x[0] - cos(x[i]) + 1.0;
        }
        else{
            A = x[i - 1] * sin(x[i] + sin(x[i + 1]));
            B = pow(x[i - 1], 2) - 2.0 * x[i] + 3.0 * x[i + 1] - cos(x[i]) + 1.0;
        }
        sum1 += (i + 1) * pow(x[i], 2) + 20.0 * (i + 1) * pow(sin(A), 2) + (i + 1) * log10(1.0 + (i + 1) * pow(B, 2));
    }
    return sum1;
}
double StyblinskiTang(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += pow(x[i], 4) - 16.0 * x[i] * x[i] + 5.0 * x[i];
    }
    return 0.5 * sum1;
}
double Step(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += floor(x[i] + 0.5) * floor(x[i] + 0.5);
    }
    return sum1;
}
double Sphere(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        sum1 += x[i] * x[i];
    }
    return sum1;
}
double Whitley(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 0; i < d; ++i){
        double sum2 = 0.0;
        for(int j = 0; j < d; ++j){
            double v = 100.0 * (x[i] * x[i] - x[j]) * (x[i] * x[i] - x[j]) + (1.0 - x[j]) * (1.0 - x[j]);
            double a = (v * v) / 4000.0;
            double b = cos(v);
            sum2 += a - b + 1.0;
        }
        sum1 += sum2;
    }
    return sum1;
}
double Discus(const double *x, const int d){
    double sum1 = 0.0;
    for(int i = 1; i < d; ++i){
        sum1 += x[i] * x[i];
    }
    return pow(10.0, 6) * x[0] * x[0] + sum1;
}
/*-------*/

