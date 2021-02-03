#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <ctime>
#define N 30
#define SIZE 20
#define AVG_ITER 30
#define FES 200000

using namespace std;

struct Particle
{
    int index;
    double pBest[N];
    double vel[N];
    double pos[N];
    double fitness;
    double distance;
    double pBest_fitness;    
};

int gBest_index;

double C1;
double C2;
double W;
double FE;


double fitness_val_list[FES];
Particle particle[SIZE];
const double pi = atan(1)*4;
int fes;
double randval(double low, double high);
double fit_fun(double pos[], int fit_func);
void initialize(double pos_max, double vel_max, int fit_func, double c1, double c2, double w);
void update_vel(int i, double vel_max);
void update_pos(int i, double pos_max, int fit_func);
void update(double pos_max, double vel_max, int fit_func);
double GPSO(double pos_max, double vel_max, int fit_func, double c1, double c2, double w);
double calstd(double results[], double avg_result, int n);
int main();

double randval(double low, double high){
    
    return (low+(high-low)*rand()*1.0/RAND_MAX);}

double calstd(double results[], double avg_result, int n){
    double d = 0;
    for (int i = 0; i < n; i++){
        d += pow(results[i]-avg_result,2);
    }
    return sqrt(d/n);
}

double fit_fun(double pos[], int fit_func){
    fes += 1;
    if (fit_func == 1){
        double val = 0;
        for (int i =0; i < N; i++){
            val += pow(pos[i],2);
        }
        return val;
    }    
    
    else if (fit_func == 2){
        double val1 = 0;
        double val2 = 1;
        for (int i = 0; i < N; i++){
            val1 += fabs(pos[i]);
            val2 = val2 * fabs(pos[i]);
        }
        return val1+val2;
    }

    else if (fit_func == 3){
        double val1 = 0;
        double val = 0;
        for (int i = 0; i < N; i++){
            val1 += pos[i];
            val += val1*val1;
        }
        return val;
    }

    else if (fit_func == 4){
        double val = 0;
        for (int i = 0; i < N-1; i++){
            val += 100*(pow(pos[i+1]-pow(pos[i],2),2) + pow(pos[i]-1,2));
        }
        return val;
    }

    else if (fit_func == 5){
        double val = 0;
        for (int i = 0; i < N; i++){
            val += pow(floor(pos[i]+0.5),2);
        }
        return val;
    }

    else if (fit_func == 6){
        double val = 0;
        for (int i = 0; i < N; i++){
            val += (i+1)*pow(pos[i],4);
        }
        return val+rand()*(1.0/RAND_MAX);
    }

    else if (fit_func == 7){
        double val = 0;
        for (int i = 0; i < N; i++){
            val += (-1)*pos[i]*sin(sqrt(fabs(pos[i])));
        }
        return val;
    }

    else if (fit_func == 8){
        double val = 0;
        for (int i = 0; i < N; i++){
            val += (pow(pos[i],2) - 10 * cos(2.0*pi*pos[i]) + 10);
        }
        return val;
    }

    else if (fit_func == 9){
        double val = 0;
        double y;
        for (int i = 0; i < N; i++){
            if (fabs(pos[i]) < 0.5) y = pos[i];
            else y = round(2*pos[i])/2;
            val += (pow(y,2) - 10 * cos(2.0*pi*y) + 10);
        }
        return val;
    }

    else if (fit_func == 10){
        double val1 = 0;
        double val2 = 0;
        for (int i = 0; i < N; i++){
            val1 += pow(pos[i],2);
            val2 += cos(2*pi*pos[i]);
        }
        double val = -20 * exp(-0.2 * sqrt(val1/N)) - exp(val2/N) + 20 + exp(1);
        return val;
    }

    else if (fit_func == 11){
        double val1 = 0;
        double val2 = 1;
        for (int i = 0; i < N; i++){
            val1 += pow(pos[i],2);
            val2 = val2 * cos(pos[i]/sqrt(i+1));
        }
        double val = (1.0/4000)*val1 - val2 + 1;
        return val;
    }

    else if (fit_func == 12){
        double val1 = 0;
        double val2 = 0;
        double u, y, y_2;
        for (int i = 0; i < N-1; i++){
            if (pos[i] > 10) u = 100*pow(pos[i] - 10, 4);
            else if (pos[i] >= -10) u = 0;
            else u = 100*pow((-1) * pos[i] - 10, 4);
            y = 1+0.25*(pos[i]+1);
            y_2 = 1+0.25*(pos[i+1]+1);
            val1 += pow(y-1, 2)*(1 + 10*pow(sin(pi*y_2),2));
            val2 += u;
        }
        if (pos[N-1] > 10) u = 100*pow(pos[N-1] - 10, 4);
        else if (pos[N-1] >= -10) u = 0;
        else u = 100*pow((-1) * pos[N-1] - 10, 4);
        val2 += u;
        y = 1 + 0.25*(pos[0]+1);
        double val = (pi/N)*(10*pow(sin(pi*y),2) + val1 + pow(y_2-1,2)) + val2;
        return val;
    }
}

void initialize(double pos_max, double vel_max, int fit_func, double c1, double c2, double w){
    C1 = c1;
    C2 = c2;
    W = w;
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < N; j++){
            particle[i].pos[j] = randval(-1.0*pos_max, pos_max);
            particle[i].vel[j] = randval(-1.0*vel_max, vel_max);
            particle[i].pBest[j] = particle[i].pos[j];
        }
        particle[i].fitness = fit_fun(particle[i].pos, fit_func);
        particle[i].pBest_fitness = particle[i].fitness;
        if (i == 0 || particle[i].fitness < particle[gBest_index].pBest_fitness) gBest_index = i;
    }
}

void update_vel(int i, double vel_max){
    for (int j = 0; j < N; j++){
        double vel_value;
        vel_value = W * particle[i].vel[j] + C1 * rand()*(1.0/RAND_MAX) * (particle[i].pBest[j] \
                    - particle[i].pos[j]) + C2 * rand()*(1.0/RAND_MAX) * (particle[gBest_index].pBest[j] - particle[i].pos[j]);
        if (vel_value > vel_max) vel_value = vel_max;
        else if (vel_value < -1.0*(vel_max)) vel_value = -1.0*vel_max;
        particle[i].vel[j] = vel_value;
    }
}

void update_pos(int i, double pos_max, int fit_func){
    for (int j = 0; j < N; j++){
        double pos_value = particle[i].pos[j] + particle[i].vel[j];
        if (pos_value > pos_max) pos_value = pos_max;
        else if (pos_value < -1.0*(pos_max)) pos_value = -1.0*pos_max;
        particle[i].pos[j] = pos_value;
    }
    double value = fit_fun(particle[i].pos, fit_func);
    if (value < particle[i].pBest_fitness){
        particle[i].pBest_fitness = value;
        memcpy(particle[i].pBest, particle[i].pos, sizeof(particle[i].pos));
    }
    if (value < particle[gBest_index].pBest_fitness) gBest_index = i;
}

void update_W(){
    W = 0.9 - 0.5*fes/FES;
}

void update(double pos_max, double vel_max, int fit_func){
    
    for (int j = 0; j < SIZE; j++){
        update_vel(j, vel_max);
        update_pos(j, pos_max, fit_func);
        } 
    update_W();
}

double GPSO(double pos_max, double vel_max, int fit_func, double c1, double c2, double w, double acceptfit){
    initialize(pos_max, vel_max, fit_func, c1, c2, w);
    int i = 0;
    int flag = 0;
    while (fes <= FES){
        update(pos_max, vel_max, fit_func);
        fitness_val_list[i] = particle[gBest_index].pBest_fitness;
        i += 1;
        if (flag == 0 && particle[gBest_index].pBest_fitness<=acceptfit){
            flag = 1;
            FE = fes;}
    }
    if (flag == 0){
        FE = FES;
        }
    return fitness_val_list[i-1];
}

double calavg(double results[], int n){
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum += results[i];
    }
    return sum/n;
}

int main()
{
    double pos_max[] = {100, 10, 100, 10, 100, 1.28, 500, 5.12, 5.12, 32, 600, 50};
    double accept[] = {0.01, 0.01, 100, 100, 0, 0.01, -10000, 50, 50, 0.01, 0.01, 0.01};
    double c1 = 2.0;
    double c2 = 2.0;
    double w = 0.9;
    for (int i = 0; i < 3; i++){
        double results[AVG_ITER];
        double FEs[AVG_ITER];
        int fit_func = i + 1;
        for (int j = 0; j < AVG_ITER; j++){
            srand((unsigned)j);
            fes = 0;
            double vel_max = 0.2 * pos_max[i] * 2;
            results[j] = GPSO(pos_max[i], vel_max, fit_func, c1, c2, w, accept[i]);
            cout<<results[j]<<endl;
            //for (int k = 0; k < N; k++) cout<<particle[gBest_index].pBest[k]<<' ';
            //cout<<endl;
            FEs[j] = FE;
        }
        double avg_results = calavg(results, AVG_ITER);
        double avg_FEs = calavg(FEs, AVG_ITER);
        double std_results = calstd(results, avg_results, AVG_ITER);
        cout << "\nResults for fitness function" << fit_func << ":" << endl;
        cout << "average FEs:" << avg_FEs << endl;
        cout << "average fitness:" << avg_results << endl;
        cout << "fitness std:" << std_results << endl;
    }
    return 0;
}
