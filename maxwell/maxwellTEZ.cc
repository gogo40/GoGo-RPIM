/*
GoGoRPIM - A RPIM implementation and 2-D Eletromagnetic simulator
    Copyright (C) 2012  Péricles Lopes Machado (LANE-UFPA)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**Copyright 2012 Péricles Lopes Machado
@file maxwellTEZ.cc
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include "maxwell/tez.h"

#ifndef M_PI
const double M_PI = acos(-1);
#endif

enum SourceType {
    SIN = 0,
    COS = 1,
    GAUSS = 2,
    SIGMOID = 3,
    MONOCICLOGAUSS = 4,
    DESCARGAATM = 5,
    USER = 6,
    MONOCICLOGAUSSBAVELIS = 7
};


int type;
double t_ant=-1;
double pico;
double w;
std::vector<double> fonte;
FILE* fpf;

double q0, W;
double fc;


double source(const double& t, const double& dt,
    const lane_maxwell::modules::point::Point* pts) {

    int i;

    double lamb, C, neper, q;
    double tau, t0;
    double f = w / (2 * M_PI);
    double sinal = 0;

    double pico1=2.0e-3;
    double tau1=400e-12;
    double tau2=2.5e-9;
    double nn=2.0e+0;
    double pico2=6.5e+3, tau12 = 2e-6, tau22=230e-6;

    double eta = exp(  - (tau1/tau2)*powf(  nn*tau2/tau1 ,1.0e0/nn )  );

    double eta2 = exp(  - (tau12/tau22)*powf(  nn*tau22/tau12 , 1.0e0/nn )  );


    double x = pts->getX(0);
    double y = pts->getX(1);

    C = 299792458.0e0;
    neper = 2.718281e0;

    lamb = C / f;

    tau = 50.0e0*dt  , t0 = 130.0e0*dt  ;

    switch (type) {
        case SIN:
            sinal = pico * std::sin(w * t);
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case COS:
            sinal = pico * std::cos(w * t);
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case GAUSS:
            sinal = pico * std::sin(w * t) * std::exp(-(t - t0) * (t - t0) / tau / tau); //gaussiana modulada em seno
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case SIGMOID:
            sinal = pico * std::sin(w * t) / (1 + std::exp(5 - (t / tau))); //sigmoid modulada em seno
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case MONOCICLOGAUSS:
            //tau = 10.0e0*dt, t0 = 75.0e0*dt;
            tau = 40.0e0*dt, t0 = 150.0e0*dt;
            sinal = pico * std::sqrt(2*neper/(tau*tau)) * (t-t0) * std::exp(-std::pow((t-t0),2)/(tau*tau)); //geracao do monociclogaussiano.
            if(t > t_ant) { // Evita repetição do mesmo valos para fontes não pontuais
                t_ant = t;
                fprintf(fpf,"%e %e\n", t_ant, sinal);
            }
            return sinal;
        break;

        case USER:
            i = static_cast<int>(t / dt);
            sinal = pico * fonte[i];
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case DESCARGAATM:
            sinal = 3.5e-4*(((pico1/eta)/(pico1/eta)*(exp(-t/tau2))*(pow((t/tau1),1.0e0/nn)/(pow((t/tau1),nn)+1.0e0)))/0.52231);
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case MONOCICLOGAUSSBAVELIS:
            q = t / dt;
            eta = - pow(q - q0 + x / (C * dt), 2) / (2 * W * W);
            tau1 = 2 * M_PI * fc * dt * (q - q0 + x / (C * dt));
            sinal = pico * std::sin(tau1)* std::exp(eta);


            if(t > t_ant) { // Evita repetição do mesmo valos para fontes não pontuais
                t_ant = t;
                fprintf(fpf,"%e %e\n", t_ant, sinal);
            }
            return sinal;
        break;
    }

    return 0;
}

int main(int argc, char** argv) {
    if (argc < 15) {
        printf("%s <input file name> <dt> <Ts> "
        "<resX> <resY> <videoGen> <nf> <isLog> "
        "<lecont> <NTHREADS> <source type> <filename or freq> <pico> <source mode>"
        "[duracao do pulso (opcional)]?\n", argv[0]);
        printf("\nSources type:\n");
        printf("0 = pico * sin(2 * pi * freq)\n");
        printf("1 = pico * cos(2 * pi * freq)\n");
        printf("2 = pico * sin(2 * pi * freq) * exp(-(t-t0)^2/tau^2),\n");
        printf("where t0 = 50dt and tau = 130dt\n");
        printf("3 = pico * fonte[floor(t/dt)], where fonte is a vector\n");
        printf("described in the file \"filename\"\n");
        printf("\nSource modes:\n");
        printf("0 = Hard\n1 = Soft\n 2 = Type 3\n");
        return 0;
    }

    fpf = fopen("fonte.dat","w+");

    double dt;
    std::sscanf(argv[2], "%lf", &dt);

    double Ts;
    std::sscanf(argv[3], "%lf", &Ts);

    double resX;
    std::sscanf(argv[4], "%lf", &resX);

    double resY;
    std::sscanf(argv[5], "%lf", &resY);

    int videoGen;
    std::sscanf(argv[6], "%d", &videoGen);

    int nf;
    std::sscanf(argv[7], "%d", &nf);

    int isLog;
    std::sscanf(argv[8], "%d", &isLog);

    int lecont;
    std::sscanf(argv[9], "%d", &lecont);

    int NTHREADS;
    std::sscanf(argv[10], "%d", &NTHREADS);

    std::sscanf(argv[11], "%d", &type);

    if (type == USER) {
        std::FILE* file = std::fopen(argv[12], "r");

        if (!file) {
            std::fprintf(stderr, "Failed to load source file %s\n", argv[12]);
            return 0;
        }

        while (!feof(file)) {
            double x;
            fscanf(file, "%lf", &x);
            fonte.push_back(x);
        }
    } else {
        std::sscanf(argv[12], "%lf", &w);
        w = 2 * M_PI * w;
    }

    std::sscanf(argv[13], "%lf", &pico);

    int sourceType;

    std::sscanf(argv[14], "%d", &sourceType);

    double duracao = 0;
    if (argc == 16) {
        std::sscanf(argv[15], "%lf", &duracao);
    }
    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;
    double cte = dt / (2*eps0);

    printf("dt = %e cte = %e eps = %e w = %e pico = %e\n", dt, cte, eps0, w, pico);

    lane_maxwell::modules::point::Points* p;
    if (!lane_maxwell::modules::rpim::load(argv[1], &p)) {
        fprintf(stderr, "Failed to load file %s\n", argv[1]);
        return 0;
    }

    lane_maxwell::modules::maxwell::tez(p, Ts, dt, cte,
    resX, resY, (videoGen)? true: false,
    nf, isLog? true: false, lecont, NTHREADS, source,
    sourceType, duracao);

    lane_maxwell::modules::point::free(p);

    return 0;
}

