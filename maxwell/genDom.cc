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
@file maxwellgendom.cc
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/
#include <algorithm>
#include <cstdio>
#include <cstring>
#include "maxwell/maxwellgendom.h"

int main(int argc, char** argv) {
    if (argc < 15) {
        std::printf("%s <input file name> <output file name>"
        " <number threads> <beta> <domain support size> <error>"
        " <useR> <rd> <isauto> <betamin> <betamax> <errbeta> <NP>\n"
        " <fmax> \n", argv[0]);
        return 0;
    }

    double error;

    std::sscanf(argv[6], "%lf", &error);
    lane_maxwell::modules::matrix::matrixModuleInit(error);

    lane_maxwell::modules::point::Points* p;
    if (!lane_maxwell::modules::point::load(argv[1], &p)) {
        fprintf(stderr, "Failed to load file %s\n", argv[1]);
        return 0;
    }

    int useR;
    double betamin, betamax, errbeta;
    int isauto;
    double rd;
    int NTHREADS;
    int NP;
    double fmax;
    double beta;
    int T;

    std::sscanf(argv[3], "%d", &NTHREADS);

    std::sscanf(argv[4], "%lf", &beta);
    std::sscanf(argv[5], "%d", &T);
    std::sscanf(argv[7], "%d", &useR);
    std::sscanf(argv[8], "%lf", &rd);
    std::sscanf(argv[9], "%d", &isauto);
    std::sscanf(argv[10], "%lf", &betamin);
    std::sscanf(argv[11], "%lf", &betamax);
    std::sscanf(argv[12], "%lf", &errbeta);
    std::sscanf(argv[13], "%d", &NP);
    std::sscanf(argv[14], "%lf", &fmax);


    std::printf("NTHREADS = %d\n"
    "beta = %e\nT = %d\n", NTHREADS, beta, T);

    bool status =
    lane_maxwell::modules::maxwell::calcDF(
    p, T, NTHREADS,
    lane_maxwell::modules::point::distE, beta,
    (useR)?true:false, rd,
    lane_maxwell::modules::rpim::phi,
    lane_maxwell::modules::rpim::dphi, (isauto)?true:false,
    betamin, betamax, errbeta, NP, fmax);

    printf("Finalizando geracao de dominios...\n");
    if (status) {
        if (!lane_maxwell::modules::rpim::save(argv[2], p)) {
            fprintf(stderr, "Failed to save file %s\n", argv[2]);
            return 0;
        }
    } else {
        fprintf(stderr, "Failed to calc RPIM coeficients!\n");
    }

    printf("Finalizando geracao de dominios...\n");
    lane_maxwell::modules::point::free(p);

    return 0;
}

