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
@file rpim.h


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#ifndef RPIM_RPIM_H_
#define RPIM_RPIM_H_

#include "domains/supportdomain.h"

namespace lane_maxwell { namespace modules { namespace rpim {
double phi(int M, const point::Point* p, const point::Point* c,
            double beta, double raio);
double dphi(int M, const point::Point* p, const point::Point* c,
            double beta, double raio);

typedef double (*RadialFunction)(int, const point::Point*,
            const point::Point*, double, double);

struct CallDF {
    explicit CallDF(size_t ps, size_t D)
        : N(ps), M(D+1) {
        Ro = new matrix::Matrix(N, N);
        Po = new matrix::Matrix(N, M);
        tPo = new matrix::Matrix(M, N);

        aux = new matrix::Matrix(M, N);
        auxPo = new matrix::Matrix(M, M);

        Sb = new matrix::Matrix(M, N);
        iRoPo = new matrix::Matrix(N, M);
        iRoPoSb = new matrix::Matrix(N, N);
        Sa = new matrix::Matrix(N, N);

        ciRo = new matrix::CallInv(N);
        ciauxPo = new matrix::CallInv(M);
    }

    ~CallDF() {
        delete ciRo;
        delete ciauxPo;

        delete Ro;
        delete Po;
        delete tPo;
        delete aux;
        delete auxPo;

        delete Sb;
        delete iRoPo;
        delete iRoPoSb;

        delete Sa;
    }

    matrix::CallInv *ciRo;
    matrix::CallInv *ciauxPo;

    size_t N, M;
    matrix::Matrix *Ro, *Po;
    matrix::Matrix *tPo, *iRo;
    matrix::Matrix *aux;
    matrix::Matrix *auxPo;
    matrix::Matrix *iauxPo;

    matrix::Matrix *Sb;
    matrix::Matrix *iRoPo;
    matrix::Matrix *iRoPoSb;

    matrix::Matrix* Sa;
};

void calcDF(bool* const st, size_t I, size_t ps, point::Points* P,
            CallDF* const cdf, double beta, RadialFunction R, RadialFunction dR,
            point::FPDist dist = point::distE);

void calcDF(int k, bool* const st, size_t I, size_t ps, point::Points* P,
            CallDF* const cdf, double beta, RadialFunction R, RadialFunction dR,
            point::FPDist dist = point::distE);

double bestBeta(int k, bool* const st, size_t I, size_t ps, point::Points* P,
    CallDF* const cdf, RadialFunction R, RadialFunction dR, double fmax,
    double betamin, double betamax, int NP, double errmax, std::vector<double>* F,
    double* err = 0);



void calcF(bool* const st, size_t I, size_t ps, point::Points* P,
           CallDF* const cdf, double beta, RadialFunction R,
           point::FPDist dist = point::distE);

double bestBetaF(bool* const st, size_t I, size_t ps, point::Points* P,
    CallDF* const cdf, RadialFunction R, double fmax,
    double betamin, double betamax, int NP, double errmax, std::vector<double>* F,
    double* err = 0);


bool save(const char* const fname, const point::Points* const p);

/**
 * This function load a points set from file "fname".
 * */
bool load(const char* const fname, point::Points** p);
}}}  // lane_maxwell::modules::rpim

#endif  // RPIM_RPIM_H_
