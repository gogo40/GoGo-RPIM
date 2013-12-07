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
#include <utility>
#include <map>
#include "maxwell/maxwellgendom.h"
#include "lane/lane.h"

namespace lane_maxwell { namespace modules { namespace maxwell {

typedef std::map<int, rpim::CallDF*> MapCDF;
typedef MapCDF::iterator CDFIt;
typedef std::pair<int, rpim::CallDF*> PairCDF;

struct CalcDFArgs {
    size_t ini;
    size_t fim;
    bool status;
    bool isauto;
    std::vector<double>* betax;
    std::vector<double>* Fbetax;
    int NP;
    double fmax;
    size_t T;
    point::Points* P;
    MapCDF* cdf;
    double beta;
    bool useR;
    double rd;
    rpim::RadialFunction R;
    rpim::RadialFunction dR;
    double betamin, betamax, errbeta;

    CalcDFArgs (
        size_t ini,
        size_t fim,
        bool status,
        bool isauto,
        std::vector<double>* betax,
        std::vector<double>* Fbetax,
        int NP,
        double fmax,
        size_t T,
        point::Points* P,
        MapCDF* cdf,
        double beta,
        bool useR, double rd,
        rpim::RadialFunction R,
        rpim::RadialFunction dR,
        double betamin, double betamax,
        double errbeta
    )
    :
    ini(ini), fim(fim), status(status),
    isauto(isauto), betax(betax), Fbetax(Fbetax),
    NP(NP), fmax(fmax),
    T(T), P(P), cdf(cdf), beta(beta),
    useR(useR), rd(rd), R(R), dR(dR), betamin(betamin),
      betamax(betamax),errbeta(errbeta){}

    ~CalcDFArgs() {
        if(isauto) {
            delete betax;
            delete Fbetax;
        }
    }

    DISALLOW_COPY_AND_ASSIGN(CalcDFArgs);
};

void* threadCalcDF(void* p_args) {
    CalcDFArgs* a = static_cast<CalcDFArgs*>(p_args);
    size_t ini = a->ini;
    size_t fim = a->fim;
    size_t T = a->T;
    int prop;
    bool useR = a->useR;
    point::Points* P = a->P;
    MapCDF* cdf = a->cdf;
    double beta = a->beta;
    bool isauto = a->isauto;
    double betamin = a->betamin;
    double betamax = a->betamax;
    double errbeta = a->errbeta;
    double err;
    int NP = a->NP;
    double fmax = a->fmax;
    std::vector<double>* F = a->Fbetax;
    std::vector<double>* betax = a->betax;
    rpim::RadialFunction R = a->R;
    rpim::RadialFunction dR = a->dR;
    double betaF;

    for (size_t i = ini; i < fim; ++i) {

        T = P->at(i)->sizeSD();
        prop = P->at(i)->getProp();
        //if (useR) {
            if (prop & P_PEC) continue;
        //}

        if (isauto) {
            for (size_t j = 0; j < P->at(i)->getDim(); ++j) {

                (*betax)[j] = rpim::bestBeta(j, &(a->status), i, T, P, (*cdf)[T], R, dR, fmax,
                betamin, betamax, NP, errbeta, F, &err);


                rpim::calcDF(j, &(a->status), i, T, P, (*cdf)[T],
                (*betax)[j] * (*betax)[j], R, dR);

                if (!a->status || err > 0.1) {
                    std::fprintf(stderr, "Nao foi possivel calcular %d - esima"
                    " derivada do ponto %d\nerr = %e\n", j, i,err);
                    //break;
                }
            }

            if (prop & maxwell::P_MED) {
                betaF = rpim::bestBetaF(&(a->status), i, T, P, (*cdf)[T], R,
                                        fmax, 0.01, 0.02, 2, errbeta, F, &err);

                if (err > 0.01)  {
                    betaF = rpim::bestBetaF(&(a->status), i, T, P, (*cdf)[T], R,
                                            fmax, 0.1, 0.2, 2, errbeta, F, &err);
                    if (err > 0.01) {
                        betaF = rpim::bestBetaF(&(a->status), i, T, P, (*cdf)[T], R,
                                                fmax, betamin, betamax, NP, errbeta, F, &err);
                    }

                }

                rpim::calcF(&(a->status), i, T, P, (*cdf)[T], betaF * betaF, R);

                if (!a->status || err > 0.1) {
                    std::fprintf(stderr, "Nao foi possivel calcula estimativa de campo para o "
                    "ponto %d\n err = %e\n",  i, err);
                    //break;
                }
            }

            if (a->ini == 0 && i % 1000 == 0) {
                for (size_t j = 0; j < P->at(i)->getDim(); ++j) {
                    std::fprintf(stderr, "beta %d = %e [%e]\n", j+1,
                                 (*betax)[j],
                                 err);

                }
                if (prop & maxwell::P_MED || err > 0.1) {
                    std::fprintf(stderr, "beta F = %e [%e]\n",
                                 betaF,
                                 err);
                }
            }
        } else {
            rpim::calcDF(&(a->status), i, T, P, (*cdf)[T], beta,  R, dR);
            if (!a->status) {
                //break;
            }

            if (prop & maxwell::P_MED) {
                rpim::calcF(&(a->status), i, T, P, (*cdf)[T], beta, R);
            }
        }

        if ( i % 1000 == 0) {
            if (a->ini == 0)
                fprintf(stderr, "df %d / %d\n", i+1, fim);
        }
    }
    //pthread_exit(NULL);
    return NULL;
}

inline double min(double a, double b) {
    return (a < b)? a : b;
}

bool calcDF(
       point::Points* pts, size_t T, size_t NTHREADS,
        domains::FPDist D, double beta,
        bool useR, double rd,
        rpim::RadialFunction R, rpim::RadialFunction dR,
        bool isauto, double betamin,
        double betamax, double errbeta, int NP, double fmax) {

    std::sort(pts->begin(), pts->end(), point::PointComp());
    double rmin;
    rmin = domains::buildSupportDomain(stdout, 'E', 'H', T, pts, 0.01, D, 1, useR, rd);
    rmin = min(rmin, domains::buildSupportDomain(stdout, 'H', 'E', T, pts, 0.01, D, 1, useR, rd));

    size_t dim = pts->at(0)->getDim();
    size_t N = pts->size();
    std::vector<MapCDF*> vcdf(NTHREADS);
    size_t mT = 0;

    for (size_t i = 0; i < NTHREADS; ++i) {
        vcdf[i] = new MapCDF;
    }

    for (size_t i = 0; i < pts->size(); ++i) {
        size_t T = pts->at(i)->sizeSD();
        for (size_t k = 0; k < NTHREADS; ++k) {
            if ( vcdf[k]->find(T) == vcdf[k]->end()) {
                vcdf[k]->insert(PairCDF(T, new rpim::CallDF(T, dim)));
            }
        }
        if (T > mT) {
            mT = T;
        }
    }


    std::vector<CalcDFArgs*> args(NTHREADS);
    std::vector<pthread_t> threads(NTHREADS);

    size_t sT = N / NTHREADS;
    int rc;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

    for (size_t i = 0; i < NTHREADS; ++i) {

        size_t ini = i * sT, fim = (i + 1) * sT;
        fim = min(N, fim);

        args[i] = new CalcDFArgs(ini, fim, true, isauto,
        (isauto)?new std::vector<double>(pts->at(0)->getDim()):0 ,
        (isauto)?new std::vector<double>(mT):0,
        NP, fmax, T,
        pts, vcdf[i], beta, useR, rd, R, dR,
        betamin, betamax, errbeta);
    }

    for (size_t i=0; i < NTHREADS; ++i) {
        pthread_create(&threads[i], &attr, threadCalcDF,
        static_cast<void*>(args[i]));
    }

    pthread_attr_destroy(&attr);
    void* status;

    for(size_t i = 0; i < NTHREADS; ++i){
        rc = pthread_join(threads[i], &status);
        if (rc) {
            fprintf(stderr, "ERROR; return code from"
            "pthread_join() is %d\n", rc);
            return false;
        }
    }

    bool ok = true;

    for (size_t i = 0; i < NTHREADS; ++i) {
        if (!args[i]->status) {
            ok = false;
        }

        for (CDFIt it = vcdf[i]->begin(); it != vcdf[i]->end(); ++it) {
            delete it->second;
        }

        delete vcdf[i];
        delete args[i];
    }

    args.clear();
    vcdf.clear();
    threads.clear();
    printf("rmin = %e\n", rmin);

    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;
    double dt = 0.4 * rmin * sqrt(eps0 * mu0);

    printf("dt = %e\n", dt);

    return ok;
}
}}}  // lane_maxwell::modules::maxwell

