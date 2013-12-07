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
@file test_rpim_3.cc

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE) 

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include <algorithm>
#include <map>
#include <cassert>
#include "maxwell/tmz.h"

using namespace lane_maxwell::modules;
using namespace lane_maxwell::modules::maxwell;
using namespace lane_maxwell::modules::matrix;
using namespace lane_maxwell::modules::rpim;
using namespace lane_maxwell::modules::point;
using namespace lane_maxwell::modules::domains;

Point* point2D(double x, double y, double eps, double mu,
               int prop, int ID, char type,
               double sigma_x, double sigma_y) {
    Point* po = new Point;
    Vector* p = new Vector(2);

    p->set(0, x);
    p->set(1, y);

    po->set(ID, eps, mu, prop, type, p);

    po->setSigma(0, sigma_x);
    po->setSigma(1, sigma_y);

    return po;
}

void genPoints(Points* const vp) {
    int N;

    std::FILE* f = std::fopen("rpim.txt", "r");
    fscanf(f, "%d\n", &N);

    double sigma_x, sigma_y;
    std::vector<std::vector<int> > idsd(N);
    std::vector<std::vector<int> > idsdm(N);
    std::vector<std::vector<Point*> > idsdp(N);
    std::map<int, Point*> mvp;

    printf("N = %d\n", N);
    for (int c = 0; c < N; ++c) {
        char type;
        int ID;
        double x, y, eps, mu;
        int prop;
        int n_DS;

        fscanf(f, "%c%d%lf%lf%lf%lf%lf%lf%d%d", &type, &ID, &x, &y, &eps, &mu, &sigma_x, &sigma_y, &prop, &n_DS);
        /*if(ID==1624)*/ printf("nDS: %d %d/%d\n", n_DS, ID, N);


        mvp[ID] = point2D(x, y, eps, mu, prop, ID, type, sigma_x, sigma_y);
        vp->push_back(mvp[ID]);

        for (int k = 0; k < n_DS; ++k) {
            int id;
            fscanf(f, "%d", &id);
            idsd[c].push_back(id);

        }
        fscanf(f,"%d", &n_DS);
        for (int k = 0; k < n_DS; ++k) {
            int id;
            fscanf(f, "%d%lf%lf", &id, &x, &y);
            idsdm[c].push_back(id);
            idsdp[c].push_back(point2D(x, y, 0, 0, 0, 0, 0, 0, 0));
        }
        fscanf(f,"\n");

    }
    printf("ok\n");
    for (int c = 0; c < N; ++c) {
        double rmax = 0;
        const Point* p = vp->at(c);
        for (int k = 0; k < idsd[c].size(); ++k) {
            const Point* q = mvp[idsd[c][k]];
            if (q == NULL) printf("faoip %d %d\n", idsd[c][k], vp->at(c)->getID());
            vp->at(c)->addSD(q);
            double d = distE(p, q);
            if(d > rmax) rmax = d;
        }

        for (int k = 0; k < idsdm[c].size(); ++k) {
            Point* q = mvp[idsdm[c][k]];
            Point* r = idsdp[c][k];
            vp->at(c)->addMirror(q, r);
        }

        vp->at(c)->setRMax(rmax);
    }
    mvp.clear();
    fclose(f);
    printf("ok\n");
}

int main() {
    lane_maxwell::modules::matrix::matrixModuleInit(1e-12);

    std::srand(time(NULL));

    int T, N;

    double C = 299792458.0e0;
    double fmax;

    printf("fmax=");
    scanf("%lf", &fmax);

    double lmb = C / fmax, fat;

    //printf("Numero de pontos=");
    //scanf("%d", &N);

    printf("Lambda / fat = ");
    scanf("%lf", &fat);

    point::Points* p = new Points;

    printf("d = %e\n", lmb/fat);
    genPoints(p);

    T = 12;

    double betamin, betamax, errbeta;
    double nd, delta = lmb / fat;
    double raio;
    int NP;

    printf("betamin betamax errbeta NP = ");
    scanf("%lf%lf%lf%d", &betamin, &betamax, &errbeta, &NP);

    printf("raio=");
    scanf("%lf", &raio);

    point::Points::iterator it = p->begin(), end = p->end();


    int mT = 0;
    std::map<int, rpim::CallDF*> mcdf;

    for (int i = 0; i < p->size(); ++i) {
        T = (*p)[i]->sizeSD();
        if (mcdf.find(T) == mcdf.end()) {
            mcdf[T] = new rpim::CallDF(T, 2);
        }
        if (T > mT) mT = T;
    }

    std::vector<double> F(mT);
    std::vector<double> beta(2);

    bool status;
    int i = 0;

    static char fnout[1000];
    static char fnout2[1000];

    std::FILE* ferr = std::fopen("err", "w+");
    std::FILE* fserr = std::fopen("serr.txt", "w+");
    std::FILE* fsd = std::fopen("sd.txt", "w+");
    std::FILE* fsdc = std::fopen("sdc.txt", "w+");

    std::FILE* fb1_E = std::fopen("fbeta1_E.txt", "w+");
    std::FILE* fb2_E = std::fopen("fbeta2_E.txt", "w+");

    std::FILE* fb1_H = std::fopen("fbeta1_H.txt", "w+");
    std::FILE* fb2_H = std::fopen("fbeta2_H.txt", "w+");

    std::FILE* fpos = std::fopen("positions.txt", "w+");

    int nI5 = 0;
    int nI1 = 0;
    int nI01 = 0;
    int nI001 = 0;
    int kp = 0;
    while (it != end) {
        double err1, err2, err3;
        bool okp = false;

        if ((*it)->sizeMirror() > 0) {
            printf("+++++++++++++++++++\n");
            okp = true;
        }

        for (int type = 02; type < 3; ++type) {
            T = (*it)->sizeSD();

            if (type == 2) {
                beta[0] = rpim::bestBeta(0, &status, i, T, p, mcdf[T], rpim::phi, rpim::dphi, fmax,
                                         betamin, betamax, NP, errbeta, &F, &err1);
                if (status) {
                    beta[1] = rpim::bestBeta(1, &status, i, T, p, mcdf[T], rpim::phi, rpim::dphi, fmax,
                                             betamin, betamax, NP, errbeta, &F, &err2);
                }
            } else if (type == 1) {
                beta[0] = rpim::bestBeta(0, &status, i, T, p, mcdf[T], rpim::phi, rpim::dphi, fmax,
                                         0.1, 0.101, 2, errbeta, &F, &err1);
                if (status) {
                    beta[1] = rpim::bestBeta(1, &status, i, T, p, mcdf[T], rpim::phi, rpim::dphi, fmax,
                                             0.1, 0.101, 2, errbeta, &F, &err2);
                }
            } else if (type == 0) {
                beta[0] = rpim::bestBeta(0, &status, i, T, p, mcdf[T], rpim::phi, rpim::dphi, fmax,
                                         0.3, 0.301, 2, errbeta, &F, &err1);
                if (status) {
                    beta[1] = rpim::bestBeta(1, &status, i, T, p, mcdf[T], rpim::phi, rpim::dphi, fmax,
                                             0.3, 0.301, 2, errbeta, &F, &err2);
                }
            }

            double rmax = (*it)->getRMax();

            int prop = (*it)->getProp();

            double betaF= 0;

            if (prop & maxwell::P_MED) {
                betaF = rpim::bestBetaF(&status, i, T, p, mcdf[T], rpim::phi,
                                        fmax, betamin, betamax, NP, errbeta, &F, &err3);

                rpim::calcF(&status, i, T, p, mcdf[T], betaF*betaF, rpim::phi);

                if (!status) {
                    std::fprintf(stderr, "Nao foi possivel calcula estimativa de campo para o "
                                 "ponto %d\n",  i);
                }
            }

            if (status) {

                rpim::calcDF(0, &status, i, T, p,
                             mcdf[T], beta[0]*beta[0], rpim::phi, rpim::dphi);

                rpim::calcDF(1, &status, i, T, p,
                             mcdf[T], beta[1]*beta[1], rpim::phi, rpim::dphi);

                if (kp % 300 == 0)
                for (size_t k = 0; k < (*it)->sizeSD(); ++k) {
                    const point::Point* pt = (*it)->getSD(k);
                    fprintf(fsd, "%e %e\n", pt->getX(0), pt->getX(1));
                }

                if (kp % 300 == 0) fprintf(fsdc, "%e %e\n", (*it)->getX(0), (*it)->getX(1));

                fprintf(ferr,"%e %e %e\n", err1, err2, err3);


                if(i%400==0 || okp) {
                    if (okp) {

                        printf("%%%\n");
                        for (size_t k = 0; k < (*it)->sizeSD(); ++k) {
                            const point::Point* pt = (*it)->getSD(k);
                            printf("%e %e\n", pt->getX(0), pt->getX(1));
                        }
                        printf("\n MIRROR\n");
                        for (size_t k = 0; k < (*it)->sizeMirror(); ++k) {
                            const point::Point* pt = (*it)->getMirror(k);
                            const point::Point* pt2 = (*it)->getFPoint(k);
                            printf("%e %e (%e %e) %c\n", pt->getX(0), pt->getX(1),
                                   pt2->getX(0), pt2->getX(1), pt->getType());
                        }
                        printf("\np = (%e %e) %c\n", (*it)->getX(0), (*it)->getX(1), (*it)->getType());
                    }
                    printf("<>%d %d rmax = %e "
                       "b1 = %e b2 = %e bf = %e %c\n", T, i, rmax, beta[0]*beta[0], beta[1]*beta[1], betaF, (*it)->getType());
                    printf("%e %e\n", err1, err2);
                }

                if ((*it)->getX(0) > 0.89 &&  (*it)->getX(0) < 1.11 &&
                     (*it)->getX(1) > 0.19 &&  (*it)->getX(1) < 0.41) {

                    if ((*it)->getType() == 'E') {
                        fprintf(fb1_E, "%e %e %e\n",(*it)->getX(0), (*it)->getX(1), beta[0]*beta[0]);
                        fprintf(fb2_E, "%e %e %e\n",(*it)->getX(0), (*it)->getX(1), beta[1]*beta[1]);
                    } else {
                        fprintf(fb1_H, "%e %e %e\n",(*it)->getX(0), (*it)->getX(1), beta[0]*beta[0]);
                        fprintf(fb2_H, "%e %e %e\n",(*it)->getX(0), (*it)->getX(1), beta[1]*beta[1]);
                    }
                }

                if ((*it)->getProp() & P_FONTE) {
                    fprintf(fpos, "%e %e FONTE\n",(*it)->getX(0), (*it)->getX(1));
                }

                if (err1 < 5 && err2 < 5) { nI5++; }
                if (err1 < 1 && err2 < 1)  { nI1++; }
                if (err1 < 0.1 && err2 < 0.1) { nI01++; }
                if (err1 < 0.01 && err2 < 0.01) { nI001++;  }



            }else {
                printf("+%d %d %e %e [%e %e] %e\n",i, T, err1, err2, beta[0]*beta[0], beta[1]*beta[1], rmax);
                err1 = err2 =10;
                if (kp % 300 == 0)
                    for (size_t k = 0; k < (*it)->sizeSD(); ++k) {
                        const point::Point* pt = (*it)->getSD(k);
                        fprintf(fsd, "%e %e\n", pt->getX(0), pt->getX(1));
                    }
                if (kp % 300 == 0)  fprintf(fsdc, "%e %e\n", (*it)->getX(0), (*it)->getX(1));
            }
        }

        if (i%100==0) fprintf(fserr, "%e %e %e\n", (*it)->getX(0), (*it)->getX(1), std::sqrt(err1*err1+err2*err2));

        ++i;
        ++kp;
        ++it;


    }

    printf("err < 5 %% = %d\n", nI5);
    printf("err < 1 %% = %d\n", nI1);
    printf("err < 0.1 %% = %d\n", nI01);
    printf("err < 0.01 %% = %d\n", nI001);

    std::fclose(fserr);
    std::fclose(ferr);
    std::fclose(fsd);
    std::fclose(fsdc);

    std::fclose(fb1_E);
    std::fclose(fb2_E);

    std::fclose(fb1_H);
    std::fclose(fb2_H);

    std::fclose(fpos);

    rpim::save("maxwell.rpim", p);
    free(p);

    for (int i = 0; i <= mT; ++i) if (mcdf.find(i) != mcdf.end()) delete mcdf[i];

    return 0;
}


