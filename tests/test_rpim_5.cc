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

#include "maxwell/maxwellgendom.h"

using namespace lane_maxwell::modules;
using namespace lane_maxwell::modules::maxwell;
using namespace lane_maxwell::modules::matrix;
using namespace lane_maxwell::modules::rpim;
using namespace lane_maxwell::modules::point;
using namespace lane_maxwell::modules::domains;

Point* point2D(double x, double y, double eps, double mu, int prop, int ID, char type) {
    Point* po = new Point;
    Vector* p = new Vector(2);

    p->set(0, x);
    p->set(1, y);

    po->set(ID, eps, mu, prop, type, p);
    return po;
}

void genPoints(Points* const vp, int nx, int ny, char type) {

    double dx = 0.01;
    double dy = 0.01;

    double eps = 8.841941e-12;
    double mu = 1.256637e-6;
    int prop = 0;
    static int ID = 0;

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            double x, y;
            if (type == 'E') {
                x = (2*ix)*dx;
                y = (2*iy)*dy;
                vp->push_back(point2D(x, y, eps, mu, prop, ID, type));
                ++ID;

                x = (2*ix+1)*dx;
                y = (2*iy+1)*dy;
                vp->push_back(point2D(x, y, eps, mu, prop, ID, type));
                ++ID;

            } else {
                x = (2*ix+1)*dx;
                y = (2*iy)*dy;
                vp->push_back(point2D(x, y, eps, mu, prop, ID, type));
                ++ID;

                x = (2*ix)*dx;
                y = (2*iy+1)*dy;
                vp->push_back(point2D(x, y, eps, mu, prop, ID, type));
                ++ID;
            }
        }
    }
}

double dumb(const Points* p, int T, char A, char B) {
    double acc = 0;

    for (size_t i = 0; i < p->size(); ++i) {
        if(p->at(i)->getType() == A) {
            std::vector<double> d;

            for (size_t j = 0; j < p->size(); ++j)
                if (i != j && p->at(j)->getType() == B) {
                    d.push_back(distE(p->at(i), p->at(j)));
                }

            sort(d.begin(), d.end());
            for (int j = 0; j < T; ++j) {
                acc+=d[j];
            }
        }
        if (i%100 == 0) std::printf("acc = %e\n", acc);
    }

    return acc;
}


double f(Point* p) {
    double x = p->getX(0);
    double y = p->getX(1);

    return 100*x*x*x*y*y+x*x+y*y+x*y+3+5*x+2*y;
}

double dfx(Point* p) {
    double x = p->getX(0);
    double y = p->getX(1);

    return 300*x*x*y*y+2*x+y+5;
}

double dfy(Point* p) {
    double x = p->getX(0);
    double y = p->getX(1);

    return 200*x*x*x*y+2*y+x+2;
}



int main() {
    lane_maxwell::modules::matrix::matrixModuleInit(1e-7);

    point::Points* p = new Points;
    int nx, ny, T;
    double beta;
    std::printf("(nx,ny)=");
    std::scanf("%d%d", &nx, &ny);

    std::printf("T=");
    std::scanf("%d", &T);

    std::printf("beta");
    std::scanf("%lf", &beta);

    genPoints(p, nx, ny, 'E');
    genPoints(p, nx, ny, 'H');

    for (size_t i = 0; i < p->size(); i++) {
        Point* pt = p->at(i);
        pt->resizeE(1);
        pt->setE(0, f(pt));
    }

    std::printf("N threads: ");
    int NTHREADS;
    scanf("%d", &NTHREADS);

    bool status = calcDF(p, T, NTHREADS, distE, beta, false, 0, phi, dphi, false, 0, 0, 0,
    0, 0);

    double err = 0;

    if (!status) {
        std::printf("Fail calcDF!\n");
    } else {
        for (size_t i = 0; i < p->size(); i++) {
            Point* pt = p->at(i);
            double dfx_rpim = pt->calcDFE(0, 0);
            double dfx_real = dfx(pt);

            //std::printf("dfx_rpim = %e\ndfx_real = %e\nerro = %e\n\n", dfx_rpim,
            //dfx_real, 100*std::fabs(dfx_rpim-dfx_real)/fabs(dfx_real));

            err += 100*std::fabs(dfx_rpim-dfx_real)/fabs(dfx_real);
        }

        for (size_t i = 0; i < p->size(); i++) {
            Point* pt = p->at(i);
            double dfy_rpim = pt->calcDFE(0, 1);
            double dfy_real = dfy(pt);

            //std::printf("dfy_rpim = %e\ndfy_real = %e\nerro = %e\n\n", dfy_rpim,
            //	dfy_real, 100*std::fabs(dfy_rpim-dfy_real)/fabs(dfy_real));

            err += 100*std::fabs(dfy_rpim-dfy_real)/fabs(dfy_real);
        }

        std::printf("Err medium = %g %g\n", err, err / (2*p->size()));

        rpim::save("test_rpim.rpim", p);
    }

    free(p);

    pthread_exit(NULL);
    return 0;
}


