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
@file test_maxwell_1.cc

Programa utilizado para calcular coeficiente de reflexão numa parede dielétrica infinita.

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include <algorithm>
#include <cmath>

#include "maxwell/tmz.h"
#include "maxwell/maxwellgendom.h"


using namespace lane_maxwell::modules;
using namespace lane_maxwell::modules::maxwell;
using namespace lane_maxwell::modules::matrix;
using namespace lane_maxwell::modules::rpim;
using namespace lane_maxwell::modules::point;
using namespace lane_maxwell::modules::domains;

Point* point2D(double x, double y, double eps = 0, double mu = 0,
               int prop = 0, int ID = 0, char type = 0,
               double sigma_x = 0, double sigma_y = 0) {
    Point* po = new Point;
    Vector* p = new Vector(2);

    p->set(0, x);
    p->set(1, y);

    po->set(ID, eps, mu, prop, type, p);

    po->setSigma(0, sigma_x);
    po->setSigma(1, sigma_y);

    return po;
}

double calcSigma(double x, double dx,
                 double Lx, int *prop) {
    double m = 4.0;
    double sigma_max_x = 15.0;
    double tam_upml_x = 20 * dx;
    double x_upml_1 = 20 * dx;
    double x_upml_2 = Lx - 20 * dx;
    double sigma_x = 0;

    sigma_max_x  = (1 + m) / (15 * M_PI * 15 * dx);
    //sigma_max_x  = 15;

    if (x <= x_upml_1) {
        sigma_x = pow(fabs(x_upml_1 - x) / tam_upml_x, m) * sigma_max_x;
        *prop |= P_UPML;
    } else if (x >= x_upml_2) {
        sigma_x = pow(fabs(x_upml_2 - x) / tam_upml_x, m) * sigma_max_x;
        *prop |= P_UPML;
    }

    return sigma_x;
}


/*
dx = dy = 2e-3
Ts = 4000 * 3e-12
dt = 3e-12

fonte = 0.1
r = 0.05
PMED = (1.06, 0.33) (1.076, 0.33) (1.088, 0.33) (1.104, 0.33) (1.12, 0.33) (1.136, 0.33)
*/

void genPoints(Points* const vp) {

    double eps = 8.854187817620e-12;
    double mu = 4 * acos(-1) * 1e-7;

    std::FILE* fpE = std::fopen("pointsE.dat", "w+");
    std::FILE* fpH = std::fopen("pointsH.dat", "w+");
    std::FILE* fPEC = std::fopen("pointsPEC.dat", "w+");
    std::FILE* fUPML = std::fopen("pointsUPML.dat", "w+");
    std::FILE* fFONTE = std::fopen("pointsFONTE.dat", "w+");
    std::FILE* frpim = std::fopen("rpim.txt", "w+");

    double co = 3e8;
    double fmax = 1e9;
    double lmb = co / fmax;
    double dx = lmb / 40;
    double dy = lmb / 40;
    double Lx = 300 * dx;
    double Ly =  300 * dy;
    int ID = 0;
    int nx = static_cast<int>(Lx / dx);
    int ny = static_cast<int>(Ly / dy);

    nx;
    ny;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            char type = (i % 2 == j % 2)?'E':'H';
            double x = i * dx, y = j * dy;
            double sx = 0, sy = 0;
            int prop = 0;
            sx = calcSigma(x, dx, Lx, &prop);

            ID = i * ny + j;

            if (i < 10 || i > nx - 10 ||
                    (j > nx / 2 + 20  &&
                     j < nx / 2 + 40)) {
                prop |= P_PEC;
            } else {
                if (j == nx / 2 and type == 'E') {
                    prop |= P_FONTE;
                }
            }

            vp->push_back(point2D(x, y, eps, mu, prop, ID, type, sx, sy));

            if (prop & P_PEC) {
                fprintf(fPEC, "%e %e\n", x, y);
            } else  if (prop & P_UPML) {
                fprintf(fUPML, "%e %e\n", x, y);
            } else if (prop & P_FONTE) {
                fprintf(fFONTE, "%e %e\n", x, y);
            } else if (type == 'E') {
                fprintf(fpE, "%e %e\n", x, y);
            } else {
                fprintf(fpH, "%e %e\n", x, y);
            }
        }
    }

    int N = nx * ny;

    fprintf(frpim, "%d\n", N);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            ID = i * ny + j;
            Point* p = vp->at(ID);
            Point* q = 0;
            int prop = 0;
            char type = p->getType();

            fprintf(frpim, "%c %d %e %e %e %e %e %e %d\n", type,
                    p->getID(),
                    p->getX(0), p->getX(1),
                    p->getEps(), p->getMu(),
                    p->getSigma(0),
                    p->getSigma(1),
                    p->getProp());

            ID = i * ny +  ((j-1 < 0)? (ny-1):(j-1));
            q = vp->at(ID);
            p->addSD(q);

            if (j - 1 < 0) {
                p->addMirror(q, point2D(0, -(ny)*dy, 0, 0, 0, 0, 0, 0, 0));
            }

            ID = i * ny +  ((j+1 < ny)? (j+1):0);
            q = vp->at(ID);
            p->addSD(q);


            if (j +1 >= ny) {
                p->addMirror(q, point2D(0, (ny)*dy, 0, 0, 0, 0, 0, 0, 0));
            }

            if (i - 1 >= 0) {
                ID = (i-1) * ny +  j;
                q = vp->at(ID);
                p->addSD(q);
            }

            if (i + 1 < nx) {
                ID = (i+1) * ny +  j;
                q = vp->at(ID);
                p->addSD(q);
            }

            fprintf(frpim, "%d\n", p->sizeSD());
            for (int i = 0 ;i < p->sizeSD(); ++i) {
                const Point* q = p->getSD(i);
                fprintf(frpim, "%d ", q->getID());
            }
            fprintf(frpim, "\n%d\n", p->sizeMirror());
            for (int i = 0; i < p->sizeMirror(); ++i) {
                Point* q = p->getMirror(i);
                Point* r = p->getFPoint(i);
                fprintf(frpim,"%d %e %e ",q->getID(), r->getX(0), r->getX(1));
            }
            fprintf(frpim, "\n");
        }
    }



    std::fclose(fpE);
    std::fclose(fpH);
    std::fclose(fPEC);
    std::fclose(fUPML);
    std::fclose(fFONTE);
    std::fclose(frpim);
}

int main(int argc, char** argv) {
    lane_maxwell::modules::matrix::matrixModuleInit(1e-7);

    point::Points* p = new Points;

    genPoints(p);

    rpim::save("maxwell.rpim", p);

    free(p);
    return 0;
}


