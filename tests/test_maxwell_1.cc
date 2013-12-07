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

double calcSigma(double x, double dx,
                 double Lx, int *prop) {
    double m = 4.0;
    double sigma_max_x = 15.0;
    double tam_upml_x = 10 * dx;
    double x_upml_1 = 10 * dx;
    double x_upml_2 = Lx - 10 * dx;
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

void genPoints(Points* const vp, double disc, double er, double ang) {
    double dx = 2e-3;
    double dy = 2e-3;
    
    double Lx = 640 * dx;
    double Ly = 320 * dy;

    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;
    double eps = 8.854187817620e-12;
    double mu = 4 * acos(-1) * 1e-7;

    int prop = 0;

    double sigma_x;
    double sigma_y;

    int ID = 0;

    int sx = 0;
    int sy = 0;

    int nx = static_cast<int>(Lx / (0.5 * dx));
    int ny = static_cast<int>(Ly / (0.5 * dy));

    printf("dt = %e <|> \n", 0.4 * dx * sqrt(eps * mu));
    printf("nx = %d ny = %d\n", nx, ny);

    for (int ix = -10; ix < nx + 10; ++ix) {
        for (int iy = -10; iy < ny + 10; ++iy) {
            double x, y;
            char type;

            sx = ix%2;
            sy = iy%2;

            if (sx == sy) {
                type = 'E';
            } else {
                type = 'H';
            }

            prop = 0;

            eps = eps0;

            if (ix < 0 || ix >= nx){ prop |= P_METAL; prop |= P_PEC;}
            if (iy < 0 || iy >= ny){ prop |= P_METAL; prop |= P_PEC;}

            x = ix * dx * 0.5;
            y = iy * dy * 0.5;

            if (type == 'E') {

                int isUPML = 0;
                int isSource = 0;

                sigma_x = calcSigma(x, dx, Lx, &isUPML);
                sigma_y = calcSigma(y, dy, Ly, &isUPML);

                double xo = (20+1000)*0.5*dx;
                double yo = (20+300)*0.5*dx;

                double ro = 0.05;
                
                if (ix == 20) {
                    isSource = P_FONTE;
                }

                double dx = x - xo;
                double dy = y - yo;

                double r = sqrt(dx * dx + dy * dy);

                if (r < ro + 1e-5) {
                    prop |= P_METAL; prop |= P_PEC;
                }

                if (iy == 320) {
                    if (ix == 1080) {
                        printf("m (%e %e)\n", x, y);
                        prop |= P_MED;
                    }

                    if (ix == 1096) {
                        printf("m (%e %e)\n", x, y);
                        prop |= P_MED;
                    }

                    if (ix == 1108) {
                        printf("m (%e %e)\n", x, y);
                        prop |= P_MED;
                    }

                    if (ix == 1124) {
                        printf("m (%e %e)\n", x, y);
                        prop |= P_MED;
                    }

                    if (ix == 1140) {
                        printf("m (%e %e)\n", x, y);
                        prop |= P_MED;
                    }

                    if (ix == 1156) {
                        printf("m (%e %e)\n", x, y);
                        prop |= P_MED;
                    }
                }


                vp->push_back(point2D(x, y, eps, mu, prop | isUPML | isSource,
                                      ID, type, sigma_x, sigma_y));
                ++ID;

            } else {
                int isUPML = 0;

                sigma_x = calcSigma(x, dx, Lx, &isUPML);
                sigma_y = calcSigma(y, dy, Ly, &isUPML);

                vp->push_back(point2D(x, y, eps, mu, prop|isUPML,
                                      ID, type, sigma_x, sigma_y));
                ++ID;
            }
        }
    }

}

int main(int argc, char** argv) {
    if (argc != 4) {
        printf("%s <disc> <er> <ang>\n",argv[0]);
        return 0;
    }
    lane_maxwell::modules::matrix::matrixModuleInit(1e-7);

    point::Points* p = new Points;

    double disc, er, ang;

    sscanf(argv[1],"%lf", &disc);
    sscanf(argv[2],"%lf", &er);
    sscanf(argv[3], "%lf", &ang);

    genPoints(p,disc,er, ang);

    rpim::save("maxwell_1.rpim", p);

    free(p);
    return 0;
}


