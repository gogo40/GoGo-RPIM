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

double calc_area(double* x, double* y) {
    double D = fabs(x[0]*y[1] + x[1]*y[2] + x[2]*y[0] -
                    x[0]*y[2] - x[1]*y[0] - x[2]*y[1]);
    return D * 0.5;
}

void genPoints(Points* const vp) {

    double eps = 8.854187817620e-12;
    double mu = 4 * acos(-1) * 1e-7;

    std::FILE* fpE = std::fopen("pointsE.dat", "w+");
    std::FILE* fpH = std::fopen("pointsH.dat", "w+");
    std::FILE* fPEC = std::fopen("pointsPEC.dat", "w+");
    std::FILE* fUPML = std::fopen("pointsUPML.dat", "w+");
    std::FILE* fFONTE = std::fopen("pointsFONTE.dat", "w+");

    double co = 3e8;
    double fmax = 1e9;
    double lmb = co / fmax;
    double dx = lmb / 40;
    double dy = lmb / 40;
    double Lx = 300 * dx;
    double Ly =  300 * dy;
    double B = sqrt(2) * lmb;
    double H = sqrt(2) * lmb;
    int ID = 0;
    int nx = static_cast<int>(Lx / dx);
    int ny = static_cast<int>(Ly / dy);

    double x[4], y[4], X[4], Y[4];

    x[0] = 0; y[0] = 0;
    x[1] = 0; y[1] = B;
    x[2] = H; y[2] = B * 0.5;
    x[3] = 0; y[3] = 0;

    /*
r1 : para y > B / 2, y < B, x > 0 e x < H
>>> y = B - 0.5 * B * x / H
r2 : para y > 0, y < B / 2, x > 0 e x < H
>>> y = 0.5 * B * x / H
    */

    for (X[3] = 0; X[3] < H + 0.005 * dx; X[3] += dx) {
        double yf1 = B - 0.5 * B * X[3] / H;
        double yf2 = 0.5 * B * X[3] / H;

        double yf3 = 2 * B;
        double yf4 = - B;

        double N1 = floor(1 + (yf1 - yf2)/ (1.2 * dx));
        double dy1 = (yf1 - yf2) / N1;

        double N2 = floor(1 + (yf2-yf4)/(1.2*dx));
        double dy2 = (yf2 - yf4) / N2;

        double N3 = floor(1 + (yf3-yf1)/(1.2*dx));
        double dy3 = (yf3 - yf1) / N3;

        if (N1 > 0 && dy1 > 0.8 * dx)
        for (Y[3] = yf2; Y[3] <= yf1 + 0.0005 * dy1; Y[3] += dy1) {
            fprintf(fpE, "%e %e\n", X[3], Y[3]);
        }else {
            fprintf(fpE, "%e %e\n", X[3], yf2);
            if (yf1 - 0.5 * B > 0.7 * dx) {
                fprintf(fpE, "%e %e\n", X[3], 0.5 * B);
            }
            fprintf(fpE, "%e %e\n", X[3], yf1);
            printf("%e %e\n", (yf1 - 0.5 * B) / dx,(yf1 - yf2)/dx);
        }

        for (Y[3] = yf4; Y[3] < yf2 - 0.999*dy2; Y[3] += dy2) {
            fprintf(fpE, "%e %e\n", X[3], Y[3]);
        }

        for (Y[3] = yf1 + dy3; Y[3] < yf3; Y[3] += dy3) {
            fprintf(fpE, "%e %e\n", X[3], Y[3]);
        }

    }

    fprintf(fpE, "%e %e\n", H, 0.5 * B);

    std::fclose(fpE);
    std::fclose(fpH);
    std::fclose(fPEC);
    std::fclose(fUPML);
    std::fclose(fFONTE);
}

int main(int argc, char** argv) {
    lane_maxwell::modules::matrix::matrixModuleInit(1e-7);

    point::Points* p = new Points;

    genPoints(p);

    rpim::save("maxwell.rpim", p);

    free(p);
    return 0;
}


