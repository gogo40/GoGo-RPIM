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


    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;
    double eps = 8.854187817620e-12;
    double mu = 4 * acos(-1) * 1e-7;

    std::FILE* fpE = std::fopen("pointsE.dat", "w+");
    std::FILE* fpH = std::fopen("pointsH.dat", "w+");
    std::FILE* fPEC = std::fopen("pointsPEC.dat", "w+");
    std::FILE* fUPML = std::fopen("pointsUPML.dat", "w+");
    std::FILE* fFONTE = std::fopen("pointsFONTE.dat", "w+");

    double co = 3e8;
    double fmax = 4.9e9;
    double lmb = co / fmax;
    double dx = lmb / 40;
    double dy = lmb / 40;

    printf("dx = %e\n", dx);
    double pi = std::acos(-1.0);
    double ro = 7.62e-2;
    double r1 = ro - 2.9e-2;
    double r2 = ro + 2.9e-2;

    double dr;

    dx = 1e-3;
    dr = dx;
     printf("dt = %e <|> \n", 0.4 * dx * sqrt(eps * mu));

    double dphi =  2 * dx / 7.62e-2;
    int ix = 0, iy = 0, key = 0;
    int ID = 0;
    for (double r = r1 - 20 * dr; r < 1.000001 * r2 + 30 * dr; r += dr ) {


         ix++;
         for (double phi = 0; phi < 0.99 *pi * 0.5; phi += dphi ) {
            double x, y;
            double x1 = -r * cos(phi);
            double y1 = r * sin(phi);

            x = x1 + 0.1053 + 50 * dx;
            y = y1 + 0.06 + 100 *dx;

            char type = (key == 0)?'E':'H';

            if (iy % 2 == ix % 2) {
                type = (key == 0)?'H':'E';
              //  if (iy % 2 == 1) type = 'I';
            }
            iy++;

            if (type == 'E') fprintf(fpE, "%e %e\n", x, y);
            else if (type == 'H') fprintf(fpH, "%e %e\n", x, y);

            vp->push_back(point2D(x, y, eps, mu, (r >= r1 && r <= 1.000001 * r2)?0:P_PEC|P_METAL,
                                  ID, type, 0, 0));

            if (r  < r1 || r > 1.000001 * r2) fprintf(fPEC, "%e %e\n", x, y);


            ++ID;

        }

        double xo = - (r), yo = (r);


        for (double y = -20e-2 - 20 * dr; y < 0; y += dr) {
            double X = xo + 0.1053 + 50 * dx;
            double Y = y + 0.06 + 100 *dx;
            int prop = (r >= r1 && r <= r2)?0:(P_PEC|P_METAL), isUPML = 0, isSource = 0;
            double sigma_x = calcSigma(X, dr, 100, &isUPML);
            double sigma_y = calcSigma(Y, dr, 100, &isUPML);



            if (r < r1 || r >1.000001 * r2) {
                prop |= P_FIXO;
                fprintf(fPEC, "%e %e\n", X, Y);
            }


            if (isUPML) prop |= P_FIXO;

            char type = (key == 0)?'E':'H';


            if (iy % 2 == ix % 2) {
                type = (key == 0)?'H':'E';
                //if (iy % 2 == 0) type = 'I';
            }


            if (fabs(y + 20 * dr + 0.1) < 1e-5 && !(prop & P_METAL) && type == 'E') {
               // printf("s %e %e\n", X, Y);
                prop |= P_FONTE;
                fprintf(fFONTE, "%e %e\n", X, Y);
            }

            iy++;



            if (fabs(xo+7.62e-2) < 1e-5 && fabs(y+7.1e-2) < 1e-5 && type == 'E') {
                printf("m %e %e\n", X, Y);
                prop |= P_MED;
            }

            if (type == 'E') fprintf(fpE, "%e %e\n", X, Y);
            else if (type == 'H') fprintf(fpH, "%e %e\n", X, Y);




            vp->push_back(point2D(X, Y, eps, mu, (prop | isUPML),
                                  ID, type, sigma_x, sigma_y));
            if (isUPML) fprintf(fUPML, "%e %e\n", X, Y);

            ++ID;
        }

        for (double x = 0; x <= 10e-2 + 30 * dr; x += dr) {
            double X = x + 0.1053 + 50 * dx;
            double Y = yo + 0.06 + 100 *dx;
            int prop = (r >= r1 && r <= r2)?0:P_PEC|P_METAL, isUPML = 0;
            double sigma_x = calcSigma(X, dr, 20e-2, &isUPML);
            double sigma_y = calcSigma(Y, dr, 100, &isUPML);


            char type = (key == 0)?'E':'H';

            if (iy % 2 == ix % 2) {
                type = (key == 0)?'H':'E';
                //if (iy % 2 == 1) type = 'I';
            }

            iy++;


            if (r < r1 || r >1.000001 * r2){
                prop |= P_FIXO;
                fprintf(fPEC, "%e %e\n", X, Y);
            }

            if (isUPML) prop |= P_FIXO;

            if (type == 'E') fprintf(fpE, "%e %e\n", X, Y);
            else if (type == 'H') fprintf(fpH, "%e %e\n", X, Y);


            vp->push_back(point2D(X, Y, eps, mu, prop | isUPML,
                                  ID, type, sigma_x, sigma_y));

            if (isUPML) fprintf(fUPML, "%e %e\n", X, Y);

            ++ID;
        }

    }

    printf("N points = %d\n", ID);
    std::fclose(fpE);
    std::fclose(fpH);
    std::fclose(fPEC);
    std::fclose(fUPML);
    std::fclose(fFONTE);

    point::coloumb2D(vp, dx, 200, 200, 5*dx);

    FILE* fColoumbE = fopen("ColoumbE.dat", "w+");
    FILE* fColoumbH = fopen("ColoumbH.dat", "w+");

    int N = vp->size();

    for (int i = 0; i < N; ++i) {
        point::Point* pt = vp->at(i);

        if (pt->getType() == 'E') {
            fprintf(fColoumbE, "%e %e\n", pt->getX(0), pt->getX(1));
        } else {
            fprintf(fColoumbH, "%e %e\n", pt->getX(0), pt->getX(1));
        }
    }

    fclose(fColoumbE);
    fclose(fColoumbH);
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

    rpim::save("maxwell.rpim", p);

    free(p);
    return 0;
}


