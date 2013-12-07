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

double calcSigma(double x,
                 int *prop,
                 double x_upml_1,
                 double x_upml_2,
                 double tam_upml,
                 double dx, double smax) {
    double m = 4;
    double sigma_max_x = 15.0;
    double tam_upml_x = tam_upml;
    double sigma_x = 0;

    sigma_max_x  = (1 + m) / (15 * M_PI * 15 * dx);
    sigma_max_x  = smax;

    if (x <= x_upml_1 + 1e-4) {
        double v = fabs(x_upml_1 - x) / tam_upml_x; 
                
        sigma_x = pow(v, m) * sigma_max_x;

        //sigma_x = (exp(v * m) - 1) * sigma_max_x;

        *prop |= P_UPML;
    } else if (x >= x_upml_2 + 1e-4) {
        double v = fabs(x_upml_2 - x) / tam_upml_x; 
                
        sigma_x = pow(v, m) * sigma_max_x;
        
        //sigma_x = (exp(v * m) - 1) * sigma_max_x;

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

double dist(double x, double y, double x1, double y1)  {
    double dx = x - x1;
    double dy = y - y1;

    return sqrt(dx * dx + dy * dy);
}

bool inTri(double x, double y, double xo, double yo, double B, double H) {
    double px = x - xo;
    double py = y - yo;
    
   // if (px < 2 * B * B and px > 0 and py < 2 * B * B and py > 0) return true; return false;
    
/*
  (0, 0) (0, B) (H, B/2)
*/
    double vx[] = {0, 0, H};
    double vy[] = {0, B, B / 2};

    double a1 = 0.5 * fabs(px * vy[0] + py * vx[1] + vx[0] * vy[1] - vx[1] * vy[0] - vy[1] * px - vx[0] * py);
    double a2 = 0.5 * fabs(px * vy[0] + py * vx[2] + vx[0] * vy[2] - vx[2] * vy[0] - vy[2] * px - vx[0] * py);
    double a3 = 0.5 * fabs(px * vy[2] + py * vx[1] + vx[2] * vy[1] - vx[1] * vy[2] - vy[1] * px - vx[2] * py);
    double A = B * H * 0.5;

    if (fabs(a1 + a2 + a3 - A) < 1e-8) return true;

    return false;
}


void genPoints(Points* vp, double dt,
               double x_upml_1, double x_upml_2,
               double y_upml_1, double y_upml_2,
               double tam_upml, const char* fe,
               const char* fh,
               double xf, double yf,
               double xm, double ym,
               double xt, double yt,
               double B, double H, 
               double smax, int okTri, double scfact) {

    FILE* fE = fopen(fe, "r");
    FILE* fH = fopen(fh, "r");


    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;
    double sx;
    double sy;
    int prop;
    double d;
    int ID = 0;

    int nupml = 20;
    double dx = tam_upml / nupml;

    FILE* f1 = fopen("E.out", "w+");
    FILE* fs = fopen("B.poly", "w+");

    std::vector<double> Px;
    std::vector<double> Py;

    double raio = 0.15;
    double xc = xt + 0.5 * B, yc = yt + 0.5 * B;

    for (double x = x_upml_1- nupml * 2 * dx; x < x_upml_2 + nupml * 2 *dx; x += dx) {
        for (double y = y_upml_1- nupml * 2 * dx; y < y_upml_2 + nupml * 2 * dx; y += dx) {

            prop = 0; sx = sy = 0;

            if (x <= x_upml_1+dx || y <= y_upml_1+dx || x > x_upml_2 || y > y_upml_2) {
                prop |= P_PEC;
                prop |= P_METAL;
            }

            double x1 = xc, y1 = yc;
            double d = dist(x, y, x1, y1);

    
            if (! ( 
		    (inTri(x, y, xt, yt, B, H) and okTri) and//or 
	               fabs(d-raio) < 2e-3
		  )) 
            {
                Px.push_back(x); Py.push_back(y);
                fprintf(f1, "%e %e\n", x, y);
            }
        }
    }

    double phi = 2 * M_PI;
    int NPART = 500;
    double dphi = phi / NPART;

    FILE* fmedmap = fopen("medmap.pos","w+");
    fprintf(fmedmap,"%lf %lf\n", xc, yc);

    for (int i = 0; i < NPART - 1; ++i) {
        double x = xc + raio * cos(i * dphi), y = yc + raio * sin(i * dphi);
        fprintf(f1, "%e %e\n", x, y);
        //Px.push_back(x); Py.push_back(y);
        fprintf(fmedmap, "Ez_%g_%g.dat med_%d.dat %lf\n", scfact*x, scfact*y, i, i*dphi / (2*M_PI) * 360);
    }

    fclose(fmedmap);

    double sigma_max_x  = (1 + 5) / (15 * M_PI * 15 * dx);
    printf("smax = %e\n", sigma_max_x);

    /*Caixa externa*/

    int C[20];

    Px.push_back(x_upml_1 - 2 * nupml *dx);
    Py.push_back(y_upml_1 - 2 * nupml *dx);
    C[0] = Px.size();

    Px.push_back(x_upml_2 + 2 * nupml *dx);
    Py.push_back(y_upml_2 + 2 * nupml *dx);
    C[1] = Px.size();


    Px.push_back(x_upml_1 - 2 * nupml *dx);
    Py.push_back(y_upml_2 + 2 * nupml *dx);
    C[2] = Px.size();

    Px.push_back(x_upml_2 + 2 * nupml *dx);
    Py.push_back(y_upml_1 - 2 * nupml *dx);
    C[3] = Px.size();



    /*Estrutura*/

    /*Ponto de fonte*/
    //Px.push_back(xf);
    //Py.push_back(yf);

    //Px.push_back(xm);
    //Py.push_back(ym);

//Triangulo

    if (okTri) {
        //Px.push_back(xt);
        //Py.push_back(yt);
        C[4] = Px.size();

        //Px.push_back(xt);
        //Py.push_back(yt + B);
        C[5] = Px.size();

        //Px.push_back(xt + H);
        //Py.push_back(yt + B * 0.5);
        C[6] = Px.size();
    }

/*
//caixa envolvendo a fonte
    Px.push_back(xf +  H );
    Py.push_back(yf -  0*B);
    C[7] = Px.size();

    Px.push_back(xf +  H);
    Py.push_back(yf +  B);
    C[8] = Px.size();

    Px.push_back(xf -  0*H );
    Py.push_back(yf -  0*B);
    C[9] = Px.size();

    Px.push_back(xf -  0*H);
    Py.push_back(yf +  B);
    C[10] = Px.size();
*/

    /*Gerando arquivo da estrutura*/

    fprintf(fs, "%d 2 0 0\n", Px.size());

    for (int i = 0; i < Px.size(); ++i) {
        fprintf(fs, "%d %e %e\n", i + 1, Px[i], Py[i]);
    }

    if (okTri) {
	fprintf(fs, "4 0\n");//se tri
        //fprintf(fs, "7 0\n");
    } else {
        fprintf(fs, "4 0\n");
    }

    fprintf(fs, "1 %d %d\n", C[0], C[2]);
    fprintf(fs, "2 %d %d\n", C[2], C[1]);
    fprintf(fs, "3 %d %d\n", C[1], C[3]);
    fprintf(fs, "4 %d %d\n\n", C[3], C[0]);

    //Triangulo
 
    if (okTri) {  
        //fprintf(fs, "5 %d %d\n", C[4], C[5]);
        //fprintf(fs, "6 %d %d\n", C[4], C[6]);
        //fprintf(fs, "7 %d %d\n", C[6], C[5]);
    }

/* //caixa envolvendo no triangulo

    

    fprintf(fs, "8 %d %d\n", C[7], C[8]);    
    fprintf(fs, "9 %d %d\n", C[8], C[10]);
    fprintf(fs, "10 %d %d\n", C[10], C[9]);
    fprintf(fs, "11 %d %d\n", C[9], C[7]);

//*/

    fprintf(fs, "0\n\n");
    fclose(fs);
    fclose(f1);

    if (!fE or !fH) return;

    printf("xc = %lf yc = %lf\n", xc, yc);
    FILE* fff = fopen("FILEFONTE.POS","w+");
    fprintf(fff,"--BEGIN--\n");
    while (!feof(fE)) {
        double x, y, x1, y1;

        if (fscanf(fE, "%lf%lf", &x, &y) != 2) {
            break;
        }

        prop = 0;
        sx = calcSigma(x, &prop, x_upml_1, x_upml_2, tam_upml, dx, smax);
        sy = calcSigma(y, &prop, y_upml_1, y_upml_2, tam_upml, dx, smax);

        if (x < x_upml_1 - nupml * dx || y < y_upml_1 - nupml * dx ||
            x > x_upml_2 + nupml * dx || y > y_upml_2 + nupml * dx) {
            prop |= P_PEC;
            prop |= P_METAL;
        }

        if (inTri(x, y, xt, yt, B, H) and okTri) {
            prop |= P_PEC;
            prop |= P_METAL;
        }

        /*****************************************************/
        x1 = xf; y1 = yf;
        d = dist(x, y, x1, y1);

        if (fabs(x-x1) < 1e-4 //and fabs(y-y1) < 1e-4
	  and ! ((prop&P_METAL) || (prop&P_PEC))) {
            fprintf(fff, "fonte: %e %e\n", x, y, d);
            prop |= P_FONTE;
        }    

        x1 = xc; y1 = yc;
        d = dist(x, y, x1, y1);


        if  (d > raio - 1e-3 and d < raio + 1e-3) {
           prop |= P_MED;
        }

        /*****************************************************/
        //sx = sy = 0;
        vp->push_back(point2D(scfact * x, scfact * y, eps0, mu0, prop, ID, 'E', sx, sy));
        ++ID;
    }
    fprintf(fff,"--BEGIN--\n");
    fclose(fff);
    fclose(fE);

    while (!feof(fH)) {
        double x, y;

        if (fscanf(fH, "%lf%lf", &x, &y) != 2) {
            break;
        }

        prop = 0;
        sx = calcSigma(x, &prop, x_upml_1, x_upml_2, tam_upml, dx, smax);
        sy = calcSigma(y, &prop, y_upml_1, y_upml_2, tam_upml, dx, smax);

        if (x < x_upml_1 - nupml * dx || y < y_upml_1 - nupml * dx ||
            x > x_upml_2 + nupml * dx || y > y_upml_2 + nupml * dx) {
            prop |= P_PEC;
            prop |= P_METAL;
        }

        if (inTri(x, y, xt, yt, B, H) and okTri) {
            prop |= P_PEC;
            prop |= P_METAL;
        }
        /*****************************************************/


        /*****************************************************/
        //sx = sy = 0;
        vp->push_back(point2D(scfact * x, scfact * y, eps0, mu0, prop, ID, 'H', sx, sy));
        ++ID;
    }


    fclose(fH);
    rpim::save("maxwell.rpim", vp);

    free(vp);
}

int main(int argc, char** argv) {

    if (argc != 20) {
        printf("%s <dx> <x_upml_1> <x_upml_2> <y_upml_1> <y_upml_2> <tam_upml> <file E> <file H> <xf> <yf>"
               " <xm> <ym> <xt> <yt> <B> <H> <smax> <okTri> <scale factor>\n", argv[0]);
        return 0;
    }

    lane_maxwell::modules::matrix::matrixModuleInit(1e-7);

    point::Points* p = new Points;
    double dt, dx, x_upml_1, x_upml_2, y_upml_1, y_upml_2, tam_upml;
    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;

    sscanf(argv[1], "%lf", &dx);

    dt = 0.4 * dx * sqrt(eps0 * mu0);

    printf("dt = %e\n", dt);

    sscanf(argv[2], "%lf", &x_upml_1);
    sscanf(argv[3], "%lf", &x_upml_2);

    sscanf(argv[4], "%lf", &y_upml_1);
    sscanf(argv[5], "%lf", &y_upml_2);

    sscanf(argv[6], "%lf", &tam_upml);

    double x, y;
    double xm, ym;
    double xt, yt;
    double B, H;
    double smax;
    int okTri;
    double scfact;

    sscanf(argv[9], "%lf", &x);
    sscanf(argv[10], "%lf", &y);

    sscanf(argv[11], "%lf", &xm);
    sscanf(argv[12], "%lf", &ym);


    sscanf(argv[13], "%lf", &xt);
    sscanf(argv[14], "%lf", &yt);

    sscanf(argv[15], "%lf", &B);
    sscanf(argv[16], "%lf", &H);


    sscanf(argv[17], "%lf", &smax);
    sscanf(argv[18], "%d", &okTri);
    sscanf(argv[19], "%lf", &scfact);

    genPoints(p, dt, x_upml_1, x_upml_2,
              y_upml_1, y_upml_2,
              tam_upml, argv[7], argv[8], 
              x, y, xm, ym, 
              xt, yt, B, H, smax, okTri, scfact);



    return 0;
}



