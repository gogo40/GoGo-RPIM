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

A New FDTD Algorithm Based
   on Alternating-Direction
      Implicit Method

EXAMPLE A - FREE SPACE

*/

#include <algorithm>
#include <ctime>
#include <cstdlib>

#include "maxwell/tmz.h"

using namespace lane_maxwell::modules;
using namespace lane_maxwell::modules::maxwell;
using namespace lane_maxwell::modules::matrix;
using namespace lane_maxwell::modules::rpim;
using namespace lane_maxwell::modules::point;
using namespace lane_maxwell::modules::domains;

Point* point2D(double x, double y, double eps, double mu,
	int prop, int ID, char type) {
	Point* po = new Point;
	Vector* p = new Vector(2);

	p->set(0, x);
	p->set(1, y);

	po->set(ID, eps, mu, prop, type, p);

	return po;
}

void genPoints(Points* const vp, int N, double lmb) {
	bool key = true;
	std::FILE* fpe = fopen("pointse.dat", "w+");
	std::FILE* fph = fopen("pointsh.dat", "w+");
	int ID = 0;
	int sx=0, sy=0;
if(0){
	for (int i = 0; i < N; ++i){ 
		for (int j = 0; j < N; ++j){
				key = (sx==sy);

				vp->push_back(point2D(i * lmb, j * lmb, 0, 0, 0, ID, (key)?'E':'H'));
				if(key) {
					fprintf(fpe, "%e %e\n", i*lmb, j*lmb);
					key = false;
				} else {
					fprintf(fph, "%e %e\n", i*lmb, j*lmb);
					key = true;
				}
				sy = (sy+1)%2;
				++ID;
		}
		sx = (sx+1)%2;
	}
}



if(1){
	int N, i;
	i = 0;
	ID = 0;

	std::FILE* FH = fopen("Malha_H.txt","r");
	//std::rewind(FH);
	fscanf(FH,"%d\n",&N);
	i=0;
	while (i < N) {
		double x, y;
		fscanf(FH,"%lf%lf\n", &x, &y);
		vp->push_back(point2D(x, y, 0, 0, 0, ID,'H'));
		++ID;
		++i;
	}
	fclose(FH);


	std::FILE* FE = fopen("Malha_E.txt","r");
	//std::rewind(FE);
	fscanf(FE,"%d\n",&N);
	i = 0;
	while (i < N) {
		double x, y;
		fscanf(FE,"%lf%lf\n", &x, &y);
		vp->push_back(point2D(x, y, 0, 0, 0, ID,'E'));
		++ID;
		++i;
	}
	fclose(FE);
	
	
}



	fclose(fpe);
	fclose(fph);
}


void plot(const point::Point* p, int resx, int resy) {
	double xmin, xmax, Lx;
	double ymin, ymax, Ly;

	std::vector<std::vector<char> > grid(resx+2);

	for (size_t k = 0; k < resx+2; ++k) grid[k].resize(resy+2);

	for (size_t k = 0; k < resx+2; ++k)
		for (size_t i = 0; i < resy+2; ++i)
			grid[k][i] = '.';

	for (int k = 0; k < p->sizeSD(); ++k) {
		const point::Point* pt = p->getSD(k);

		if (k == 0 || xmin > pt->getX(0)) {
			xmin = pt->getX(0);
		}

		if (k == 0 || ymin > pt->getX(1)) {
			ymin = pt->getX(1);
		}

		if (k == 0 || xmax < pt->getX(0)) {
			xmax = pt->getX(0);
		}

		if (k == 0 || ymax < pt->getX(1)) {
			ymax = pt->getX(1);
		}
	}

	Lx = xmax - xmin;
	Ly = ymax - ymin;

	for (int k = 0; k < p->sizeSD(); ++k) {
		const point::Point* pt = p->getSD(k);
		double x = pt->getX(0);
		double y = pt->getX(1);
		int ix = static_cast<int>(abs( (x - xmin) / Lx) * resx);
		int iy = static_cast<int>(abs( (y - ymin) / Ly) * resy);
		grid[ix][iy] = pt->getType();
	}

	for (size_t k = 0; k < resx+2; ++k) {
		for (size_t i = 0; i < resy+2; ++i)
			std::printf("%c", grid[k][i]);
		std::printf("\n");
	} 
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
	
	printf("Numero de pontos=");
	scanf("%d", &N);

	printf("Lambda / fat = ");
	scanf("%lf", &fat);

	point::Points* p = new Points;

	printf("d = %e\n", lmb/fat);
	genPoints(p, N, lmb / fat);

	T = 12;

	double betamin, betamax, errbeta;
	double nd, delta = lmb / fat;
	double raio;
	int NP;

	printf("betamin betamax errbeta NP = ");
	scanf("%lf%lf%lf%d", &betamin, &betamax, &errbeta, &NP);

	printf("raio=");
	scanf("%lf", &raio);

	std::sort(p->begin(), p->end(), point::PointComp());

	lane_maxwell::modules::domains::buildSupportDomain(
	stdout, 'E', 'H', T, p, 0.01,
	lane_maxwell::modules::point::distE, 1.1, true, raio);

	lane_maxwell::modules::domains::buildSupportDomain(
	stdout, 'H', 'E', T, p, 0.01,
	lane_maxwell::modules::point::distE, 1.1, true, raio);

	point::Points::iterator it = p->begin(), end = p->end();


	int mT = 0;
	
	for (int i = 0; i < p->size(); ++i) {
		if ((*p)[i]->sizeSD() > mT) mT = (*p)[i]->sizeSD();
	}

	std::vector<double> F(mT);
	std::vector<double> beta(2);

	rpim::CallDF* cdf = new rpim::CallDF(mT, 2);

	bool status;
	int i = 0;

	static char fnout[1000];
	static char fnout2[1000];

	sprintf(fnout, "fbeta_N=%d_fmax=%e_fat=%e_T=%d.dat", N, fmax, fat, T);
	sprintf(fnout2, "fbeta_N=%d_fmax=%e_fat=%e_T=%d_dist.dat", N, fmax, fat, T);

	std::FILE* out = std::fopen(fnout, "w+");
	std::FILE* out2 = std::fopen(fnout2, "w+");
	std::FILE* ferr = std::fopen("err", "w+");
	std::FILE* fserr = std::fopen("serr.txt", "w+");
	std::FILE* fsd = std::fopen("sd.txt", "w+");
	std::FILE* fsdc = std::fopen("sdc.txt", "w+");

	int nI5 = 0;
	int nI1 = 0;
	int nI01 = 0;
	int nI001 = 0;
	int kp = 0;
 	while (it != end) {
		double err1, err2;
		T = (*it)->sizeSD();
		beta[0] = rpim::bestBeta(0, &status, i, T, p, cdf, rpim::phi, rpim::dphi, fmax,
			betamin, betamax, NP, errbeta, &F, &err1);

		if (status) {
			beta[1] = rpim::bestBeta(1, &status, i, T, p, cdf, rpim::phi, rpim::dphi, fmax,
			betamin, betamax, NP, errbeta, &F, &err2);
		}

		double rmax = (*it)->getRMax();
			
		if (status) {
/*			double fn = rmax * std::sqrt(2);

			for (size_t k = 0; k < (*it)->sizeSD(); ++k) {
				const point::Point* pt = (*it)->getSD(k);
				fprintf(out,"%e %e ", pt->getX(0) / fn, pt->getX(1) / fn);
			}


			for (size_t k = 0; k < (*it)->sizeSD(); ++k) {
				const point::Point* pt = (*it)->getSD(k);
				
				if (kp % 300 == 0) {
					fprintf(fsd, "%e %e\n", pt->getX(0), pt->getX(1));
				}
			}
			fprintf(fsdc, "%e %e\n", (*it)->getX(0), (*it)->getX(1));
*/


			fprintf(out,"%e %e\n", beta[0], beta[1]);
			fprintf(out2,"%d %e %e\n", i, beta[0], beta[1]);

			fprintf(ferr,"%e %e\n", err1, err2);

			printf("<>%d %d %e %c\n", T, i, rmax, (*it)->getType());

			if (err1 < 5 && err2 < 5) nI5++;
			if (err1 < 1 && err2 < 1) nI1++;
			if (err1 < 0.1 && err2 < 0.1) nI01++;
			if (err1 < 0.01 && err2 < 0.01) nI001++;

		}else {
			printf("+%d %d %e %e [%e %e]\n",i, T, err1, err2, beta[0], beta[1]);
			err1 = err2 =10;

			/*
			printf("%d %d %e %e %e %c\n", T, i, rmax, (*it)->getX(0), (*it)->getX(1), (*it)->getType());
			for (size_t k = 0; k < (*it)->sizeSD(); ++k) {
				const point::Point* pt = (*it)->getSD(k);
				printf("(%e,%e) ", pt->getX(0), pt->getX(1));

				if (kp % 300 == 0) {
					fprintf(fsd, "%e %e\n", pt->getX(0), pt->getX(1));
				}
			}
			fprintf(fsdc, "%e %e\n", (*it)->getX(0), (*it)->getX(1));
			printf("\n");*/
			//plot((*it), 100, 100);
		}
		fprintf(fserr, "%e %e %e\n", (*it)->getX(0), (*it)->getX(1), std::sqrt(err1*err1+err2*err2));
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
	std::fclose(out);
	std::fclose(out2);
	std::fclose(fsd);
	std::fclose(fsdc);

	free(p);

	return 0;
}


