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

#include "test_rpim_3.h"

using namespace lane_maxwell::modules;
using namespace lane_maxwell::modules::matrix;
using namespace lane_maxwell::modules::rpim;
using namespace lane_maxwell::modules::point;
using namespace lane_maxwell::modules::domains;


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

	point::Points* p;
	
	rpim::load("test_rpim.rpim", &p);

	for (size_t i = 0; i < p->size(); i++) {
		Point* pt = p->at(i);
		pt->resizeE(1);
		pt->setE(0, f(pt));
	}


	double err = 0;

	for (int i = 0; i < p->size(); i++) {
		Point* pt = p->at(i);
		double dfx_rpim = pt->calcDFE(0, 0);
		double dfx_real = dfx(pt);

		std::printf("dfx_rpim = %e\ndfx_real = %e\nerro = %e\n\n", dfx_rpim,
		dfx_real, 100*std::fabs(dfx_rpim-dfx_real)/fabs(dfx_real));

		err += 100*std::fabs(dfx_rpim-dfx_real)/fabs(dfx_real);
	}

	for (int i = 0; i < p->size(); i++) {
		Point* pt = p->at(i);
		double dfy_rpim = pt->calcDFE(0, 1);
		double dfy_real = dfy(pt);

		std::printf("dfy_rpim = %e\ndfy_real = %e\nerro = %e\n\n", dfy_rpim,
			dfy_real, 100*std::fabs(dfy_rpim-dfy_real)/fabs(dfy_real));

		err += 100*std::fabs(dfy_rpim-dfy_real)/fabs(dfy_real);
	}

	std::printf("Err medium = %g %g\%\n", err, err / (2*p->size()));

	free(p);
	return 0;
}


