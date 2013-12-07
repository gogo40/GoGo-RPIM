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
@file test_rpim_2.cc

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE) 

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include <algorithm>

#include "test_rpim_2.h"

using namespace lane_maxwell::modules::matrix;
using namespace lane_maxwell::modules::point;
using namespace lane_maxwell::modules::rpim;
using namespace lane_maxwell::modules::domains;

Point* point3D(double x, double y, double z, double eps, double mu, int prop, int ID, char type) {
	Point* po = new Point;
	Vector* p = new Vector(3);

	p->set(0, x);
	p->set(1, y);
	p->set(2, z);

	po->set(ID, eps, mu, prop, type, p);
	return po;
}

void genPoints(Points* const vp, int nx, int ny, int nz, char type) {

	double dx = 0.1;
	double dy = 0.1;
	double dz = 0.1;

	double eps = 8.841941e-12;
	double mu = 1.256637e-6; 
	int prop = 0;
	static int ID = 0;
	//void set(int ID, double eps, double mu, int prop, char type,
				//matrix::Vector* const p)
	for (int ix = 0; ix < nx; ++ix) {
		for (int iy = 0; iy < ny; ++iy) { 
			for (int iz = 0; iz < nz; ++iz) {
				double x, y, z;
				if (type == 'E') {
					x = (2*ix+0)*dx;
					y = (2*iy+0)*dy;
					z = (2*iz+0)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;

					x = (2*ix+0)*dx;
					y = (2*iy+1)*dy;
					z = (2*iz+1)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;

					x = (2*ix+1)*dx;
					y = (2*iy+1)*dy;
					z = (2*iz+0)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;

					x = (2*ix+1)*dx;
					y = (2*iy+0)*dy;
					z = (2*iz+1)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;
				} else {
					x = (2*ix+1)*dx;
					y = (2*iy+1)*dy;
					z = (2*iz+1)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;

					x = (2*ix+1)*dx;
					y = (2*iy+0)*dy;
					z = (2*iz+0)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;

					x = (2*ix+0)*dx;
					y = (2*iy+0)*dy;
					z = (2*iz+1)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;

					x = (2*ix+0)*dx;
					y = (2*iy+1)*dy;
					z = (2*iz+0)*dz;
					vp->push_back(point3D(x, y, z, eps, mu, prop, ID, type));
					++ID;


				}
			}
		}
	}
}

double dumb(const Points* p, int T, char A, char B) {
	double acc = 0;

	for (int i = 0; i < p->size(); ++i) {
		if(p->at(i)->getType() == A) {
			std::vector<double> d;

			for (int j = 0; j < p->size(); ++j) 
				if (i != j && p->at(j)->getType() == B) {
					d.push_back(distE(p->at(i), p->at(j)));
				}

			sort(d.begin(), d.end());
			for (int j = 0; j < T; ++j) {
				acc+=d[j];
			}
		}
		if (i%1000 == 0) std::printf("acc = %e\n", acc);
	}

	return acc;
}

int main() {
	lane_maxwell::modules::matrix::matrixModuleInit(1e-7);

	Points* p = new Points;
	int nx, ny, nz, T;

	std::printf("(nx,ny,nz)=");
	std::scanf("%d%d%d", &nx, &ny, &nz);

	std::printf("T=");
	std::scanf("%d", &T);

	genPoints(p, nx, ny, nz, 'E');
	genPoints(p, nx, ny, nz, 'H');

	std::printf("Points...\n");
	//print(stdout, p);
	
	std::sort(p->begin(), p->end(),PointComp());

	//printAll(stdout, p);

	buildSupportDomain(stdout, 'E', 'H', T, p, 0.01, distE, 1);

	std::printf("r acc dumb = %e\n", dumb(p, T, 'E', 'H'));
	std::printf("r acc line sweep = %e\n", accSD(p));
	
	free(p);
	return 0;
}


