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
@file maxwell.cc
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include "maxwell/maxwell.h"

namespace lane_maxwell { namespace modules { namespace maxwell {
/**
@brief inicializa modo tmz
*/
void initTEZ(point::Points* const pts) {
    point::Points::iterator it = pts->begin(), end = pts->end();

    while (it != end) {
        point::Point* p = *it;
        char type = p->getType();

        if (type == 'E') {
            p->resizeE(2);
            p->resizeD(2);

            p->setE(0, 0);
            p->setD(0, 0);

            p->setE(1, 0);
            p->setD(1, 0);
        } else {
            p->resizeH(1);
            p->resizeB(1);

            p->setH(0, 0);
            p->setB(0, 0);
        }

        ++it;
    }
}

/**
@brief atualiza Hz de ini até fim
*/
void updateHTEZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt) {
    for (size_t i = ini; i < fim; ++i) {
        point::Point* p = pts->at(i);
        int prop = p->getProp();
        char type = p->getType();

        if (type != 'H') continue;

        if (prop & P_METAL) {
            p->setH(0, 0);
        } else if ((prop & P_PEC) == 0) {
            double sum, old_Bz, Bz, Hz;
            double sigma_x = p->getSigma(0);
            double sigma_y = p->getSigma(1);
            double mu = p->getMu();

            old_Bz = p->getB(0);
            sum = p->calcDFE(0, 1) - p->calcDFE(1, 0);

            if (prop & P_UPML) {
                Hz = p->getH(0);
                Bz = p->getB(0);

                Bz = Bz * (1 - sigma_y * cte) / (1 + sigma_y * cte) +
                sum * dt / (1 + sigma_y * cte);

                Hz = Hz * (1-sigma_x * cte) / (1 + sigma_x * cte) +
                (Bz - old_Bz) * 1 / (mu * (1 + sigma_x * cte));

                p->setH(0, Hz);
                p->setB(0, Bz);
            } else {
                Hz = p->getH(0);

                Hz += sum * dt / mu;
                p->setH(0, Hz);
            }
        }
    }
}

/**
@brief atualiza Ex e Ey de ini até fim
*/
void updateETEZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt) {
    for (size_t i = ini; i < fim; ++i) {
        point::Point* p = pts->at(i);
        int prop = p->getProp();
        char type = p->getType();

        if (type != 'E') continue;

        if (prop & P_PECX) {
            p->setE(0, 0);
            p->setD(0, 0);
        }

        if (prop & P_PECY) {
            p->setE(1, 0);
            p->setD(1, 0);
        }

        if ((prop & P_PEC) == 0) {
            double sum, old_Dx, old_Dy, Dx, Ex, Dy, Ey;
            double sigma_x = p->getSigma(0);
            double sigma_y = p->getSigma(1);
            double eps = p->getEps();

            old_Dx = p->getD(0);
            sum = p->calcDFH(0, 1);

            if (prop & P_UPML) {
                Dx = p->getD(0);
                Ex = p->getE(0);

                Dx = Dx + dt*sum;
                Ex = Ex * (1 - sigma_y * cte) / (1 + sigma_y * cte) +
                (Dx * (1 + sigma_x * cte) - old_Dx * (1 - sigma_x* cte)) *
                (1 / (eps * (1 + sigma_y * cte)));

                p->setD(0, Dx);
                p->setE(0, Ex);
            } else {
                Ex = p->getE(0);
                Ex += sum * dt / eps;

                p->setE(0, Ex);
            }

            old_Dy = p->getD(1);
            sum = p->calcDFH(0, 0);

            if (prop & P_UPML) {
                Dy = p->getD(1);
                Ey = p->getE(1);

                Dy = Dy * (1 - sigma_x * cte) / (1 + sigma_x * cte) -
                sum* dt / (1 + sigma_x * cte);
                Ey += (Dy * (1 + sigma_y * cte) - old_Dy * (1 - sigma_y * cte)) / eps;

                p->setD(1, Dy);
                p->setE(1, Ey);
            } else {
                Ey = p->getE(1);
                Ey -= sum * dt / eps;

                p->setE(1, Ey);
            }
        }
    }
}

void initTMZ(point::Points* const pts) {
    point::Points::iterator it = pts->begin(), end = pts->end();

    while (it != end) {
        point::Point* p = *it;
        char type = p->getType();

        if (type == 'E') {
            p->resizeE(1);
            p->resizeD(1);

            p->setE(0, 0);
            p->setD(0, 0);
        } else {
            p->resizeH(2);
            p->resizeB(2);

            p->setH(0, 0);
            p->setH(1, 0);

            p->setB(0, 0);
            p->setB(1, 0);
        }

        ++it;
    }
}

void updateETMZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt) {
    for (size_t i = ini; i < fim; ++i) {
        point::Point* p = pts->at(i);
        int prop = p->getProp();
        char type = p->getType();

        if (type != 'E') continue;

        if (prop & P_METAL) {
            p->setE(0, 0);
        } else if ((prop & P_PEC) == 0) {
            double sum, old_Dz, Dz, Ez;
            double sigma_x = p->getSigma(0);
            double sigma_y = p->getSigma(1);
            double eps = p->getEps();

            old_Dz = p->getD(0);
            sum = p->calcDFH(1, 0) - p->calcDFH(0, 1);

            if (prop & P_UPML) {
                Ez = p->getE(0);
                Dz = p->getD(0);

                Dz = Dz * (1 - sigma_y * cte) / (1 + sigma_y * cte) +
                sum * dt / (1 + sigma_y * cte);

                Ez = Ez * (1-sigma_x * cte) / (1 + sigma_x * cte) +
                (Dz - old_Dz) * 1 / (eps * (1 + sigma_x * cte));

                p->setE(0, Ez);
                p->setD(0, Dz);
            } else {
                Ez = p->getE(0);

                Ez += sum * dt / eps;
                p->setE(0, Ez);
            }
        }
    }
}


void updateHTMZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt) {
    for (size_t i = ini; i < fim; ++i) {
        point::Point* p = pts->at(i);
        int prop = p->getProp();
        char type = p->getType();

        if (type != 'H') continue;

        if (prop & P_PECX) {
            p->setH(0, 0);
            p->setB(0, 0);
        } else if (prop & P_PECY) {
            p->setH(1, 0);
            p->setB(1, 0);
        } else if ((prop & P_PEC) == 0) {
            double sum, old_Bx, old_By, Bx, Hx, By, Hy;
            double sigma_x = p->getSigma(0);
            double sigma_y = p->getSigma(1);
            double mu = p->getMu();

            old_Bx = p->getB(0);
            sum = p->calcDFE(0, 1);

            if (prop & P_UPML) {
                Bx = p->getB(0);
                Hx = p->getH(0);

                Bx = Bx - dt*sum;
                Hx = Hx * (1 - sigma_y * cte) / (1 + sigma_y * cte) +
                (Bx * (1 + sigma_x * cte) - old_Bx * (1 - sigma_x* cte)) *
                (1 / (mu * (1 + sigma_y * cte)));

                p->setB(0, Bx);
                p->setH(0, Hx);
            } else {
                Hx = p->getH(0);
                Hx -= sum * dt / mu;

                p->setH(0, Hx);
            }

            old_By = p->getB(1);
            sum = p->calcDFE(0, 0);

            if (prop & P_UPML) {
                By = p->getB(1);
                Hy = p->getH(1);

                By = By * (1 - sigma_x * cte) / (1 + sigma_x * cte) +
                sum* dt / (1 + sigma_x * cte);
                Hy += (By * (1 + sigma_y * cte) - old_By * (1 - sigma_y * cte)) / mu;

                p->setB(1, By);
                p->setH(1, Hy);
            } else {
                Hy = p->getH(1);
                Hy += sum * dt / mu;

                p->setH(1, Hy);
            }
        }
    }
}
}}}  // lane_maxwell::modules::maxwell

