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
@file maxwellgendom.h
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/


#ifndef MAXWELL_MAXWELLGENDOM_H_
#define MAXWELL_MAXWELLGENDOM_H_

#include "maxwell/maxwell.h"
#include "lane/lane.h"
namespace lane_maxwell { namespace modules { namespace maxwell {

bool calcDF(
            point::Points* pts, size_t T, size_t NTHREADS,
            domains::FPDist D, double beta,
            bool useR, double rd,
            rpim::RadialFunction R,
            rpim::RadialFunction dR, bool isauto, double betamin,
            double betamax, double errbeta, int NP, double fmax);

}}}  // lane_maxwell::modules::maxwell
#endif  // MAXWELL_MAXWELLGENDOM_H_
