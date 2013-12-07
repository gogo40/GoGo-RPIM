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
@file supportdomain.h

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#ifndef DOMAINS_SUPPORTDOMAIN_H_
#define DOMAINS_SUPPORTDOMAIN_H_

#include "domains/domains.h"

namespace lane_maxwell { namespace modules { namespace domains {
//Se A é igual a 'a' então todos pontos proximos são considerados
double buildSupportDomain(std::FILE* f, char B, char A, size_t T,
                            const point::Points* P,  double fator, FPDist d, double fator_osc,
                            bool useR = false, double rd = 0);
}}}  // lane_maxwell::modules::domains

#endif  // DOMAINS_SUPPORTDOMAIN_H_

