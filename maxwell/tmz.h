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
@file tmz.h
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#ifndef MAXWELL_TMZ_H_
#define MAXWELL_TMZ_H_

#include "maxwell/maxwell.h"

namespace lane_maxwell { namespace modules { namespace maxwell {
typedef double (*PFSource)(const double& t, const double& dt, const point::Point* pts);

enum SourceType {
    HARD = 0,
    SOFT = 1,
    TYPE3 = 2
};

bool tmz(point::Points* pts, const double& Ts, const double& dt,
    const double& cte, double resX, double resY,
    bool videoGen, int nf, bool isLog, int lecont,
    size_t NTHREADS, PFSource source, int sourceType, double duracaoPulso, double wfreq);
}}}  // lane_maxwell::modules::tmz
#endif  // MAXWELL_TMZ_H_

