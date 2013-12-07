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
@file maxwell.h
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/


#ifndef MAXWELL_MAXWELL_H_
#define MAXWELL_MAXWELL_H_

#include "rpim/rpim.h"

namespace lane_maxwell { namespace modules { namespace maxwell {
enum MaxwellProp {
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de PEC
*/
    P_PEC = 1<<0,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de medição
*/
    P_MED = 1<<1,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de fonte
*/
    P_FONTE = 1<<2,
    /**
@brief Palavra binária usada para indicar que um ponto é metal
*/
    P_METAL = 1<<3,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de PEC na direção X
*/
    P_PECX = 1<<4,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de PEC na direção Y
*/
    P_PECY = 1<<5,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de FONTE na direção X
*/
    P_DIRFX = 1<<6,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de FONTE na direção Y
*/
    P_DIRFY = 1<<7,
    /**
@brief Palavra binária usada para indicar que um ponto é ponto de FONTE na direção Z
*/
    P_DIRFZ = 1<<8,
    /**
@brief Palavra binária usada para indicar que um ponto pertence a linha de visualizacao de campo X
*/
    P_LCX = 1<<9,
    /**
@brief Palavra binária usada para indicar que um ponto pertence a linha de visualizacao de campo Y
*/
    P_LCY = 1<<10,
    /**
@brief Palavra binária usada para indicar que um ponto pertence a linha de visualizacao de campo Z
*/
    P_LCZ = 1<<11,

    P_COR = 1<<12,
    P_PECZ = 1<<13,
    P_FDTD = 1<<14,
    /**
@brief Fonte excita E?
*/
    P_FONTEE = 1<<15,

    /**
@brief Palavra binária usada para indicar que ponto utiliza coordenadas locais
*/
    P_LC = 1<<16,

    /**
@brief Palavra binária usada para indicar que ponto é de UPML
*/
    P_UPML = 1 << 17,

    /**
@brief Palavra binária usada para indicar que é preciso gravar
dados sobre o dominio de suporte do ponto atual
*/
    P_SDEBUG = 1 << 18,

    /**
@brief Palavra binária usada para indicar que o ponto é fixo
*/
    P_FIXO = 1 << 19
};

/**
@brief inicializa modo tmz
*/
void initTMZ(point::Points* const pts);

/**
@brief atualiza Ez de ini até fim
*/
void updateETMZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt);

/**
@brief atualiza Hx e Hy de ini até fim
*/
void updateHTMZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt);

/**
@brief inicializa modo tmz
*/
void initTEZ(point::Points* const pts);

/**
@brief atualiza Hz de ini até fim
*/
void updateHTEZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt);

/**
@brief atualiza Ex e Ey de ini até fim
*/
void updateETEZ(point::Points* const pts, size_t ini, size_t fim, const double& cte, const double& dt);
}}}  // lane_maxwell::modules::maxwell
#endif  // MAXWELL_MAXWELL_H_
