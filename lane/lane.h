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
@file lane.h
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/


#ifndef LANE_LANE_H_
#define LANE_LANE_H_

#include <pthread.h>
#include <cctype>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstdio>
#include <cstdlib>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

namespace lane_maxwell {
typedef std::vector<std::FILE*> Files;

inline size_t min(const size_t& a, const size_t& b) {
    return (a < b)? a: b;
}

class Mutex {
public:
    Mutex(): mut(), mattr() {}

    void init() {
        pthread_mutex_init (&mut,&mattr);
    }

    void lock() {
        pthread_mutex_lock (&mut);
    }

    void unlock() {
        pthread_mutex_unlock (&mut);
    }

private:
    pthread_mutex_t mut;
    pthread_mutexattr_t mattr;

    DISALLOW_COPY_AND_ASSIGN(Mutex);
};


void interpolImage(float **image, int nx, int ny);

void levelCurve(
     int im,		/* Limite da Matriz (i)           */
     int jm,		/* Limite da Matriz (j)           */
     float **a,		/* Matriz                         */
     const char* fname,	/* Nome do Arquivo de saida       */
     const char* comment,	/* Texto a ser inserido na imagem */
     int LOG,		/* usar escala log    ?           */
     int PNG,		/* converter em PNG / JPG   ?     */
     int DAT		/* salvar arquivo DAT ?	          */
     );

extern Mutex printMutex;
}  // lane_maxwell

#endif  // LANE_LANE_H_
