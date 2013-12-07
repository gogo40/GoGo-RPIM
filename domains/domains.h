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
@file domains.h
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#ifndef DOMAINS_DOMAINS_H_
#define DOMAINS_DOMAINS_H_

#include <map>
#include <set>
#include <vector>
#include <utility>

#include "math/point.h"

namespace lane_maxwell { namespace modules { namespace domains {

struct Double {
    double v;
    explicit Double(const double& v = 0):v(v) {}
};

struct Box;

typedef std::pair<double, const point::Point*> PairDoublePoint;
typedef std::set<int> SetInt;
typedef std::set<const point::Point*> SetPoint;
typedef std::map<int, Box*> MapBox;
typedef std::map<const point::Point*, int> MapPoint;

struct Box {
        explicit Box()
            : w(0) {}

        void setPoint(SetPoint* nw) {
            w = nw;
        }

        ~Box() {
            if (w) {
                w->clear();
                delete w;
            } else {
                MapBox::iterator it = next.begin(), end = next.end();

                while (it != end) {
                    if(it->second) {
                        delete it->second;
                    }
                    ++it;
                }
            }
        }
        MapBox next;
        SetPoint* w;
};

typedef double (*FPDist) (const point::Point* const, const point::Point* const);
typedef std::vector<PairDoublePoint> Cube;

bool operator<(const Double& a, const Double& b);

void insertElement(Box* const box, const double& v, const int l,
                const point::Point* const x, const int M,
                const double& delta);
void removeElement(Box* const box, const double& v, const int l,
                const point::Point* const x, const int M,
                const double& delta);
void buildCube(Box* const box, const int l, const point::Point* const x,
                const double& h, Cube* cube, const int M, const double& delta,
                FPDist d, MapPoint& marc, const int p);
void printBox(std::FILE* const f, const Box* const box, const int l,
                const int M, const double& delta);
}}}  // lane_maxwell::modules::domains

#endif  // DOMAINS_DOMAINS_H_

