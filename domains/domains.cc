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
@file domains.cc

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação foi feita por Péricles Lopes Machado

*/

#include "domains/domains.h"

namespace lane_maxwell { namespace modules { namespace domains {

bool operator<(const Double& a, const Double& b ) {
    return matrix::cmpD(a.v, b.v) < 0;
}


void printBox(std::FILE* const f, const Box* const box, const int l,
                const int M, const double& delta) {
    if ( l == M-1 ) {
        for (int i = 0; i < l; ++i) {
             std::fprintf(f, "\t");
         }

        SetPoint::const_iterator it = box->w->begin(), end = box->w->end();

        while (it != end) {
            fprintf(f, "[");
            fprintf(f, " %p ", *it);
            point::printAll(f, *it);
            fprintf(f, " ]");
            ++it;
        }

        fprintf(f, "\n___________________\n");
    } else {
        MapBox::const_iterator it = box->next.begin(), end = box->next.end();

        while (it != end) {
            Box* b = it->second;

            for (int i = 0; i < l; ++i) {
                std::fprintf(f, "\t");
            }

            fprintf(f, " [ %e, %e]:\n", it->first * delta, (it->first + 1) * delta);
            printBox(f, b, l + 1, M, delta);

            ++it;
        }
    }
}

void insertElement(Box* const box, const double& w, const int l,
                    const point::Point* const x, const int M,
                    const double& delta) {
    int v = static_cast<int>(round(w/delta));

    if (l == M - 1) {
        box->w->insert(x);
    } else {
        Box* nB;
        MapBox::iterator it = box->next.find(v);

        if (it == box->next.end()) {
             nB = box->next[v] = new Box();

             if (l == M-2) {
                 nB->setPoint(new SetPoint);
             }
        } else {
            nB = it->second;
        }

        if (l <= M - 3) {
            insertElement(nB, x->getX(l + 2), l+1, x, M, delta);
        } else {
            insertElement(nB, 0, l+1, x, M, delta);
        }
    }
}

void removeElement(Box* const box, const double& w, const int l,
                    const point::Point* const x, const int M,
                    const double& delta) {
    int v = static_cast<int>(round(w/delta));

    if (l == M-1) {
        box->w->erase(x);
    } else {
        Box* nB;
        MapBox::iterator it = box->next.find(v);

        if (it != box->next.end()) {
            nB = it->second;

            if (l <= M-3) {
                removeElement(nB, x->getX(l+2), l+1, x, M, delta);
            } else {
                removeElement(nB, 0, l+1, x, M, delta);
            }

            if (nB != 0) {
                if (nB->next.size() == 0 && l+1 < M-1) {
                    delete nB;
                    nB = 0;
                    box->next.erase(v);
                } else if (l == M-2) {
                    if (nB->w->size() == 0) {
                        delete nB;
                        nB = 0;
                        box->next.erase(v);
                    }
                }
            }
        }
    }
}

void buildCube(Box* const box, const int l, const point::Point* const x,
                const double& h, Cube* cube, const int M,
                const double& delta, FPDist d, MapPoint& marc,
                const int p) {
    if ( l == M - 1 ) {
        SetPoint::iterator it = box->w->begin(), end = box->w->end();

        while (it != end) {
            if (x != *it) {
                double dist = d(x, *it);

                if ( marc[*it] != p ) {
                    marc[*it] = p;
                    cube->push_back(PairDoublePoint(dist, *it));
                }
            }
            ++it;
        }
    } else {
        int ini = static_cast<int>(round((x->getX(l+1) - h - delta) / delta));
        int fim = static_cast<int>(round((x->getX(l+1) + h + delta) / delta));

        std::map<int, Box*>::iterator bg = box->next.lower_bound(ini);
        std::map<int, Box*>::iterator ed = box->next.lower_bound(fim);
        std::map<int, Box*>::iterator it = bg;

        while (it != ed) {
            buildCube(it->second, l + 1, x, h, cube, M, delta, d, marc, p);
            ++it;
        }
    }
}
}}}  // lane_maxwell::modules::domains
