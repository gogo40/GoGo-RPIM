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
@file supportdomain.cc

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação foi feita por Péricles Lopes Machado

*/

#include <algorithm>

#include "domains/domains.h"

namespace lane_maxwell { namespace modules { namespace domains {

double buildSupportDomain(std::FILE* f, char B, char A, size_t T,
    const point::Points* P, double fator, FPDist d,
    double fator_osc, bool useR, double rd) {
    int N = P->size();
    size_t M = (*P)[0]->size();

    int ip;

    double h, ho, rmin;
    double den, vol, delta;

    matrix::Vector* xmin = new matrix::Vector(M);
    matrix::Vector* xmax = new matrix::Vector(M);

    MapPoint marc;

    point::Points::const_iterator it = P->begin(), end = P->end();

    while (it != end) {
        if ((*it)->getType() == A ||  A == 'a') {
            for (size_t i = 0; i < M; ++i) {
                xmin->set(i, (*it)->getX(i));
                xmax->set(i, (*it)->getX(i));
            }
            break;
        }

        ++it;
    }

    it = P->begin();
    end = P->end();

    while (it != end) {
        if ((*it)->getType() == A ||  A == 'a') {
            marc[*it]=-1;
            for (size_t i = 0; i < M; ++i) {
                if ( xmin->get(i) > (*it)->getX(i) ) {
                    xmin->set(i, (*it)->getX(i));
                }

                if (xmax->get(i) < (*it)->getX(i)) {
                    xmax->set(i, (*it)->getX(i));
                }
            }
        }

        ++it;
    }

    vol = 1;

    for(size_t j = 0; j < M; ++j) {
        vol *= (xmax->get(j) - xmin->get(j));
    }

    delete xmin;
    delete xmax;

    fprintf(f, "\nvol = %e\nfator_osc = %e\n", vol, fator_osc);

    Box* box = new Box;

    den = N / vol;
    ho = std::pow(T / den, 1.0 / M);
    delta = fator_osc * std::pow(2.0 / den, 1.0 / M);

    ip = 0;
    rmin = -1;

    point::Points::const_iterator begin = P->begin();
    point::Points::const_iterator it_p, left, right;

    it_p = left = right = begin;

    fprintf(f, "Gerando domínio de suporte...\n");
    fprintf(f, "delta=%e\n", delta);

    while ( it_p != end ) {
        point::Point* p = *it_p;

        if (p->getType() == B ||  A == 'a') {
            h = ho;
            if (useR) h = rd;

            Cube* cube = new Cube;

            while (cube->size() <= T) {

                while (matrix::cmpD(p->getX(0) - (*left)->getX(0) , h+delta) > 0) {
                    if (left == it_p) break;

                    if( (*left)->getType() == A  ||  A == 'a') {
                        removeElement(box, (*left)->getX(1), 0, *left, M, delta);
                    }

                    ++left;
                }

                while (matrix::cmpD((*right)->getX(0) - p->getX(0), h+delta) > 0) {
                    if (right == begin || right == it_p) break;

                    if ((*right)->getType() == A  ||  A == 'a') {
                        removeElement(box, (*right)->getX(1), 0 , *right, M, delta);
                    }

                    --right;
                }

                while (matrix::cmpD((*right)->getX(0) - p->getX(0), h+delta) <= 0) {
                    if ( (*right)->getType() == A  ||  A == 'a') {
                        insertElement(box, (*right)->getX(1), 0, *right, M, delta);
                    }

                    ++right;

                    if (right == end) {
                        --right;
                        break;
                    }
                }

                while (matrix::cmpD(p->getX(0) - (*left)->getX(0), h+delta) <= 0) {
                    if( (*left)->getType() == A ||  A == 'a') {
                        insertElement(box, (*left)->getX(1), 0, *left, M, delta);
                    }

                    if(left == begin) break;

                    left--;
                }

                buildCube(box, 0, *it_p, h, cube, M, delta, d, marc, ip);

                h = (1 + fator) * h;
                if (useR) break;
            }

            p->setRMax(0);
            std::sort(cube->begin(), cube->end());
            Cube::iterator it = cube->begin(), cend = cube->end();

            size_t i = 0;

            if (useR) {
                while (it != cend) {

                    if (it->first > rd) break;

                    p->addSD(it->second);

                    if (p->getRMax() < it->first) {
                        p->setRMax(it->first);
                    }

                    if (it->first < rmin || rmin < 0) {
                        rmin = it->first;
                    }
                    i++;
                    it++;
                }
            } else {
                while (i < T) {
                    p->addSD(it->second);

                    if (p->getRMax() < it->first) {
                        p->setRMax(it->first);
                    }

                    if (it->first < rmin || rmin < 0) {
                        rmin = it->first;
                    }

                    it++;
                    i++;
                }
            }

            ho = p->getRMax();

            delete cube;

            if (ip % 10000 == 0) {
                fprintf(f, "< %d / %d > %e %d\n", ip, N, ho, i);
            }
        }


        ++ip;
        ++it_p;
    }

    delete box;

    fprintf(f, "Domínio de suporte gerado!\n");

    return rmin;
}
}}}  // lane_maxwell::modules::domains

