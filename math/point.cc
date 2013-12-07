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
@file point.cc


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação foi feita por Péricles Lopes Machado

Baseado na implementação  do PROF. DR. RODRIGO M.S. DE OLIVEIRA
*/

#include "math/point.h"
#include "maxwell/maxwell.h"
#include "domains/domains.h"
#include <cstdio>
#include <iostream>

namespace lane_maxwell { namespace modules { namespace point {

void Point::addSD(const Point* const p) {
    m_sd.push_back(p);
    m_df.push_back(new matrix::Vector(m_D));
}

double accSD(const Points* p) {
    double acc = 0;
    Points::const_iterator it = p->begin(), end = p->end();

    while (it != end) {

        for (size_t i = 0; i < (*it)->sizeSD(); ++i) {
            acc += distE(*it, (*it)->getSD(i));
        }

        ++it;
    }

    return acc;
}

void free(Points* p) {
    Points::iterator it = p->begin(), end = p->end();

    while (it != end) {
        delete *it;
        ++it;
    }

    delete p;
}

void print(std::FILE* const f, const Point* const p) {
    fprintf(f, "[ id = %d, \t", p->getID());
    matrix::print(f, p->getVector());
    fprintf(f, "]");
}

void Point::printSD(std::FILE* const f) {
    fprintf(f, "\n{\n");
    for (size_t i = 0; i < m_sd.size(); ++i) {
        print(f, m_sd[i]);
    }
    fprintf(f, "}\n");
}

void printAll(std::FILE* const f, const Point* const p) {
    fprintf(f, "[id = %d, \t", p->getID());
    fprintf(f, "(");
    for (size_t i = 0; i <  p->size(); ++i) {
        fprintf(f, " sigma[%d] = %e, ", i, p->getSigma(i));
    }
    fprintf(f, "eps = %e, mu = %e, prop = %d, type = %c", p->getEps(), p->getMu(),
            p->getProp(), p->getType());
    fprintf(f, "), \t(");
    matrix::print(f, p->getVector());
    fprintf(f, ")]");
}

void print(std::FILE* const f, const Points* const p) {
    Points::const_iterator it = p->begin(), end = p->end();
    int i = 1;

    while (it != end) {
        print(f, *it);
        fprintf(f, " [%d]\n", i);
        ++i;
        ++it;
    }
}

void printAll(std::FILE* const f, const Points* const p) {
    Points::const_iterator it = p->begin(), end = p->end();
    int i = 1;

    while (it != end) {
        printAll(f, *it);
        fprintf(f, " [%d]\n", i);
        ++i;
        ++it;
    }
}


double mult(const Point* const a, const Point* const b) {
    double ans = 0;
    Point::const_iterator it = a->begin(), end = a->end();
    Point::const_iterator itb = b->begin();

    while ( it != end ) {
        ans += (*it) * (*itb);
        ++it;
        ++itb;
    }

    return ans;
}

double norm(const Point* const a) {
    double ans = 0;
    Point::const_iterator it = a->begin(), end = a->end();

    while ( it != end ) {
        ans += (*it) * (*it);
        ++it;
    }

    return ans;
}

void mult(Point* const p, const double* const a, const Point* const b) {
    p->copyFrom(b);

    Point::iterator it = p->begin(), end = p->end();

    while ( it != end ) {
        (*it) *= (*a);
        ++it;
    }
}

void mult(Point* const p, const Point* const a, const double* const b) {
    p->copyFrom(a);

    Point::iterator it = p->begin(), end = p->end();

    while ( it != end ) {
        (*it) *= (*b);
        ++it;
    }
}

void add(Point* const p, const Point* const a, const Point* const b) {
    p->copyFrom(a);

    Point::const_iterator it = b->begin(), end = b->end();
    Point::iterator itp = p->begin();

    while ( it != end ) {
        (*itp) +=  (*it);
        ++it;
        ++itp;
    }
}

void sub(Point* const p, const Point* a, const Point* b) {
    p->copyFrom(a);

    Point::const_iterator it = b->begin(), end = b->end();
    Point::iterator itp = p->begin();

    while ( it != end ) {
        (*itp) -=  (*it);
        ++it;
        ++itp;
    }
}

void dot(Point* const p, const Point* a, const Point* b) {
    p->copyFrom(a);

    p->setX(0, a->getX(1)*b->getX(3)-a->getX(3)*b->getX(1));
    p->setX(1, a->getX(0)*b->getX(3)-a->getX(3)*b->getX(0));
    p->setX(2, a->getX(0)*b->getX(1)-a->getX(1)*b->getX(0));
}

bool operator<(const Point& a, const Point& b) {
    Point::const_iterator ita = a.begin(), end = a.end();
    Point::const_iterator itb = b.begin();

    while ( ita != end ) {
        int c = matrix::cmpD(*ita, *itb);

        if (c < 0) {
            return true;
        }

        if (c > 0) {
            return false;
        }

        ++ita;
        ++itb;
    }

    return false;
}

double distE(const Point* const a, const Point* const b) {
    double d = 0;
    Point::const_iterator it = a->begin(), end = a->end();
    Point::const_iterator itp = b->begin();

    while ( it != end ) {
        double de = *it - *itp;

        d += de * de;

        ++it;
        ++itp;
    }

    return sqrt(d);
}

double distMax(const Point* const a, const Point* const b) {
    Point::const_iterator it = a->begin(), end = a->end();
    Point::const_iterator itp = b->begin();
    double d = fabs(*it-*itp);

    ++it;
    ++itp;

    while ( it != end ) {
        double de = fabs(*it - *itp);

        if ( de > d ) {
            d = de;
        }

        ++it;
        ++itp;
    }

    return d;
}

double distM(const Point* const a, const Point* const b) {
    Point::const_iterator it = a->begin(), end = a->end();
    Point::const_iterator itp = b->begin();
    double d = 0;

    while ( it != end ) {
        d += fabs(*it - *itp);

        ++it;
        ++itp;
    }

    return d;
}

bool save(const char* const fname, const Points* const p) {
    std::FILE* F = fopen(fname, "wb+");

    if (!F) {
        return false;
    }

    size_t N = p->size();
    int I = 0;
    Points::const_iterator it = p->begin(), end = p->end();

    size_t dim;

    dim = (*it)->getDim();

    rewind(F);
    fwrite(&N, sizeof(N), 1, F);
    fwrite(&dim, sizeof(dim), 1, F);

    while (it != end) {
        char type = (*it)->getType();
        double X;
        double eps = (*it)->getEps();
        double mu = (*it)->getMu();
        int prop = (*it)->getProp();

        I = (*it)->getID();

        fwrite(&type, sizeof(type), 1, F);
        fwrite(&I, sizeof(I), 1, F);

        for (size_t i = 0; i < dim; ++i) {
            X = (*it)->getX(i);
            fwrite(&X, sizeof(X), 1, F);
        }

        fwrite(&eps, sizeof(eps), 1, F);
        fwrite(&mu, sizeof(mu), 1, F);
        fwrite(&prop, sizeof(prop), 1, F);

        for (size_t i = 0; i < dim; ++i) {
            X = (*it)->getSigma(i);
            fwrite(&X, sizeof(X), 1, F);
        }

        ++it;
    }

    fclose(F);
    return true;
}

/* Adaptação da implementação do professor Rodrigo oliveira */
void coloumb2D(Points* p, double dx, int NP, double K, double Rmax, double Rmin) {

    size_t N = p->size();

    std::vector<double> x(N), y(N), xo(N), yo(N);
    std::vector<bool> f(N, true);

    int I;


    for (int i = 0; i < N; ++i) {
        int prop = p->at(i)->getProp();

        xo[i] = x[i] = p->at(i)->getX(0);
        yo[i] = y[i] = p->at(i)->getX(1);

        if (prop != 0) {
            f[i] = false;
        }
    }


    printf("Gerando dominios de suporte...\n");
    double rmin = domains::buildSupportDomain(stdout, 'a', 'a', 20, p, 0.01, distE, 1, true, Rmax);

    fprintf(stdout, "rmin = %e\n", rmin);

    printf("Atualizando a posicao dos pontos..\n");

    double ix, iy, R;

    for (int k = 0; k <= NP; ++k) {

        if ((k%10) == 0) printf("%d / %d\n", k, NP);

        for (int i = 0; i < N; ++i) { //ponto atual
            Point* pt = (*p)[i];
            int I = pt->getID();

            if (f[I]) { //se os pontos forem moveis

                ix = 0.0e0; //contribuicoes de cada ponto no ponto i
                iy = 0.0e0;

                int M = pt->sizeSD();

                for (int j = 0; j < M ; ++j) {
                    const Point* pivot = pt->getSD(j);
                    if( pivot->getID() != pt->getID()) {

                        int ip = pivot->getID();

                        R = sqrt(pow(xo[I]-xo[ip], 2.0e0) + pow(yo[I]-yo[ip], 2.0e0));
                        R = pow(R, 3.0e0);

                        if (R > Rmin){ //evita explosoes de forca
                            ix += (xo[i]-xo[ip]) / R;  // Lei de Coulomb (soma vetorial)
                            iy += (yo[i]-yo[ip]) / R;
                        }
                    }
                }

                R = sqrt(ix * ix + iy * iy);

                ix = ix / R;
                iy = iy / R;

                x[I] += ix * dx / K; // isto garante que os pontos não se movem muito rapidamente
                y[I] += iy * dx / K;

                pt->setX(0, x[I]);
                pt->setX(1, y[I]);
            }
        }


        for (int j = 0; j < N; ++j){
            xo[j] = x[j];
            yo[j] = y[j];
        }
    }
}

bool load(const char* const fname, Points** p) {
    std::FILE* F = fopen(fname, "rb");

    if (!F) {
        return false;
    }

    size_t N, dim;

    rewind(F);

    fread(&N, sizeof(N), 1, F);
    fread(&dim, sizeof(dim), 1, F);

    size_t I = 0;

    *p = new Points;

    while (I < N) {
        char type;

        double X;
        double eps;
        double mu;

        int ID;
        int prop;

        Point* np = new Point;
        matrix::Vector* v = new matrix::Vector(dim);

        fread(&type, sizeof(type), 1, F);
        fread(&ID, sizeof(ID), 1, F);

        for (size_t i = 0; i < dim; ++i) {
            fread(&X, sizeof(X), 1, F);
            v->set(i, X);
        }

        fread(&eps, sizeof(eps), 1, F);
        fread(&mu, sizeof(mu), 1,  F);
        fread(&prop, sizeof(prop), 1, F);


        np->set(ID, eps, mu, prop, type, v);

        for (size_t i = 0; i < dim; ++i) {
            fread(&X, sizeof(X), 1, F);
            np->setSigma(i, X);
        }

        (*p)->push_back(np);

        ++I;
    }

    fclose(F);
    return true;
}
}}}  // lane_maxwell::modules::point
