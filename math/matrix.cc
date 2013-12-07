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
@file matrix.cc


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação foi feita por Péricles Lopes Machado

Baseado na implementação  do PROF. DR. RODRIGO M.S. DE OLIVEIRA
*/

#include "math/matrix.h"
#include "math/point.h"

namespace lane_maxwell { namespace modules { namespace matrix {

void print(std::FILE* const f, const Vector* const v) {
    Vector::const_iterator it = v->begin(), end = v->end();

    while (it != end) {
        fprintf(f, "%e\t", *it);
        ++it;
    }
}

void print(std::FILE* const f, const Matrix* const v) {
    Matrix::const_iterator it = v->begin(), end = v->end();

    while (it != end) {
        print(f, *it);
        fprintf(f, "\n");
        ++it;
    }
}

void matrix(Matrix* const A, const Vector* const x) {
    Vector::const_iterator it = x->begin(), end = x->end();

    size_t i = 0;

    while ( it!= end ) {
        A->set(0, i, *it);
        ++it;
        ++i;
    }
}

void trans(Matrix* const C, const Matrix* const A) {
    size_t nA = A->getN(), mA = A->getM();

    for (size_t i = 0; i < nA; ++i) {
        for (size_t j = 0; j < mA; ++j) {
            C->set(j, i, A->get(i, j));
        }
    }
}

/**
 * Warning: The matrix A must be a copy! Gauss change values!
 */

void inv(CallInv* const cA, Matrix* const A) {
    cA->newCall(A);
    inv(cA->Ai, cA->A, cA->I, cA->c, cA->y, cA->p, &cA->status);
}

void inv(Matrix*  Ai, Matrix* A, Matrix* I,
            Matrix*  c, Matrix* y, VectorInt*  p,
            bool*  status) {
    size_t n = A->getN();

    I->init(0);

    for (size_t i = 0; i < n; ++i) {
        I->set(i, i, 1.0);
    }

    gauss(n, A, Ai, I, c, y, p,  status);
}

void mult(Matrix* const C, const Matrix* const A, const Matrix* const B) {
    size_t nA = A->getN(), mA = A->getM();
    size_t mB = B->getM();

    C->init(0);

    for (size_t i = 0; i < nA; ++i) {
        for (size_t j = 0; j < mB; ++j) {
            for (size_t k = 0; k < mA; ++k) {
                C->set(i, j, C->get(i, j) + A->get(i, k) * B->get(k, j));
            }
        }
    }
}

void mult(Matrix* const C, const double& A, const Matrix* const B) {
    size_t nC = B->getN(), mC = B->getM();

    for(size_t i = 0; i < nC; ++i) {
        for(size_t j = 0; j < mC; ++j) {
            C->set(i, j, B->get(i, j) * A);
        }
    }
}

void mult(Matrix* const C, const Matrix* B, const double& A) {
    size_t nC = B->getN(), mC = B->getM();

    for(size_t i = 0; i < nC; ++i) {
        for(size_t j = 0; j < mC; ++j) {
            C->set(i, j, B->get(i, j) * A);
        }
    }
}

void add(Matrix* const C, const Matrix* A, const Matrix* B) {
    size_t nA = A->getN(), mA = A->getM();

    for (size_t i = 0; i < nA; ++i) {
        for (size_t j = 0; j < mA; ++j) {
            C->set(i, j, A->get(i, j) + B->get(i, j));
        }
    }
}

void sub(Matrix* const C, const Matrix* A, const Matrix* B) {
    size_t nA = A->getN(), mA = A->getM();

    for (size_t i = 0; i < nA; ++i) {
        for (size_t j = 0; j < mA; ++j) {
            C->set(i, j, A->get(i, j) - B->get(i, j));
        }
    }
}

double EPS = 1e-12;

void setTol(const double& eps) {
    EPS = eps;
}

void matrixModuleInit(const double& tol) {
    setTol(tol);
}

int cmpD(const double& x, const double& y) {
    if (fabs(x-y) < EPS) {
        return 0;
    }

    if (x > y + EPS) {
        return 1;
    }

    return -1;
}

void gauss(size_t n, Matrix* const A, Matrix* const X, Matrix* const B,
            Matrix* const c, Matrix* const y, VectorInt* const p,
            bool* const st) {
    double pv, t, soma;
    size_t r, h;

    *st = true;

    for (size_t i = 0; i < n; ++i) {
        (*p)[i] = i;
    }

    for (size_t k = 0; k < n - 1; ++k) {
        pv = std::fabs(A->get(k, k));
        r = k;

        for (size_t i = k + 1; i < n; ++i) {
            if (cmpD(fabs(A->get(i, k)), pv) > 0) {
                pv = std::fabs(A->get(i, k));
                r = i;
            }
        }

        if (cmpD(pv, 0.0) == 0) {
            *st = false;
            return;
        }

        if (r != k) {
            h = (*p)[k];
            (*p)[k] = (*p)[r];
            (*p)[r] = h;

            for (size_t j = 0; j < n; ++j) {
                t = A->get(k, j);
                A->set(k, j, A->get(r, j));
                A->set(r, j, t);
            }
        }

        for (size_t i = k + 1; i < n; ++i) {
            t = A->get(i, k) / A->get(k, k);
            A->set(i, k, t);

            for(size_t j = k+1; j < n; ++j) {
                 A->set(i, j, A->get(i, j) - t * A->get(k, j));
             }
        }
    }

    pv = std::fabs(A->get(n - 1, n - 1));

    if (cmpD(pv, 0.0) == 0) {
        *st = false;
        return;
    }

    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            r = (*p)[i];
            c->set(i, k, B->get(r, k));
        }

        for (size_t i = 0; i < n; ++i) {
            soma = 0;

            for (size_t j = 0; j < i; ++j) {
                soma += A->get(i, j) * y->get(j, k);
            }

            y->set(i, k, c->get(i, k) - soma);
        }

        for (size_t i = n-1; i < n; --i) {
            soma = 0;

            for (size_t j = i+1; j < n; ++j) {
                 soma += A->get(i, j) * X->get(j, k);
             }

            X->set(i, k, (y->get(i, k) - soma) / A->get(i, i));
        }
    }
}

double dot(const Matrix* const A, const Matrix* const B) {
    double ans = 0;
    size_t N = A->getN();

    for (size_t i = 0; i < N; ++i) {
         for (size_t j = 0; j < N; ++j) {
             ans+=A->get(i, j)*B->get(i, j);
         }
     }
    return ans;
}
}}}  // lane_maxwell::modules::matrix
