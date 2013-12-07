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
@file test_matrix_1.cc

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include <cstdio>
#include "tests/test_matrix_1.h"

using namespace lane_maxwell::modules::matrix;

int main() {
    lane_maxwell::modules::matrix::matrixModuleInit(1e-6);

    size_t N, M;

    printf("(N,M)=");
    scanf("%d%d",&N,&M);

    Matrix* A = new Matrix(N, M);
    Matrix* B = new Matrix(N, M);

    bool status;


    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            A->set(i, j, i*M+i*j+std::sin(j + i * j));
        }
    }

    printf("\nA=\n");

    print(stdout, A);

    B->copyFrom(A);

    Matrix* Ai;
    CallInv* cAi = new CallInv(N);

    cAi->newCall(A);
    inv(cAi, A);

    Ai = cAi->Ai;
    status = cAi->status;

    Matrix* C = new Matrix(N, N);
    mult(C, Ai, B);

    Matrix* D = new Matrix(N, N);
    mult(D, B, Ai);

    printf("\nAi*A=\n");

    print(stdout, C);

    printf("\nA*Ai=\n");

    print(stdout, D);

    Matrix* E = new Matrix(N, N);
    add(E, Ai, B);

    printf("\nAi+A=\n");
    print(stdout, E);

    Matrix* F = new Matrix(N, N);
    mult(F, 0.5, B);

    Matrix* G = new Matrix(N, N);
    mult(G, B, 0.5);

    Matrix* H = new Matrix(N, N);
    trans(H, B);

    if (status) {
        printf("\ninv(A)=\n");
        print(stdout, Ai);
    } else {
        printf("Não existe inversa de A!\n");
    }

    printf("Feito!\n");

    delete A;
    delete B;
    delete C;
    delete D;
    delete E;
    delete F;
    delete G;
    delete H;
    delete cAi;

    return 0;
}

