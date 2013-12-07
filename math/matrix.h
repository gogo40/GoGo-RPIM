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
@file matrix.h


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#ifndef MATH_MATRIX_H_
#define MATH_MATRIX_H_

#include <vector>

#include "lane/lane.h"

namespace lane_maxwell { namespace modules { namespace matrix {

/**
@class Vector
@brief Um vetor de doubles
*/

class Vector {
 public:
    typedef double* iterator;
    typedef const double* const_iterator;

/**
@brief Construtor da classe Vector.
@param size é o numero de linhas do vetor, o valor default de size é 1.
*/
    explicit Vector(size_t size)
        : m_size(size) {
        m_vector = static_cast<double*> (std::malloc(sizeof(double)*size));
    }

    void init(const double& v) {
        iterator it = begin(), end = this->end();

        while (it != end) {
            *it = v;
            ++it;
        }
    }

    void set(size_t i, const double& v) {
        m_vector[i] = v;
    }

    double& get(size_t i) const {
        return m_vector[i];
    }

    double& operator[](size_t i) {
        return m_vector[i];
    }

    iterator begin() const {
        return m_vector;
    }

    iterator end() const {
        return m_vector + m_size;
    }

    size_t size() const {
        return m_size;
    }

    void copyFrom(const Vector* const v) {
        iterator it = begin();
        const_iterator itv = v->begin(), end = v->end();

        while (itv != end) {
            *it = *itv;
            ++it;
            ++itv;
        }
    }

    ~Vector() {
        std::free(m_vector);
    }

 private:
    size_t m_size;
    double* m_vector;

    DISALLOW_COPY_AND_ASSIGN(Vector);
};

/**
@class Matrix
@brief Matrix é uma classe que define uma matriz de doubles.
*/

class Matrix {
 public:
    typedef Vector** iterator;
    typedef Vector** const_iterator;

/**
@brief Construtor da classe Matrix.
@param n é o numero de linhas da matriz, o valor default de n é 1.
@param v é o valor que deve ser usado para inicializar a matriz, o valor default de v é 0.
*/
    explicit Matrix(size_t n)
        : m_n(n), m_m(0) {
        m_matrix = new Vector*[n];
    }

/**
@brief Construtor da classe Matrix.
@param n é o numero de linhas da matriz, o valor default de n é 1.
@param m é o numero de colunas da matriz, o valor default de m é 1.
@param v é o valor que deve ser usado para inicializar a matriz, o valor default de v é 0.
*/
    Matrix(size_t n, size_t m)
        : m_n(n), m_m(m) {
        m_matrix = static_cast<Vector**>(malloc(sizeof(Vector*) * n));

        for (size_t i = 0; i < n ; ++i) {
            m_matrix[i] = new Vector(m);
        }
    }

    void copyFrom(const Matrix* const m) {
        const_iterator it = begin(), end = this->end();
        const_iterator jt = m->begin();

        while (it != end) {
            (*it)->copyFrom(*jt);
            ++it;
            ++jt;
        }
    }

    size_t getN() const {
        return m_n;
    }

    size_t getM() const {
        return m_m;
    }

    void init(const double& v) {
        iterator it = begin(), end = this->end();

        while (it != end) {
            (*it)->init(v);
            ++it;
        }
    }

    Vector& operator[](size_t i) {
        return *m_matrix[i];
    }

    void setVector(size_t i, Vector* v) {
        m_matrix[i] = v;
    }

    void set(size_t i, size_t j, const double& v) {
        m_matrix[i]->set(j, v);
    }

    double& get(size_t i, size_t j) const {
        return m_matrix[i]->get(j);
    }

    iterator begin() const {
        return m_matrix;
    }

    iterator end() const {
        return m_matrix+m_n;
    }

    ~Matrix() {
        for (size_t i = 0; i < m_n; ++i) {
            delete m_matrix[i];
        }
        free(m_matrix);
    }

 private:
    size_t m_n, m_m;
    Vector** m_matrix;

    DISALLOW_COPY_AND_ASSIGN(Matrix);
};

void matrixModuleInit(const double& tol = 1e-12);

void print(std::FILE* const f, const Vector* const v);

void print(std::FILE* const f, const Matrix* const v);

/**
@fn void matrix(Matrix* res, const Vector* x)
@brief Função que converte um vetor Vector em uma matriz Matrix.
@see Matrix
@see Vector
@param x é vetor  com n elementos que será convertido.
*/
void matrix(Matrix* const res, const Vector* const x);

/**
@fn void trans(Matrix* const res, const Matrix* const A)
@brief Função que gera a transposta de uma matriz.
@see Matrix
@param A é matriz que será transposta.
*/
void trans(Matrix* const res, const Matrix* const A);

/**
@fn void inv(Matrix* const res, const Matrix* const A, bool& status)
@brief Função que gera a inversa de uma matriz
@see Matrix
@see gauss
@param A é matriz a ser invertida.

A matriz A deve ser quadrada e inversível, caso contrário, uma mensagem de erro
será gerada e o programa será abortado (a função exit(1) será chamada). Esta função
utiliza eliminação de Gauss com pivoteamento parcial para calcular a matriz inversa de A.
*/

typedef std::vector<size_t> VectorInt;

struct CallInv {
    explicit CallInv(size_t N) {
        n = N;
        Ai = new Matrix(n, n);
        A =  new Matrix(n, n);
        I = new Matrix(n, n);
        c = new Matrix(n, n);
        y = new Matrix(n, n);
        p = new VectorInt(n);
    }

    void newCall(const Matrix* const nA) {
        A->copyFrom(nA);
    }

    void clear() {
        delete Ai;
        delete A;
        delete I;
        delete c;
        delete y;
        delete p;
    }

    ~CallInv() {
        delete Ai;
        delete A;
        delete I;
        delete c;
        delete y;
        delete p;
    }

    size_t n;
    Matrix* Ai;
    Matrix* A;
    Matrix* I;
    Matrix* c;
    Matrix* y;
    VectorInt* p;
    bool status;
};

void inv(Matrix*  Ai, Matrix* A, Matrix*  I,
            Matrix*  c, Matrix*  y, VectorInt*  p,
            bool*  status);

void inv(CallInv* const cA, Matrix* const A);

/**
@fn Matrix* mult(const Matrix*& A, const Matrix*& B)
@brief Função que realiza o produto de duas matrizes
@see Matrix
@param A é uma matriz Matrix com n linhas e m colunas.
@param B é uma matriz Matrix com m linhas e p colunas.
@return Uma matriz Matrix com n linhas e p colunas contendo o produto de A e B, ou seja, C=A*B.

m tem de ser igual a k para que o produto seja realizado, caso contrário, poderá haver erro
de tempo de execução.
*/
void mult(Matrix* const res, const Matrix* const A, const Matrix* B);

/**
@fn Matrix* mult(const double& A, const Matrix* B)
@brief função que realiza o produto de um escalar por uma matriz
@see Matrix
@param A é um escalar.
@param B é uma matriz Matrix com k linhas e p colunas.
@return Uma matriz Matrix com k linhas e p colunas contendo o produto de A e B, ou seja, C=A*B.
*/
void mult(Matrix* const res, const double& A, const Matrix* B);

/**
@fn Matrix* mult(const Matrix* B, const double& A)
@brief função que realiza o produto de um escalar por uma matriz
@see Matrix
@param B é uma matriz Matrix com k linhas e p colunas.
@param A é um escalar.
@return Uma matriz Matrix com k linhas e p colunas contendo o produto de A e B, ou seja, C=A*B.
*/

void mult(Matrix* const res, const Matrix* B, const double& A);

/**
@fn Matrix* add(const Matrix* A, const Matrix* B)
@brief função que realiza a soma de duas matrizes
@see Matrix
@param A é uma matriz Matrix com n linhas e m colunas.
@param B é uma matriz Matrix com k linhas e p colunas.
@return Uma matriz Matrix com n linhas e m colunas contendo a soma de A e B, ou seja, C=A+B.

m tem de ser igual a p e n tem de ser igual a k  para que a soma seja realizada, caso contrário,
poderá haver erro de tempo de execução.
*/
void add(Matrix* const res, const Matrix* const A, const Matrix* const B);

/**
@fn Matrix sub(const Matrix* A, const Matrix* B)
@brief função que realiza a soma de duas matrizes
@see Matrix
@param A é uma matriz Matrix com n linhas e m colunas.
@param B é uma matriz Matrix com k linhas e p colunas.
@return Uma matriz Matrix com n linhas e m colunas contendo a subtração de A e B, ou seja, C=A-B.

m tem de ser igual a p e n tem de ser igual a k  para que a subtração seja realizada,
caso contrário, poderá haver erro de tempo de execução.
*/
void sub(Matrix* const res, const Matrix* const A, const Matrix* const B);

/**
@fn inline void setTol(double eps=1e-20)
@brief Função usada para alterar a tolerância EPS usadas em algumas operações de ponto flutuante.
*/
void setTol(const double& eps);

/**
@fn  inline int cmpD(double x, double y, double eps=EPS)
@brief Função usada para comparar dois doubles com uma tolerância de EPS.
@see EPS
@param x é o primeiro valor a ser comparado.
@param y é o segundo valor a ser comparado.
@param eps é a tolerância usada na comparação.
@return 0 caso \f$  \left| x-y \right|  < eps \f$ , 1 caso \f$ x>y+eps \f$ e
retorna -1 caso \f$ x<y-eps \f$ .
*/
int cmpD(const double& x, const double& y);

/**
@fn void gauss(int n, Matrix* A, Matrix* X, Matrix* B)
@brief Função que resolve a equação matricial linear AX=B, caso exista solução.
@see Matrix
@param n é a dimensão do problema, ou seja, o numero de linhas da matriz quadrada a.
@param A é uma matriz com pelo menos n linhas e pelo menos n colunas, caso a matriz a
tenha menos de n linhas e colunas haverá erro de tempo de execução.
@param X é uma matriz icognita com pelo menos n linhas e pelo menos p colunas, caso a
matriz X tenha menos de n linhas e p colunas haverá erro de tempo de execução. Nessa
matriz será jogado o resultado da equação AX=B.
@param B é uma matiz com pelo menos n linhas e p colunas, caso a matriz a tenha menos de
n linhas haverá erro de tempo de execução.

Esta função soluciona a equação matricial AX=B, ou seja, retorna \f$ X = A^{-1}B \f$  .
Para resolver a equação, esta função implementa a eliminação de Gauss com uma estrátegia
de pivoteamento parcial para reduzir erros de arrendondamento.
*/
void gauss(size_t n, Matrix* const A, Matrix* const X, Matrix* const B,
            Matrix* const c, Matrix* const y, VectorInt* const p,
            bool* const st);
/**
@fn double dot(const Matrix* A, const Matrix* B)
@brief Função que computa o produto interno de duas matrizes A e B
@see Matrix
@param A é uma matriz quadrada NxN
@param B é uma matriz quadrada NxN
@return o produto interno de A e B, ou seja, \f$ \sum_{i=0}^{N} { \sum_{j=0}^{N}{A_{ij} B_{ij}}}
\f$, onde N é o tamanho da matriz.
*/
double dot(const Matrix* const A, const Matrix* const B);
}}}  // lane_maxwell::modules::matrix

#endif  // MATH_MATRIX_H_
