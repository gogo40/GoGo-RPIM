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
@file point.h


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#ifndef MATH_POINT_H_
#define MATH_POINT_H_

#include <map>
#include <set>
#include <vector>

#include "math/matrix.h"

namespace lane_maxwell { namespace modules { namespace point {

class Point;

typedef double (*FPDist) (const point::Point* const, const point::Point* const);
double distE(const Point* const a, const Point* const b);
double distMax(const Point* const a, const Point* const b);
double distM(const Point* const a, const Point* const b);

class Point {
 public:
    typedef double* iterator;
    typedef const double* const_iterator;

    explicit Point()
        : m_sd(), m_df(), m_f(), m_sigma(0), m_p(0), m_mirror(), m_fpoint(),
          m_mutexMirror() { }

    void set(int ID, double eps, double mu, int prop, char type,
                matrix::Vector* const p) {
        m_ID = ID;
        m_D = p->size();
        m_sigma = static_cast<double*>(malloc(sizeof(m_eps)*m_D));

        for (size_t i = 0 ; i < m_D; ++i) {
            m_sigma[i] = 0;
        }

        m_eps = eps;
        m_mu = mu;
        m_prop = prop;
        m_type = type;
        m_p = p;

    }

    void copyFrom(const Point* const p) {
        m_p = new matrix::Vector(p->size());
        m_D = p->size();
        m_sigma = static_cast<double*>(malloc(sizeof(m_eps)*m_D));

        if(p->m_sigma != 0) {
            for (size_t i = 0 ; i < m_D; ++i) {
                m_sigma[i] = p->m_sigma[i];
            }
        }

        m_eps = p->m_eps;
        m_mu = p->m_mu;
        m_prop = p->m_prop;
        m_type = p->m_type;
        m_p->copyFrom(p->m_p);
    }

    void setSigma(size_t i, const double& v) {
        m_sigma[i] = v;
    }

    double& getSigma(size_t i) const {
        return m_sigma[i];
    }

    void setEps(const double& eps)  {
        m_eps = eps;
    }

    double getEps() const {
        return m_eps;
    }

    void setMu(const double& mu)  {
        m_mu = mu;
    }

    double getMu() const  {
        return m_mu;
    }

    int getID() const {
        return m_ID;
    }

    void setID(int ID) {
        m_ID = ID;
    }

    size_t size() const {
        return m_p->size();
    }

    void setProp(int prop) {
        m_prop = prop;
    }

    int getProp() const {
        return m_prop;
    }

    void setType(char type) {
        m_type = type;
    }

    char getType() const {
        return m_type;
    }

    void setX(size_t i, const double& v) {
        m_p->set(i, v);
    }

    double& getX(size_t i) const {
        return m_p->get(i);
    }

    double& operator[](size_t i)  {
        return m_p->get(i);
    }

    matrix::Vector* getVector() {
        return m_p;
    }

    matrix::Vector* getVector() const {
        return m_p;
    }

    iterator begin() const {
        return m_p->begin();
    }

    iterator end() const {
        return m_p->end();
    }

    void addSD(const Point* const p);

    size_t sizeSD() const {
        return m_sd.size();
    }

    void printSD(std::FILE* const f);

    const Point* getSD(size_t i) const {
        return m_sd[i];
    }

    void setDF(size_t i, size_t j, const double& v) {
        m_df[i]->set(j, v);
    }

    void setF(size_t i, const double & v) {
        m_f[i] =  v;
    }

    double& getDF(size_t i, size_t j) const {
        return m_df[i]->get(j);
    }

    double getF(size_t i) {
        return m_f[i];
    }

    void setRMax(const double& r) {
        m_r_max = r;
    }

    double getRMax() const {
        return m_r_max;
    }

    double calcDF(size_t j, const std::vector<double>& F) const {
        double sum = 0;

        for (size_t i = 0; i < m_df.size(); ++i) {
            sum += F[i] * m_df[i]->get(j);
        }

        return sum;
    }

    double calcF(const std::vector<double>& F) const {
        double sum = 0;

        for (size_t i = 0; i < m_f.size(); ++i) {
            sum += F[i] * m_f.at(i);
        }

        return sum;
    }

    //Calc d [E_k][j]
    double calcDFE(size_t k, size_t j) const {
        double sum = 0;

        for (size_t i = 0; i < m_df.size(); ++i) {
            sum += m_sd[i]->E[k] * m_df[i]->get(j);
        }

        return sum;
    }

    double calcE(size_t k) {
        double sum = 0;

        for (size_t i = 0; i < m_f.size(); ++i) {
            sum += m_sd[i]->E[k] * m_f[i];
        }

        return sum;
    }

    double calcDFD(size_t k, size_t j) const {
        double sum = 0;

        for (size_t i = 0; i < m_df.size(); ++i) {
            sum += m_sd[i]->D[k] * m_df[i]->get(j);
        }

        return sum;
    }

    double calcD(size_t k) {
        double sum = 0;

        for (size_t i = 0; i < m_f.size(); ++i) {
            sum += m_sd[i]->D[k] * m_f[i];
        }

        return sum;
    }


    double calcDFH(size_t k, size_t j) const {
        double sum = 0;

        for (size_t i = 0; i < m_df.size(); ++i) {
            sum += m_sd[i]->H[k] * m_df[i]->get(j);
        }

        return sum;
    }

    double calcH(size_t k) {
        double sum = 0;

        for (size_t i = 0; i < m_f.size(); ++i) {
            sum += m_sd[i]->H[k] * m_f[i];
        }

        return sum;
    }


    double calcDFB(size_t k, size_t j) const {
        double sum = 0;

        for (size_t i = 0; i < m_df.size(); ++i) {
            sum += m_sd[i]->B[k] * m_df[i]->get(j);
        }

        return sum;
    }

    double calcB(size_t k) {
        double sum = 0;

        for (size_t i = 0; i < m_f.size(); ++i) {
            sum += m_sd[i]->B[k] * m_f[i];
        }

        return sum;
    }


    void resizeE(size_t k) {
        E.resize(k);
    }

    void resizeD(size_t k) {
        D.resize(k);
    }

    void resizeH(size_t k) {
        H.resize(k);
    }

    void resizeB(size_t k) {
        B.resize(k);
    }

    void setE(size_t k, const double& v) {
        E[k] = v;
    }

    void setD(size_t k, const double& v) {
        D[k] = v;
    }

    void setH(size_t k, const double& v) {
        H[k] = v;
    }

    void setB(size_t k, const double& v) {
        B[k] = v;
    }

    double getE(size_t k) {
        return E[k];
    }

    double getD(size_t k) {
        return D[k];
    }

    double getH(size_t k) {
        return H[k];
    }

    double getB(size_t k) {
        return B[k];
    }

    size_t getDim() const {
        return m_p->size();
    }

    ~Point() {
        if ( m_sigma != 0 ) {
            free(m_sigma);
        }

        if ( m_p != 0 ) {
            delete m_p;
        }

        for (size_t i = 0, j = m_df.size(); i < j ; ++i) {
            delete m_df[i];
        }

        for (size_t i = 0; i < m_fpoint.size(); ++i) {
            delete m_fpoint[i];
        }

        m_sd.clear();
    }

    void addMirror(Point* p, Point* fp) {
        if (m_mirror.size() == 0) {
            m_mutexMirror.init();
        }
        m_mirror.push_back(p);
        m_fpoint.push_back(fp);
    }

    size_t sizeMirror() {
        return m_mirror.size();
    }

    void setVector(matrix::Vector* p) {
        m_p = p;
    }

    Point* getMirror(size_t i) {
        return m_mirror[i];
    }

    Point* getFPoint(size_t i) {
        return m_fpoint[i];
    }

    void startMirror(FPDist D) {
        m_mutexMirror.lock();
        for (size_t i = 0; i < m_mirror.size(); ++i) {
            size_t dim = m_p->size();
            Point* p1 = m_mirror[i];
            Point* p2 = m_fpoint[i];

            p1->m_mutexMirror.lock();
            p2->m_mutexMirror.lock();
            for (size_t k = 0; k < dim; ++k) {
                double x1 = p1->getX(k);
                double x2 = p2->getX(k);
                p1->setX(k, x1 + x2);
            }
        }

        m_r_max = 0;
        for (size_t i = 0; i < m_mirror.size(); ++i) {
            Point* p1 = m_mirror[i];
            double d = D(p1, this);
            if (d > m_r_max) {
                m_r_max = d;
            }
        }
    }

    void endMirror() {
        for (size_t i = 0; i < m_mirror.size(); ++i) {
            size_t dim = m_p->size();
            Point* p1 = m_mirror[i];
            Point* p2 = m_fpoint[i];

            for (size_t k = 0; k < dim; ++k) {
                double x1 = p1->getX(k);
                double x2 = p2->getX(k);
                p1->setX(k, x1 - x2);
            }

            p1->m_mutexMirror.unlock();
            p2->m_mutexMirror.unlock();
        }
        m_mutexMirror.unlock();
    }


 private:
    std::vector<const Point*> m_sd;
    std::vector<matrix::Vector*> m_df;
    std::map<size_t, double> m_f;
    std::vector<double> E, D, H, B;

    int m_ID;
    size_t m_D;
    double* m_sigma;
    double m_eps;
    double m_mu;
    double m_r_max;
    int m_prop;
    char m_type;


    matrix::Vector* m_p;

    std::vector<Point*> m_mirror;
    /*
        fpoint é o conjunto de pontos usados
        para corrigir o espelho.
        mirror'[k] = mirror[k] + fpoint[k]
    */
    std::vector<Point*> m_fpoint;
    Mutex m_mutexMirror;

    DISALLOW_COPY_AND_ASSIGN(Point);
};

bool operator<(const Point& a, const Point& b);

struct PointComp {
    bool operator()(const Point * const a, const Point * const b) {
        return *a < *b;
    }
};

typedef std::vector<Point*> Points;
typedef std::vector<size_t> VectorInt;
typedef std::vector<double> VectorDouble;

double accSD(const Points* p);
void free(Points* p);
void print(std::FILE* const f, const Point* const p);
void printAll(std::FILE* const f, const Point* const p);
void print(std::FILE* const  f, const Points* const p);
void printAll(std::FILE* const f, const Points* const p);

double mult(const Point* const a, const Point* const b);
double norm(const Point* const a);
void mult(Point* const  res, const double* const a, const Point* const b);
void mult(Point* const  res, const Point* const a, const double* const b);
void add(Point* const  res, const Point* const a, const Point* const b);
void sub(Point* const  res, const Point* const a, const Point* const b);
void dot(Point* const  res, const Point* const a, const Point* const b);

bool save(const char* const fname, const Points* const p);

void coloumb2D(Points* p, double dx, int NP = 200,  double K = 200, double Rmax = 1e-3, double Rmin = 1e-10);

/**
 * This function load a points set from file "fname".
 * */
bool load(const char* const fname, Points** p);
}}}  // lane_maxwell::modules::point

#endif  // MATH_POINT_H_
