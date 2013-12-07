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
@file rpim.cc


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include <map>
#include "rpim/rpim.h"
#include "maxwell/maxwell.h"

#ifndef M_PI
const double M_PI = acos(-1);
#endif

namespace lane_maxwell { namespace modules { namespace rpim {
double phi(int M, const point::Point* p, const point::Point* c,
            double beta, double raio) {
    double r;

    M = p->size();
    r = 0;

    for (int i = 0; i < M; ++i) {
        double d = p->getX(i) - c->getX(i);
        r += d * d;
    }

    r = sqrt(r) / raio;

    return exp(-beta * r * r);
}

double dphi(int k, const point::Point* p, const point::Point* c,
            double beta, double raio) {
    return -2.0 * beta * (p->getX(k) - c->getX(k)) *
            phi(k, p, c, beta, raio) / (raio * raio);
}


void calcF(bool* const st, size_t I, size_t ps, point::Points* P,
           CallDF* const cdf, double beta, RadialFunction R,
           point::FPDist dist) {
    point::Point* p = (*P)[I];

    size_t D = p->size();

    size_t N = ps, M = D+1;

    if (p->sizeMirror() > 0) {
        p->startMirror(dist);
    }

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            cdf->Ro->set(i, j, R(D, p->getSD(i), p->getSD(j), beta, p->getRMax()));
        }
    }

    for (size_t i = 0; i < N; ++i) {
        cdf->Po->set(i, 0, 1.0);

        for (size_t j = 1; j < M; ++j) {
            cdf->Po->set(i, j, p->getSD(i)->getX(j-1));
        }
    }

    matrix::trans(cdf->tPo, cdf->Po);

    cdf->ciRo->newCall(cdf->Ro);
    matrix::inv(cdf->ciRo, cdf->ciRo->A);
    cdf->iRo = cdf->ciRo->Ai;

    *st = cdf->ciRo->status;

    if (!(*st)) {

        if (p->sizeMirror() > 0) {
            p->endMirror();
        }

        return;
    }

    matrix::mult(cdf->aux, cdf->tPo, cdf->iRo);
    matrix::mult(cdf->auxPo, cdf->aux, cdf->Po);

    cdf->ciauxPo->newCall(cdf->auxPo);
    matrix::inv(cdf->ciauxPo, cdf->ciauxPo->A);
    cdf->iauxPo = cdf->ciauxPo->Ai;

    *st = cdf->ciauxPo->status;

    if (!(*st)) {

        if (p->sizeMirror() > 0) {
            p->endMirror();
        }

        return;
    }

    matrix::mult(cdf->Sb, cdf->iauxPo, cdf->aux);
    matrix::mult(cdf->iRoPo, cdf->iRo, cdf->Po);
    matrix::mult(cdf->iRoPoSb, cdf->iRoPo, cdf->Sb);
    matrix::sub(cdf->Sa, cdf->iRo, cdf->iRoPoSb);


    for (size_t k = 0; k < N; ++k) {
        p->setF(k, 0);

        for (size_t i = 0; i < N; ++i) {
            p->setF(k, p->getF(k) + R(0, p, p->getSD(i),
                 beta, p->getRMax()) * cdf->Sa->get(i, k));
        }

        for (size_t i = 0; i < M; ++i) {
            if (i == 0) {
                p->setF(k, p->getF(k) + cdf->Sb->get(i, k));
            } else {
                p->setF(k, p->getF(k) + cdf->Sb->get(i, k) * p->getX(i - 1) );
            }
        }
    }


    if (p->sizeMirror() > 0) {
        p->endMirror();
    }
}


void calcDF(bool* const st, size_t I, size_t ps, point::Points* P,
            CallDF* const cdf, double beta, RadialFunction R,
            RadialFunction dR, point::FPDist dist) {
    point::Point* p = (*P)[I];

    size_t D = p->size();

    size_t N = ps, M = D+1;


    if (p->sizeMirror() > 0) {
        p->startMirror(dist);
    }

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            cdf->Ro->set(i, j, R(D, p->getSD(i), p->getSD(j), beta, p->getRMax()));
        }
    }

    for (size_t i = 0; i < N; ++i) {
        cdf->Po->set(i, 0, 1.0);

        for (size_t j = 1; j < M; ++j) {
            cdf->Po->set(i, j, p->getSD(i)->getX(j-1));
        }
    }

    matrix::trans(cdf->tPo, cdf->Po);

    cdf->ciRo->newCall(cdf->Ro);
    matrix::inv(cdf->ciRo, cdf->ciRo->A);
    cdf->iRo = cdf->ciRo->Ai;

    *st = cdf->ciRo->status;

    if (!(*st)) {

        if (p->sizeMirror() > 0) {
            p->endMirror();
        }

        return;
    }

    matrix::mult(cdf->aux, cdf->tPo, cdf->iRo);
    matrix::mult(cdf->auxPo, cdf->aux, cdf->Po);

    cdf->ciauxPo->newCall(cdf->auxPo);
    matrix::inv(cdf->ciauxPo, cdf->ciauxPo->A);
    cdf->iauxPo = cdf->ciauxPo->Ai;

    *st = cdf->ciauxPo->status;

    if (!(*st)) {

        if (p->sizeMirror() > 0) {
            p->endMirror();
        }

        return;
    }

    matrix::mult(cdf->Sb, cdf->iauxPo, cdf->aux);
    matrix::mult(cdf->iRoPo, cdf->iRo, cdf->Po);
    matrix::mult(cdf->iRoPoSb, cdf->iRoPo, cdf->Sb);
    matrix::sub(cdf->Sa, cdf->iRo, cdf->iRoPoSb);


    for (size_t k = 0; k < N; ++k) {
        for (size_t j = 0; j < D; ++j) {
            p->setDF(k, j, 0);

            for (size_t i = 0; i < N; ++i) {
                p->setDF(k, j, p->getDF(k, j) + dR(j, p, p->getSD(i),
                     beta, p->getRMax())*cdf->Sa->get(i, k));
            }

            p->setDF(k, j, p->getDF(k, j) + cdf->Sb->get(j + 1, k));
        }
    }


    if (p->sizeMirror() > 0) {
        p->endMirror();
    }
}



void calcDF(int kD, bool* const st, size_t I, size_t ps, point::Points* P,
            CallDF* const cdf, double beta, RadialFunction R, RadialFunction dR,
            point::FPDist dist) {

    point::Point* p = (*P)[I];

    size_t D = p->size();

    size_t N = ps;
    size_t M = D+1;


    if (p->sizeMirror() > 0) {
        p->startMirror(dist);
    }

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            cdf->Ro->set(i, j, R(D, p->getSD(i), p->getSD(j), beta, p->getRMax()));
        }
    }

    for (size_t i = 0; i < N; ++i) {
        cdf->Po->set(i, 0, 1.0);

        for (size_t j = 1; j < M; ++j) {
            cdf->Po->set(i, j, p->getSD(i)->getX(j-1) - P->at(I)->getX(j-1));
        }
    }

    matrix::trans(cdf->tPo, cdf->Po);

    cdf->ciRo->newCall(cdf->Ro);
    matrix::inv(cdf->ciRo, cdf->ciRo->A);
    cdf->iRo = cdf->ciRo->Ai;

    *st = cdf->ciRo->status;

    if (!(*st)) {

        if (p->sizeMirror() > 0) {
            p->endMirror();
        }

        return;
    }

    matrix::mult(cdf->aux, cdf->tPo, cdf->iRo);
    matrix::mult(cdf->auxPo, cdf->aux, cdf->Po);

    cdf->ciauxPo->newCall(cdf->auxPo);
    matrix::inv(cdf->ciauxPo, cdf->ciauxPo->A);
    cdf->iauxPo = cdf->ciauxPo->Ai;

    *st = cdf->ciauxPo->status;

    if (!(*st)) {

        if (p->sizeMirror() > 0) {
            p->endMirror();
        }

        return;
    }

    matrix::mult(cdf->Sb, cdf->iauxPo, cdf->aux);
    matrix::mult(cdf->iRoPo, cdf->iRo, cdf->Po);
    matrix::mult(cdf->iRoPoSb, cdf->iRoPo, cdf->Sb);
    matrix::sub(cdf->Sa, cdf->iRo, cdf->iRoPoSb);


    for (size_t k = 0; k < N; ++k) {
        p->setDF(k, kD, 0);

        for (size_t i = 0; i < N; ++i) {
            p->setDF(k, kD, p->getDF(k, kD) + dR(kD, p, p->getSD(i),
                 beta, p->getRMax())*cdf->Sa->get(i, k));
        }

        p->setDF(k, kD, p->getDF(k, kD) + cdf->Sb->get(kD + 1, k));
    }


    if (p->sizeMirror() > 0) {
        p->endMirror();
    }
}

double calFunction(const point::Point* const p, const double& cnst, const point::Point* const c) {
    double acc = 0;
    for (size_t i = 0; i < p->getDim(); ++i) {
            acc += std::sin(cnst * (p->getX(i) - c->getX(i)));
    }
    return acc;
}

double dcalFunction(int k, const point::Point* const p, const double& cnst, const point::Point* const c) {
    return cnst * std::cos(cnst *  (p->getX(k) - c->getX(k)));
}



double bestBetaF(bool* const st, size_t I, size_t ps, point::Points* P,
    CallDF* const cdf, RadialFunction R, double fmax,
    double betamin, double betamax, int NP, double errmax, std::vector<double>* F,
    double* err) {
    double C = 299792458.0e0;
    double cnst = 2 * M_PI * fmax / C;
    double dbeta = (betamax - betamin) / NP;
    double DF;
    bool fib = false;

    if (P->at(I)->sizeMirror() > 0) {
        P->at(I)->startMirror(point::distE);
    }

    for (size_t i = 0; i < ps; ++i) {
        (*F)[i] = calFunction(P->at(I)->getSD(i), cnst, P->at(I));
    }
    DF = calFunction(P->at(I), cnst, P->at(I));


    if (P->at(I)->sizeMirror() > 0) {
        P->at(I)->endMirror();
    }

    bool ok = false;
    double bestBeta = betamin;
    double errmin = 100;
    double a = 0, b = 0, fa = 0, fb = 0 , fx = 0, fxo = 0;

    for (double beta = betamin; beta < betamax; beta += dbeta) {
        double xa = beta;
        double xb = beta + 0.5 * dbeta;

        calcF(st, I, ps, P, cdf, xa * xa, R);

        double dfa = P->at(I)->calcF(*F) - DF;
        double err = std::fabs(dfa / (fabs(DF) + 1));

        if (!(*st)) err = 100;

        if (matrix::cmpD(err,errmin) < 0 && (*st)) {
            bestBeta = xa;
            fib =true;
            errmin = err;
        }


        calcF(st, I, ps, P, cdf, xb * xb, R);

        double dfb = P->at(I)->calcF(*F) - DF;
        err = std::fabs(dfb / (fabs(DF) + 1));

        if (!(*st)) err = 100;

        if (matrix::cmpD(err,errmin) < 0 && (*st)) {
            bestBeta = xb;
            errmin = err;
            fib = true;
        }

        if (dfa * dfb < 0 && (*st) ) {
            ok = true;
            fib = true;
            fa = dfa;
            fb = dfb;
            a = xa;
            b = xb;
        }

    }

    if (ok && errmin < 10) {
        int ni = 0;
        double err;
        double x = a;

        fx = fa;

        while (ni < NP && errmin > errmax) {
            x = (a * fb - b * fa)/(fb - fa);

            calcF(st, I, ps, P, cdf, x * x, R);

            fxo = fx;
            fx = P->at(I)->calcF(*F) - DF;

            if (matrix::cmpD(fx * fa, 0) >= 0 ) {
                a = x;

                fa = fx;
                if (matrix::cmpD(fx * fxo, 0) >= 0) {
                    fb = fb * 0.5;
                }
            } else {
                b = x;
                fb = fx;
                if (matrix::cmpD(fx * fxo, 0) >= 0) {
                    fa = fa * 0.5;
                }
            }

            err= std::fabs(fx / (fabs(DF) + 1));
            if (err < errmin && (*st)) {
                errmin = err;
                bestBeta = x;
                fib = true;
            }

            ++ni;
        }
    }

    if (!fib) *st = false;

    if (!fib) *err = 100;
    else *err = 100*errmin;

    return bestBeta;

}


double bestBeta(int k, bool* const st, size_t I, size_t ps, point::Points* P,
    CallDF* const cdf, RadialFunction R, RadialFunction dR, double fmax,
    double betamin, double betamax, int NP, double errmax, std::vector<double>* F,
    double* err) {
    double C = 299792458.0e0;
    double cnst = 2 * M_PI * fmax / C;
    double dbeta = (betamax - betamin) / NP;
    double DF;
    bool fib = false;


    if (P->at(I)->sizeMirror() > 0) {
        P->at(I)->startMirror(point::distE);
    }

    for (size_t i = 0; i < ps; ++i) {
        (*F)[i] = calFunction(P->at(I)->getSD(i), cnst, P->at(I));
    }
    DF = dcalFunction(k, P->at(I), cnst, P->at(I));


    if (P->at(I)->sizeMirror() > 0) {
        P->at(I)->endMirror();
    }

    bool ok = false;
    double bestBeta = betamin;
    double errmin = 100;
    double a = 0, b = 0, fa = 0, fb = 0 , fx = 0, fxo = 0;

    for (double beta = betamin; beta < betamax; beta += dbeta) {
        double xa = beta;
        double xb = beta + 0.5 * dbeta;

        calcDF(k, st, I, ps, P, cdf, xa * xa, R, dR);

        double dfa = P->at(I)->calcDF(k, *F) - DF;
        double err = std::fabs(dfa / (fabs(DF) + 1));

        if (!(*st)) err = 100;

        if (matrix::cmpD(err,errmin) < 0 && (*st)) {
            bestBeta = xa;
            fib =true;
            errmin = err;
        }


        calcDF(k, st, I, ps, P, cdf, xb * xb, R, dR);

        double dfb = P->at(I)->calcDF(k, *F) - DF;
        err = std::fabs(dfb / (fabs(DF) + 1));

        if (!(*st)) err = 100;

        if (matrix::cmpD(err,errmin) < 0 && (*st)) {
            bestBeta = xb;
            errmin = err;
            fib = true;
        }

        if (dfa * dfb < 0 && (*st) ) {
            ok = true;
            fib = true;
            fa = dfa;
            fb = dfb;
            a = xa;
            b = xb;
        }

    }

    if (ok && errmin < 10) {
        int ni = 0;
        double err;
        double x = a;

        fx = fa;

        while (ni < NP && errmin > errmax) {
            x = (a * fb - b * fa)/(fb - fa);

            calcDF(k, st, I, ps, P, cdf, x * x, R, dR);

            fxo = fx;
            fx = P->at(I)->calcDF(k, *F) - DF;

            if (matrix::cmpD(fx * fa, 0) >= 0 ) {
                a = x;

                fa = fx;
                if (matrix::cmpD(fx * fxo, 0) >= 0) {
                    fb = fb * 0.5;
                }
            } else {
                b = x;
                fb = fx;
                if (matrix::cmpD(fx * fxo, 0) >= 0) {
                    fa = fa * 0.5;
                }
            }

            err= std::fabs(fx / (fabs(DF) + 1));
            if (err < errmin && (*st)) {
                errmin = err;
                bestBeta = x;
                fib = true;
            }

            ++ni;
        }
    }

    if (!fib) *st = false;

    if (!fib) *err = 100;
    else *err = 100*errmin;

    return bestBeta;
}


bool save(const char* const fname, const point::Points* const p) {
    std::FILE* F = fopen(fname, "wb+");

    if (!F) {
        return false;
    }

    size_t N = p->size();
    int I = 0;
    point::Points::const_iterator it = p->begin(), end = p->end();

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

    it = p->begin();
    while (it != end) {
        size_t ps = (*it)->sizeSD();
        int prop = (*it)->getProp();

        fwrite(&ps, sizeof(ps), 1, F);
        for (size_t i = 0; i < ps; ++i) {
            I = (*it)->getSD(i)->getID();
            fwrite(&I, sizeof(I), 1, F);
            for (size_t k = 0; k < dim; ++k) {
                double coef = (*it)->getDF(i, k);
                fwrite(&coef, sizeof(coef), 1, F);
            }

            if (prop & maxwell::P_MED) {
                double coef = (*it)->getF(i);
                fwrite(&coef, sizeof(coef), 1, F);
            }
        }
        point::Point* pv = *it;
        ps = pv->sizeMirror();

        fwrite(&ps, sizeof(ps), 1, F);

        for (size_t i = 0; i < ps; ++i) {
            I = pv->getMirror(i)->getID();
            fwrite(&I, sizeof(I), 1, F);
            for (size_t k = 0; k < dim; ++k) {
                double x = pv->getFPoint(i)->getX(k);
                fwrite(&x, sizeof(x), 1, F);
            }
        }
        ++it;
    }

    fclose(F);
    return true;
}

bool load(const char* const fname, point::Points** p) {
    std::FILE* F = fopen(fname, "rb");

    if (!F) {
        return false;
    }

    std::map<int, point::Point*> mp;
    size_t N, dim;

    rewind(F);

    fread(&N, sizeof(N), 1, F);
    fread(&dim, sizeof(dim), 1, F);

    size_t I = 0;

    *p = new point::Points;

    while (I < N) {
        char type;
        double X;
        double eps;
        double mu;

        int ID;
        int prop;

        point::Point* np = new point::Point;
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

        mp[ID] = np;
        (*p)->push_back(np);

        ++I;
    }

    I = 0;

    while (I < N) {
        size_t ps;
        int ID;
        point::Point* pv = (*p)->at(I);
        int prop = pv->getProp();

        fread(&ps, sizeof(ps), 1, F);
        for (size_t i = 0; i < ps; ++i) {
            fread(&ID, sizeof(ID), 1, F);
            point::Point* pt = mp[ID];

            pv->addSD(pt);
            for(size_t k = 0; k < dim; ++k) {
                double coef;
                fread(&coef, sizeof(coef), 1, F);
                pv->setDF(i, k, coef);
            }

            if (prop & maxwell::P_MED) {
                double coef;
                fread(&coef, sizeof(coef), 1, F);
                pv->setF(i, coef);
            }
        }

        /*Le pontos espelho*/
        fread(&ps, sizeof(ps), 1, F);
        for (size_t i = 0; i < ps; ++i) {
            fread(&ID, sizeof(ID), 1, F);

            point::Point* pt = mp[ID];
            matrix::Vector* vp = new matrix::Vector(dim);

            for(size_t k = 0; k < dim; ++k) {
                double x;
                fread(&x, sizeof(x), 1, F);
                vp->set(k, x);
            }
            point::Point* np = new point::Point();
            np->setVector(vp);
            pv->addMirror(pt, np);
        }

        ++I;
    }

    fclose(F);
    return true;
}

}}}  // lane_maxwell::modules::rpim

