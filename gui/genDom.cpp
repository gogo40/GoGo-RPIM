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

#include "genDom.h"

#include <algorithm>
#include <utility>
#include <map>


void CalcDFT::run() {
    QTime t;
    size_t T = data->domainSize;
    int prop;
    bool useR = data->useR;
    point::Points* P = data->p;
    MapCDF* cdf = data->vcdf[id];
    double beta = data->beta;
    bool isauto = data->isAuto;
    double betamin = data->betamin;
    double betamax = data->betamax;
    double errbeta = data->errbeta;
    double err;
    int NP = data->NP;
    double fmax = data->fmax;
    std::vector<double>* F = data->vFbetax[id];
    std::vector<double>* betax = data->vbetax[id];
    rpim::RadialFunction R = data->R;
    rpim::RadialFunction dR = data->dR;

    t.start();
    for (size_t i = ini; i < fim; ++i) {

        T = P->at(i)->sizeSD();

        if (useR) {
            prop = P->at(i)->getProp();
            if (prop & P_PEC) continue;
        }

        if (isauto) {
            for (size_t j = 0; j < P->at(i)->getDim(); ++j) {
                (*betax)[j] = rpim::bestBeta(j, &(status),
                                             i, T, P,
                                             (*cdf)[T], R, dR, fmax,
                                             betamin, betamax, NP,
                                             errbeta, F, &err);

                rpim::calcDF(j, &(status), i, T, P, (*cdf)[T],
                (*betax)[j] * (*betax)[j], R, dR);

                if (!status) {
                    data->mutex.lock();
                    data->s << "[" << id <<"]";
                    data->s << "Error to calc derivate "
                            << j << " p "
                            << i << "\n";
                    emit setLog(data->s.readAll());
                    data->mutex.unlock();
                    //break;
                }
            }

            if (t.elapsed()/1000.0 > 2.0) {
                data->mutex.lock();
                data->s << "[" <<id <<"]\n";
                for (size_t j = 0; j < P->at(i)->getDim(); ++j) {
                    data->s << "beta " << j+1 << " = "
                            << (*betax)[j] << "err = "
                            << err << "\n";
                }
                data->mutex.unlock();
            }
        } else {
            rpim::calcDF(&(status), i, T, P, (*cdf)[T], beta,  R, dR);
            if (!status) {
                data->mutex.lock();
                data->s << "[" <<id
                        <<"] Erro  to calc derivate " << i <<"\n";
                data->mutex.unlock();
            }
        }

        if ( t.elapsed()/1000.0 > 2.0) {
            t.restart();
            data->mutex.lock();
            data->s << "[" <<
                       id << "] " <<
                       i+1 <<"/" << fim <<"\n";
            emit setLog(data->s.readAll());
            data->mutex.unlock();

            if (data->isStoped()) {
                emit setLog("Stopped!");
                break;
            }
        }
    }
}

void RunDomainGeneration::setTlog(QString s) {
    emit setLog(s);
}

bool RunDomainGeneration::load() {
    online();
    double C = 299792458.0e0;
    double f = this->fmax;
    double lmb = C / f;

    s << "data = " << this->data << "\n";
    s << "N Threads = " << this->NTHREADS << "\n";
    s << "metric = ";
    switch (this->metric) {
    case 0:
        s << "Euclidian\n";
    break;

    case 1:
        s << "Manhattan\n";
    break;

    case 2:
        s << "Max\n";
    break;
    }

    s << "\nRBF = ";
    switch (this->radialBase) {
    case 0:
         s << "Gaussian\n";
    break;
    }
    s<<"\n";

    if (this->isAuto) {
        s << "Mode = Automatic Beta\n";
        s << "Beta = [" << this->betamin
          << ", " << this->betamax << "]\n";
        s << "ErrBeta = " << this->errbeta << "\n";
        s << "NP = " << this->NP << "\n";
    } else {
        s << "Mode = Fixed Beta\n";
        s << "Beta = "<<this->beta <<"\n";
    }
    s << "Maximum frequency = " << this->fmax << "\n";
    s << "Lmb = " << lmb << "\n";
    if (this->useR) {
        s << "Use ratio\n";
        s << "r = "<<  this->ratio * lmb / 20
          << " [" << this->ratio << "]\n";
        this->ratio = this->ratio * lmb / 20;
    } else {
        s << "Domain size = " << this->domainSize << "\n";
    }


    emit setLog(s.readAll());

    threads.resize(NTHREADS);
    vcdf.resize(NTHREADS);
    vbetax.resize(NTHREADS);
    vFbetax.resize(NTHREADS);


    for (size_t i = 0; i < threads.size(); ++i) {
        threads[i] = new CalcDFT;
        vcdf[i] = new MapCDF;
        threads[i]->data = this;
        threads[i]->id = i;
        QObject::connect(
                    threads[i], SIGNAL(setLog(QString)),
                    this, SLOT(setTlog(QString)));
    }


    return true;

}

void RunDomainGeneration::run() {

    s << "Loading file ...\n";


    emit setLog(s.readAll());

    if (!lane_maxwell::modules::point::load(this->data, &p)) {
        p = 0;
        s << "Error: Failed to open " << this->data <<"\n";
        emit setLog(s.readAll());
        emit setEnabled();
        return;
    }

    s << "File Loaded!\n";

    QTime t, tt;

    s << "Building support domain E\n";
    emit setLog(s.readAll());

    t.start();
    tt.start();

    std::sort(p->begin(), p->end(), point::PointComp());
    domains::buildSupportDomain(stdout, 'E', 'H', domainSize, p,
                                0.01, D, 1, useR, ratio);
    s << "Time elapsed = " << t.elapsed()/1000.0 << "s\n";

    s << "Building support domain H\n";
    emit setLog(s.readAll());

    t.restart();
    domains::buildSupportDomain(stdout, 'H', 'E', domainSize, p,
                                0.01, D, 1, useR, ratio);

    s << "Time elapsed = " << t.elapsed()/1000.0 << "s\n";
    s << "Calculating coeficients...\n";
    emit setLog(s.readAll());

    size_t dim = p->at(0)->getDim();
    size_t N = p->size();
    size_t mT = 0;

    for (size_t i = 0; i < p->size(); ++i) {
        size_t T = p->at(i)->sizeSD();
        for (int k = 0; k < NTHREADS; ++k) {
            if ( vcdf[k]->find(T) == vcdf[k]->end()) {
                vcdf[k]->insert(PairCDF(T, new rpim::CallDF(T, dim)));
            }
        }
        if (T > mT) {
            mT = T;
        }
    }

    for (int i = 0; i < NTHREADS; ++i) {
        vbetax[i] = (isAuto)?
                    new std::vector<double>(p->at(0)->getDim()):0;
        vFbetax[i] = (isAuto)?
                    new std::vector<double>(mT):0;
    }


    int ft = N / threads.size();
    t.restart();
    for (size_t i = 0; i < threads.size(); ++i) {
        threads[i]->ini = i * ft;
        threads[i]->fim = min((i+1)*ft, N);
        threads[i]->start();
    }

    for (size_t i = 0; i < threads.size(); ++i) {
        threads[i]->wait();
    }

    s << "Time elapsed = " << t.elapsed()/1000.0 << "s\n";
    emit setLog(s.readAll());


    bool status = true;

    if (status) {
        if (!lane_maxwell::modules::rpim::save(this->data, p)) {
            p=0;
            s << "Failed to save file "<<this->data<<"\n";
        }
    } else {
        s << "Failed to calc RPIM coeficients!\n";
    }
    s << "Finished!\n";

    s << "Time total elapsed = " << tt.elapsed()/1000.0 << "s\n";

    emit setLog(s.readAll());
    emit setEnabled();
}

