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

#ifndef GENDOM_H
#define GENDOM_H

#include <QDebug>
#include <QThread>
#include <QProcess>
#include <QObject>
#include <QMutex>
#include <vector>
#include <QString>
#include <QFileDialog>
#include <QMainWindow>
#include <QTextStream>
#include <QTime>

#include "maxwell/maxwellgendom.h"
#include "lane/lane.h"

class RunDomainGeneration;

using namespace lane_maxwell;
using namespace lane_maxwell::modules::maxwell;
using namespace lane_maxwell::modules;

typedef std::map<int, rpim::CallDF*> MapCDF;
typedef MapCDF::iterator CDFIt;
typedef std::pair<int, rpim::CallDF*> PairCDF;

class CalcDFT : public QThread
{
     Q_OBJECT

public:


    void run();


    int id;
    size_t ini, fim;
    bool status;

    RunDomainGeneration* data;

signals:
    void setLog(QString);

};


class RunDomainGeneration : public QThread
{
    Q_OBJECT

public:

    RunDomainGeneration(QString* buffer, FILE* st, FILE* serr)
        : mbuffer(buffer), s(buffer), sout(st), serr(serr), p(0),
          stopNow(false) {
        mbuffer->clear();
    }

    void run();

    bool load();

    void free() {

        if (p != 0){
            lane_maxwell::modules::point::free(p);
        }

        for (size_t i = 0; i < threads.size(); ++i) {
            delete threads[i];
            delete vcdf[i];
            delete vbetax[i];
            delete vFbetax[i];
        }
    }

    void online() {
        mutex.lock();
        stopNow = false;
        mutex.unlock();
    }

    void stop() {
        mutex.lock();
        stopNow = true;
        mutex.unlock();
    }

    bool isStoped() {
        bool s;
        mutex.lock();
        s = stopNow;
        mutex.unlock();
        return s;
    }


    QString* mbuffer;
    QTextStream s;
    QTextStream sout;
    QTextStream serr;


    char data[1000];
    int domainSize;
    int NTHREADS;
    int metric;
    double beta;
    std::vector<double> betax;
    std::vector<CalcDFT*> threads;
    std::vector<MapCDF*> vcdf;
    std::vector<std::vector<double>* > vbetax;
    std::vector<std::vector<double>* > vFbetax;
    bool useR;
    double ratio;
    int radialBase;
    bool isAuto;
    double betamin;
    double betamax;
    double errbeta;
    int NP;
    double fmax;
    lane_maxwell::modules::rpim::RadialFunction R;
    lane_maxwell::modules::rpim::RadialFunction dR;
    lane_maxwell::modules::domains::FPDist D;
    lane_maxwell::modules::point::Points* p;
    bool stopNow;
    QMutex mutex;
public slots:
    void setTlog(QString);

signals:
    void setLog(QString);
    void setEnabled();
};


#endif // GENDOM_H
