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

#ifndef RUNSIMULATION_H
#define RUNSIMULATION_H

#include "genDom.h"
#include "maxwell/maxwell.h"
#include "lane/lane.h"

typedef double (*PFSource)(const double& t, const double& dt,
                           const point::Point* pts);

enum SourceMode {
    HARD = 0,
    SOFT = 1,
    TYPE3 = 2
};

enum SourceType {
    SIN = 0,
    COS = 1,
    GAUSS = 2,
    SIGMOID = 3,
    MONOCICLOGAUSS = 4,
    DESCARGAATM = 5,
    USER = 6
};

class RunSimulation;

class MaxwellT: public QThread
{
    Q_OBJECT
public:
    void run();
    bool updateE;
    bool isTMZ;
    size_t ini, fim;
    RunSimulation* data;
};

class RunSimulation : public QThread
{
    Q_OBJECT
public:
    RunSimulation()
        :pts(0), status(false), stopNow(false){
        resImgs = new QStringList;
    }
    ~RunSimulation() {
        delete resImgs;
    }

    bool load();
    void run();

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

    void free() {
        if (pts != 0) {
            lane_maxwell::modules::point::free(pts);
        }

        for (size_t i = 0; i < threads.size(); ++i) {
            delete threads[i];
        }
    }

    QString pwd;
    QString buffer;
    char data[1000];
    point::Points* pts;
    bool isTMZ;
    double Ts;
    double dt;
    double cte;
    double resX;
    double resY;
    bool videoGen;
    double w;
    double freq;
    int nf;
    bool isLog;
    int lecont;
    size_t NTHREADS;
    PFSource source;
    int sourceType;
    double duracaoPulso;
    double pico;
    int ttype;
    std::vector<double> fonte;
    std::vector<MaxwellT*> threads;
    bool status;
    bool stopNow;

    QStringList* resImgs;
    QMutex mutex;

signals:
    void setLog(QString);
    void setEnabled();
};

#endif // RUNSIMULATION_H
