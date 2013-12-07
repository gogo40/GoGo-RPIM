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

#ifndef GENVIDEO_H
#define GENVIDEO_H

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

class GenVideoT : public QThread {
    Q_OBJECT

public:
    GenVideoT()
        : iswait(false), isstop(false), ispause(false) {}
    void run();
    void init();

    bool isPaused() {
        bool p;
        mutex.lock();
        p = ispause;
        mutex.unlock();
        return p;
    }

    bool isStoped() {
        bool s;
        mutex.lock();
        s = isstop;
        mutex.unlock();
        return s;
    }

    double getTime() {
        double t;
        mutex.lock();
        t = dt;
        mutex.unlock();
        return t;
    }

    void setTime(double t) {
        mutex.lock();
        clock = t;
        mutex.unlock();
    }

    void pause() {
        mutex.lock();
        ispause=true;
        mutex.unlock();
    }

    void replay() {
        mutex.lock();
        ispause = false;
        iswait = false;
        isstop = false;
        mutex.unlock();
    }

    void stop() {
        mutex.lock();
        isstop = true;
        mutex.unlock();
    }


    void waitSignal() {
        mutex.lock();
        iswait = true;
        mutex.unlock();
    }

    bool isContinue() {
        bool w;
        mutex.lock();
        w = iswait;
        mutex.unlock();
        return !w;
    }

    double dt;
    bool iswait;
    bool isstop;
    bool ispause;
    int N;
    QMutex mutex;
    double clock;

public slots:
    void nextVideoImage() {
        mutex.lock();
        iswait = false;
        ispause = false;
        mutex.unlock();
    }

signals:
    void nextImage();
    void finishVideo();
};

#endif // GENVIDEO_H
