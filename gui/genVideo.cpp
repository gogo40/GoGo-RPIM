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

#include "genVideo.h"

void GenVideoT::init() {
    mutex.lock();

    ispause = false;
    isstop = false;
    iswait = false;

    clock = 0;
    start();
    mutex.unlock();
}

void GenVideoT::run() {

    int t = static_cast<int>(1000*dt);
    int i = 0;
    while(true) {
        while (isPaused()) sleep(1);
        if (isStoped()) break;

        setTime(i * (t/1000.0));
        emit nextImage();

        waitSignal();
        while (true) {
            QThread::msleep(t);
            if (isContinue()) break;
            if (isStoped()) goto END;
        }
        ++i;
        if (i == N) break;
    }
    END:
    mutex.lock();
    isstop = false;
    mutex.unlock();
    emit finishVideo();
}
