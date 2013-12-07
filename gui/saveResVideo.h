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

#ifndef SAVERESVIDEO_H
#define SAVERESVIDEO_H

#include <QGraphicsPixmapItem>
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

class SaveResVideoT : public QThread {
    Q_OBJECT

public:

    SaveResVideoT()
    : ResImgs(0), fileName(0) {}
    ~SaveResVideoT(){}

    void run();
    void init();

    bool imgLoaded(){
        bool c;
        m.lock();
        c = loadedi;
        m.unlock();
        return c;
    }

    void setLoaded(bool c) {
        m.lock();
        loadedi = c;
        m.unlock();
    }

    QStringList* ResImgs;
    QString* fileName;
    QPixmap* img;
    QMutex m;
    bool loadedi;

signals:
    void sendLog(QString);
    void loadImg(QString);

};

#endif // SAVERESVIDEO_H
