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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include "runSimulation.h"
#include "genVideo.h"
#include "savevideoform.h"

namespace Ui {
class MainWindow;
}



class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void init();

public slots:
    void enableGen();
    void enableSim();
    void setGenLog(QString);
    void setSimLog(QString);
    void setPWD();
    void runSimulation();
    void setAuto();
    void setUseRatio();
    void runDomainGeneration();
    void setInputNameSim();
    void setInputNameDom();
    void printHelpDom();
    void printHelpSim();
    void setSourceData();
    void drawResImage();
    void nextResImage();
    void prevResImage();
    void scaleResImage(double);
    void scaleResImageX(double);
    void scaleResImageY(double);
    void stopSimulation();
    void stopGeneration();

    void playResVideo();
    void pauseResVideo();
    void stopResVideo();
    void renderResVideo();
    void saveResVideo();
    void loadResVideo();
    void resetSaveVideo();

    QPixmap getNextImg();

private:
    Ui::MainWindow *ui;
    RunSimulation* simT;
    RunDomainGeneration genT;
    GenVideoT genVideoT;
    saveVideoForm* uiSaveVideo;
    QString pwd;

    int currImg;
    double resScale;
    double resScaleX;
    double resScaleY;
    char buffer[1000];

    bool mloadResVideo;
    QString fileName;
    int NFrames;
    int ip;


    QFile* fileVideo;
    QDataStream* videoStream;

    QMutex mRes;

signals:
    void sendSimLog(QString&);
    void nextResVideo();
};



#endif // MAINWINDOW_H
