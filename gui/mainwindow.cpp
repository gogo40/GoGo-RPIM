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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "maxwell/tmz.h"
#include <cstdio>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    simT(new RunSimulation()),
    genT(new QString, stdout, stderr),
    genVideoT(),
    uiSaveVideo(new saveVideoForm),
    pwd(""),
    fileName("")
{

    ui->setupUi(this);
    resScale = 1;
    resScaleX = 1;
    resScaleY = 1;
    currImg = 0;
    ui->nextResImage->setEnabled(false);
    ui->prevResImage->setEnabled(false);

    QObject::connect(ui->actionOpen_Project, SIGNAL(triggered()),
                     this, SLOT(setPWD()));
    ui->actionOpen_Project->setShortcut(tr("Ctrl+A"));

    QObject::connect(simT, SIGNAL(setLog(QString)),
                     this, SLOT(setSimLog(QString)));
    QObject::connect(simT, SIGNAL(setEnabled()),
                     this, SLOT(enableSim()));

    QObject::connect(&genT, SIGNAL(setLog(QString)),
                     this, SLOT(setGenLog(QString)));
    QObject::connect(&genT, SIGNAL(setEnabled()),
                     this, SLOT(enableGen()));

    QObject::connect(ui->getSourceData, SIGNAL(clicked()),
                     this, SLOT(setSourceData()));

    QObject::connect(ui->getInputNameDom, SIGNAL(clicked()),
                     this, SLOT(setInputNameDom()));

    QObject::connect(ui->helpSim, SIGNAL(clicked()),
                     this, SLOT(printHelpSim()));

    QObject::connect(ui->helpDom, SIGNAL(clicked()),
                     this, SLOT(printHelpDom()));

    QObject::connect(ui->prevResImage, SIGNAL(clicked()),
                     this, SLOT(prevResImage()));

    ui->prevResImage->setShortcut(tr("Left"));

    QObject::connect(ui->nextResImage, SIGNAL(clicked()),
                     this, SLOT(nextResImage()));

    ui->nextResImage->setShortcut(tr("Right"));

    QObject::connect(ui->stopGen, SIGNAL(clicked()),
                     this, SLOT(stopGeneration()));


    QObject::connect(ui->stopSim, SIGNAL(clicked()),
                     this, SLOT(stopSimulation()));
    ui->stopSim->setEnabled(false);
    ui->stopGen->setEnabled(false);

    QObject::connect(ui->ResScale, SIGNAL(valueChanged(double)),
                     this, SLOT(scaleResImage(double)));

    QObject::connect(ui->ResScaleX, SIGNAL(valueChanged(double)),
                     this, SLOT(scaleResImageX(double)));

    QObject::connect(ui->ResScaleY, SIGNAL(valueChanged(double)),
                     this, SLOT(scaleResImageY(double)));

    ui->ResScale->setValue(1);
    ui->ResScaleX->setValue(1);
    ui->ResScaleY->setValue(1);

    QObject::connect(ui->ResPlay, SIGNAL(clicked()),
                     this, SLOT(playResVideo()));
    QObject::connect(ui->ResPause, SIGNAL(clicked()),
                     this, SLOT(pauseResVideo()));
    QObject::connect(ui->ResStop, SIGNAL(clicked()),
                     this, SLOT(stopResVideo()));

    QObject::connect(&genVideoT, SIGNAL(nextImage()),
                     this, SLOT(renderResVideo()));
    QObject::connect(&genVideoT, SIGNAL(finishVideo()),
                     this, SLOT(stopResVideo()));

    QObject::connect(this, SIGNAL(nextResVideo()),
                     &genVideoT, SLOT(nextVideoImage()));

    QObject::connect(ui->saveResVideo, SIGNAL(clicked()),
                     this, SLOT(saveResVideo()));


    QObject::connect(uiSaveVideo, SIGNAL(closed()),
                     this, SLOT(resetSaveVideo()));

    ui->ResPlay->setEnabled(false);
    ui->ResPause->setEnabled(false);
    ui->ResStop->setEnabled(false);
    ui->saveResVideo->setEnabled(false);

    mloadResVideo = false;
    this->fileVideo = 0;
    videoStream = 0;
    NFrames = 0;

    QObject::connect(ui->loadResVideo, SIGNAL(clicked()),
                     this, SLOT(loadResVideo()));

    ui->loadResVideo->setEnabled(true);
}

void MainWindow::resetSaveVideo() {
    ui->runSimulationButton->setEnabled(true);
    ui->saveResVideo->setEnabled(true);
}

void MainWindow::saveResVideo() {
    ui->runSimulationButton->setEnabled(false);
    ui->saveResVideo->setEnabled(false);
    uiSaveVideo->saveVideoT->ResImgs = this->simT->resImgs;
    uiSaveVideo->pwd = this->pwd;
    uiSaveVideo->show();

    uiSaveVideo->exec();
}

void MainWindow::loadResVideo() {
    fileName =  QFileDialog::getOpenFileName(
                this,
                "Open simulation data",
                pwd,
                "*.rpmv"
                );
    if (fileName != "") {

        if (fileVideo != 0) {
            fileVideo->close();
            delete fileVideo;
            fileVideo = 0;
        }

        if (videoStream != 0) {
            delete videoStream;
        }


        fileVideo = new QFile(fileName);
        fileVideo->open(QIODevice::ReadOnly);
        videoStream = new QDataStream(fileVideo);

        this->mloadResVideo = true;

        *videoStream >> NFrames;

        ip = 0;
        currImg = 0;

        ////qDebug() << "NFrames = " <<NFrames <<"\n";

        ui->nextResImage->setEnabled(true);
        ui->prevResImage->setEnabled(false);

        ui->ResPlay->setEnabled(true);
        ui->ResPause->setEnabled(false);
        ui->ResStop->setEnabled(false);
    }
}

void MainWindow::playResVideo() {
    mRes.lock();


    genVideoT.N = simT->resImgs->size();
    genVideoT.dt = 1.0 / ui->ResFPS->value();

    ui->loadResVideo->setEnabled(false);
    ui->nextResImage->setEnabled(false);
    ui->prevResImage->setEnabled(false);
    ui->ResPlay->setEnabled(false);
    ui->ResFPS->setEnabled(false);
    ui->ResPause->setEnabled(true);
    ui->ResStop->setEnabled(true);
    if (genVideoT.isPaused()) {
        genVideoT.replay();
    } else {
        currImg = 0;
        genVideoT.init();
    }
    mRes.unlock();

}

void MainWindow::pauseResVideo(){
    mRes.lock();

    ui->ResPlay->setEnabled(true);

    ui->ResFPS->setEnabled(true);

    ui->nextResImage->setEnabled(true);
    ui->prevResImage->setEnabled(true);
    ui->ResPause->setEnabled(false);

    genVideoT.pause();
    mRes.unlock();

}

void MainWindow::stopResVideo(){

    ////qDebug() << "STOP\n";
    genVideoT.stop();

    mRes.lock();

    ui->ResFPS->setEnabled(true);


    ui->nextResImage->setEnabled(true);
    ui->prevResImage->setEnabled(false);

    ui->ResPlay->setEnabled(true);
    ui->ResPause->setEnabled(false);
    ui->ResStop->setEnabled(false);
    ui->loadResVideo->setEnabled(true);
    currImg = 0;
    mRes.unlock();

}

QPixmap MainWindow::getNextImg() {
    QPixmap img;

    if (mloadResVideo && currImg < NFrames && currImg > -1) {

        if (currImg >= ip) {
            sprintf(buffer,"%s/tmp%dxxx",
                    pwd.toAscii().data(),
                    currImg);
            ////qDebug() << buffer << "\n";

            *videoStream >> img;

            QFile file(buffer);
            file.open(QIODevice::WriteOnly);
            QDataStream out(&file);

            out << img;

            ip++;

            file.close();
            ////qDebug() << buffer << "\n";
        } else {
            sprintf(buffer,"%s/tmp%dxxx",
                    pwd.toAscii().data(),
                    currImg);

            ////qDebug() << buffer << "\n";

            QFile file(buffer);
            file.open(QIODevice::ReadOnly);
            QDataStream in(&file);

            in >> img;

            file.close();

            ////qDebug() << buffer << "\n";
        }

    } else img.load(simT->resImgs->at(currImg));

    return img;
}

void MainWindow::renderResVideo(){
    mRes.lock();

    ////qDebug() <<">>>>" << currImg <<"\n";
    ////qDebug() << "NFrames = " <<NFrames <<"\n";

    if ( (simT->resImgs->size() < 1 && !mloadResVideo)
             ||
             (mloadResVideo && NFrames < 1)) {
        currImg = 0;
        emit nextResVideo();
        mRes.unlock();

        stopResVideo();

        return;
    }



    if (ui->resultView->scene() == 0){
        ui->resultView->setScene(new QGraphicsScene);
    }

    int t = static_cast<int>(currImg * genVideoT.getTime());
    int th = t /3600;
    int tm = (t/60)%60;
    int ts = t%60;

    sprintf(buffer,"%02d:%02d:%02d",th,tm,ts);
    ui->resName->setText(buffer);
    ui->resultView->scene()->clear();

    QPixmap img = getNextImg();

    int w = resScaleX * resScale * img.width();
    int h = resScaleY * resScale * img.height();

    ui->resultView->scene()->addPixmap(img.scaled(w,h));
    ui->resultView->scene()->setSceneRect(0.0,0.0,w,h);
    ui->resultView->show();

    currImg++;

    ////qDebug() <<">>>>" << currImg <<"\n";
    ////qDebug() << "NFrames = " <<NFrames <<"\n";


    if ( (currImg >= simT->resImgs->size()  && !mloadResVideo)
            ||
            (mloadResVideo && currImg >= NFrames)){
        currImg = 0;

        emit nextResVideo();
        mRes.unlock();

        stopResVideo();
        return;
    }

    emit nextResVideo();

    mRes.unlock();
}

void MainWindow::scaleResImage(double d) {
    resScale = d;
    this->drawResImage();
}

void MainWindow::scaleResImageX(double d) {
    resScaleX = d;
    this->drawResImage();
}

void MainWindow::scaleResImageY(double d) {
    resScaleY = d;
    this->drawResImage();
}

void MainWindow::nextResImage() {
    currImg++;

    if (currImg > 0) ui->prevResImage->setEnabled(true);

    if ( (currImg >= simT->resImgs->size() - 1 && !mloadResVideo)
          ||
            (mloadResVideo && (currImg >= NFrames - 1))) {
        ui->nextResImage->setEnabled(false);
    }

    if ( (currImg >= simT->resImgs->size() && !mloadResVideo)
         ||
            (mloadResVideo && currImg >= NFrames)){
        if (!mloadResVideo)
            currImg = simT->resImgs->size() - 1;
        else
            currImg = NFrames - 1;
        ui->nextResImage->setEnabled(false);
    }else {
        drawResImage();
    }
}

void MainWindow::prevResImage() {
    currImg--;
    if ( (currImg < simT->resImgs->size() && !mloadResVideo)
         ||
            (mloadResVideo && currImg < NFrames ))
        ui->nextResImage->setEnabled(true);

    if (currImg == 0) {
        ui->prevResImage->setEnabled(false);
    }

    if (currImg <= -1 ){
        currImg = 0;
        ui->prevResImage->setEnabled(false);
    }else {
        drawResImage();
    }
}

void MainWindow::drawResImage() {

    mRes.lock();
    if ( (simT->resImgs->size() < 1 && !mloadResVideo)
            ||
            (mloadResVideo && NFrames < 1)) return;

    if (ui->resultView->scene() == 0){
        ui->resultView->setScene(new QGraphicsScene);
    }
    QPixmap img = getNextImg();

    if (!mloadResVideo)
        ui->resName->setText(simT->resImgs->at(currImg));
    ui->resultView->scene()->clear();




    int w = resScaleX * resScale * img.width();
    int h = resScaleY * resScale * img.height();

    ui->resultView->scene()->addPixmap(img.scaled(w,h));
    ui->resultView->scene()->setSceneRect(0.0,0.0,w,h);
    ui->resultView->show();

    mRes.unlock();

}

void MainWindow::printHelpDom() {
    ui->domGenLog->setText(
                "Help\n"
                "___________________________\n"
                "Metric\n"
                "___________________________\n"
                "0 - Euclidian\n"
                "1 - Manhattan\n"
                "2 - Max\n"
                "____________________________\n"
                "RBF\n"
                "____________________________\n"
                "0 - Gaussian\n"
                );
}

void MainWindow::printHelpSim() {

    ui->simLog->setText(
                "Help\n"
                "___________________________\n"
                "Source\n"
                "___________________________\n"
                "0 - Sin\n"
                "1 - Cos\n"
                "2 - Gaussian\n"
                "3 - Sigmoid\n"
                "4 - Monocicle gaussian\n"
                "5 - Atmospheric pulse\n"
                "6 - Use source data file\n\n"
                "____________________________\n"
                "Mode\n"
                "____________________________\n"
                "0 - Hard\n"
                "1 - Soft\n"
                "2 - Type 3\n"
                );
}

void MainWindow::init()
{
    setPWD();
}


void MainWindow::setInputNameSim() {
    ui->inputSim->setText(
                QFileDialog::getOpenFileName(
                    this,
                    "Open simulation data",
                    pwd,
                    "*.rpim"
                    )
                );
}

void MainWindow::setInputNameDom() {
    ui->inputGen->setText(
                QFileDialog::getOpenFileName(
                    this,
                    "Open simulation data",
                    pwd,
                    "*.rpim"
                    )
                );
}

void MainWindow::setSourceData(){
    ui->sourceData->setText(
                QFileDialog::getOpenFileName(
                    this,
                    "Open source data",
                    pwd,
                    "*.dat"
                    )
                );
}

void MainWindow::setPWD() {
    QString title("Lane Maxwell 2 ");
    pwd = QFileDialog:: getExistingDirectory(this,
                                             "Open a project",
                                             pwd);

    setWindowFilePath(pwd);

#if _WIN32
    pwd.replace("\\","/");
#endif

    title += " [" + pwd +"] ";
    setWindowTitle(title);

    ////qDebug() << title;
}

void MainWindow::enableGen() {
    genT.free();
    ui->runGeneration->setEnabled(true);
    ui->runSimulationButton->setEnabled(true);
}

void MainWindow::enableSim() {
    simT->free();
    ui->runSimulationButton->setEnabled(true);
    ui->runGeneration->setEnabled(true);
    if (simT->videoGen && simT->status) {
        currImg = 0;
        mloadResVideo = false;
        ui->nextResImage->setEnabled(true);
        ui->prevResImage->setEnabled(false);

        ui->loadResVideo->setEnabled(true);
        ui->saveResVideo->setEnabled(true);
        ui->ResPlay->setEnabled(true);
        ui->ResPause->setEnabled(false);
        ui->ResStop->setEnabled(false);

        drawResImage();
    }
}

void MainWindow::stopSimulation() {
    ui->stopSim->setEnabled(false);
    simT->stop();
}

void MainWindow::stopGeneration() {
    ui->stopGen->setEnabled(false);
    genT.stop();
}

void MainWindow::setGenLog(QString s) {
    ui->domGenLog->append(s);
}

void MainWindow::setSimLog(QString s) {
    ui->simLog->append(s);
}

void MainWindow::setAuto() {
    ui->beta0->setEnabled(ui->isBetaAuto->isChecked());
    ui->beta1->setEnabled(ui->isBetaAuto->isChecked());
    ui->NP->setEnabled(ui->isBetaAuto->isChecked());
    ui->errDomain->setEnabled(ui->isBetaAuto->isChecked());

    ui->beta->setEnabled(!ui->isBetaAuto->isChecked());
}

void MainWindow::setUseRatio() {
    ui->ratioDomain->setEnabled(ui->useRatio->isChecked());

    ui->domainSize->setEnabled(!ui->useRatio->isChecked());
}

void MainWindow::runSimulation() {
    simT->resImgs->clear();
    if (ui->createVideo->isChecked()) {
        currImg = 0;
        ui->loadResVideo->setEnabled(false);
        ui->ResPlay->setEnabled(false);
        ui->nextResImage->setEnabled(false);
        ui->prevResImage->setEnabled(false);
    }
    ui->simLog->clear();
    ui->simLog->setText("Running simulation...\n");
    ui->runSimulationButton->setEnabled(false);
    ui->runGeneration->setEnabled(false);

    double C = 299792458.0e0;

    simT->freq = ui->sourceFrequency->value() * 1e6;

    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;
    double lmb = C / (simT->freq);

    simT->resX = ui->resX->value();
    simT->resY = ui->resY->value();
    simT->dt = 0.01 * ui->courantC->value() * lmb * sqrt(mu0 * eps0) / ui->discret->value();
    simT->cte = simT->dt / (2* eps0);
    sprintf(simT->data, "%s", ui->inputSim->text().toAscii().data());
    simT->pico = ui->sourcePick->value();
    simT->isLog = ui->isLog->isChecked();
    simT->videoGen = ui->createVideo->isChecked();
    simT->isTMZ = ui->TMz->isChecked();
    simT->lecont = ui->contrast->value();
    simT->nf = ui->df->value();
    simT->NTHREADS = ui->threadNumberSim->value();
    simT->w = 2 * acos(-1) * simT->freq;
    simT->duracaoPulso = ui->sourceDuration->value() * simT->dt;
    simT->sourceType = ui->sourceMode->value();
    simT->ttype = ui->sourceType->value();
    if (ui->sourceType->value() == USER) {
        FILE* f = fopen (ui->sourceData->text().toAscii().data(),
                         "r");
        if (!f) {
            ui->simLog->setText("Failed to open " +
                                ui->sourceData->text().toAscii()
                                + "\n");
            return;
        }

        while (!feof(f)) {
            double t, x;
            fscanf(f, "%lf%lf", &t, &x);
            simT->fonte.push_back(x);
        }
    }

    simT->Ts = ui->simulationTime->value() * simT->dt;
    simT->pwd = this->pwd;
    if (simT->load()){
        simT->start();
    }
    ui->stopSim->setEnabled(true);
}

void MainWindow::runDomainGeneration() {
    ui->domGenLog->clear();
    ui->domGenLog->setText("Domain generation...\n");
    ui->runSimulationButton->setEnabled(false);
    ui->runGeneration->setEnabled(false);
    genT.beta = ui->beta->value();
    genT.useR = ui->useRatio->isChecked();
    genT.betamin = ui->beta0->value();
    genT.betamax = ui->beta1->value();

    if (genT.betamax < genT.betamin) {
        ui->domGenLog->append("Error: beta0 > beta1\n");
        return;
    }

    genT.metric = ui->metric->value();
    genT.radialBase = ui->RBF->value();
    genT.errbeta = ui->errDomain->value()*1e-6;
    genT.fmax = ui->frequency->value()*1e6;
    genT.domainSize = ui->domainSize->value();
    genT.isAuto = ui->isBetaAuto->isChecked();
    sprintf(genT.data,"%s", ui->inputGen->text().toAscii().data());
    genT.NP = ui->NP->value();
    genT.NTHREADS = ui->threadNumberDG->value();
    genT.ratio = ui->ratioDomain->value();
    lane_maxwell::modules::matrix::matrixModuleInit(1e-11);
    switch (ui->metric->value()) {
    case 0:
        genT.D = lane_maxwell::modules::point::distE;
    break;

    case 1:
        genT.D = lane_maxwell::modules::point::distM;
    break;

    case 2:
        genT.D = lane_maxwell::modules::point::distMax;
    break;
    }

    switch (ui->RBF->value()) {
    case 0:
        genT.R = lane_maxwell::modules::rpim::phi;
        genT.dR = lane_maxwell::modules::rpim::dphi;
    break;
    }

    if (genT.load()) {
        genT.start();
    }
    ui->stopGen->setEnabled(true);

}

MainWindow::~MainWindow()
{
    QProcess call_p;

#ifdef _WIN32
    call_p.execute("del "+pwd+"/tmp*");
#else
    call_p.execute("rm "+pwd+"/tmp*");
#endif

    if (ui) {
        delete ui;
    }
    if (simT) {
        delete simT;
    }

    if (uiSaveVideo) {
        delete uiSaveVideo;
    }

    if (fileVideo) {
        fileVideo->close();
        delete fileVideo;
    }

    if (videoStream) {
        delete videoStream;
    }
}
