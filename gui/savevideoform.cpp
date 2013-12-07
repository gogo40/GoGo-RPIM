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

#include "savevideoform.h"
#include "ui_savevideoform.h"

saveVideoForm::saveVideoForm(QWidget *parent) :
    QDialog(parent),
    saveVideoT(new SaveResVideoT),
    pwd(""),
    ui(new Ui::saveVideoForm),
    fileName("movie.rpmv")
{

    ui->setupUi(this);
    QObject::connect(ui->close, SIGNAL(clicked()),
                     this, SLOT(sendClose()));
    QObject::connect(this, SIGNAL(finished(int)),
                     this, SLOT(sendClose(int)));
    QObject::connect(ui->file, SIGNAL(clicked()),
                     this, SLOT(getFileName()));
    QObject::connect(ui->save, SIGNAL(clicked()),
                     this, SLOT(saveVideo()));

    QObject::connect(saveVideoT, SIGNAL(terminated()),
                     this, SLOT(videoSaved()));

    QObject::connect(saveVideoT, SIGNAL(sendLog(QString)),
                     this, SLOT(setLog(QString)));
    QObject::connect(saveVideoT, SIGNAL(loadImg(QString)),
                     this, SLOT(loadImg(QString)));
    ui->save->setEnabled(false);
}

void saveVideoForm::loadImg(QString f) {
    ////qDebug() <<"1 Loading..." <<f <<"\n";
    img.load(f);
    ////qDebug() <<" 2 Loading..." <<f <<"\n";
    saveVideoT->setLoaded(true);
}

void saveVideoForm::setLog(QString t) {
    ui->log->setText(t);
}

void saveVideoForm::videoSaved() {
    ui->save->setEnabled(true);
    ui->close->setEnabled(true);
}

void saveVideoForm::saveVideo() {
    ui->log->clear();
    saveVideoT->fileName = &fileName;

    ui->save->setEnabled(false);
    ui->close->setEnabled(false);
    saveVideoT->img = &this->img;

    saveVideoT->start();

}

void saveVideoForm::getFileName() {

    fileName =  QFileDialog::getSaveFileName(
                this, "Name movie file",pwd,"*.rpmv");
    if (fileName != "") {
        ui->save->setEnabled(true);
        this->ui->fileName->setText(fileName);
    }
}

void saveVideoForm::sendClose(int a) {
    this->close();
    emit closed();
}

void saveVideoForm::sendClose()
{
    this->close();
    emit closed();
}

saveVideoForm::~saveVideoForm()
{
    delete saveVideoT;
    delete ui;
}
