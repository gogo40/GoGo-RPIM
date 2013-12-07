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

#ifndef SAVEVIDEOFORM_H
#define SAVEVIDEOFORM_H

#include <QDialog>

#include "saveResVideo.h"

namespace Ui {
class saveVideoForm;
}

class saveVideoForm : public QDialog
{
    Q_OBJECT

public:
    explicit saveVideoForm(QWidget *parent = 0);
    ~saveVideoForm();

    SaveResVideoT* saveVideoT;

    QString pwd;

private:
    Ui::saveVideoForm *ui;
    QString fileName;
    QPixmap img;


public slots:
    void sendClose();
    void sendClose(int);
    void getFileName();
    void saveVideo();
    void videoSaved();
    void setLog(QString);
    void loadImg(QString);
signals:
    void closed();
};

#endif // SAVEVIDEOFORM_H
