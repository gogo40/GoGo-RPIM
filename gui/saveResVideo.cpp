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

#include "saveResVideo.h"

void SaveResVideoT::init(){

}

void SaveResVideoT::run(){

    QString buffer;
    QTextStream log(&buffer);

    QFile file(*this->fileName);
    file.open(QIODevice::WriteOnly);

    if (!file.isWritable()) {
        log << "Failed to save!\n";
        emit sendLog(log.readAll());
        emit sendLog("Finished!");
        emit terminated();
        return;
    }
    QDataStream out(&file);


    int L = ResImgs->size();

    out << L;

    log << "Starting generation\n";
    emit sendLog(log.readAll());

    for (int i = 0; i < L; ++i) {

        setLoaded(false);

        emit loadImg(ResImgs->at(i));

        while(true) {
            if (imgLoaded()) break;
            msleep(10);
        }

        log <<"File: " << ResImgs->at(i) <<"\n";
        emit sendLog(log.readAll());

        out << *img;


        log << "Saving frame "<<i << "/" << L <<"...\n";
        emit sendLog(log.readAll());
    }

    file.close();

    emit sendLog("Finished!");
    emit terminated();
}

