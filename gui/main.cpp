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

#include <QLocale>
#include <clocale>
#include <QtGui/QApplication>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    setlocale(LC_ALL,"C");
    QLocale::setDefault(QLocale::C);


    QApplication a(argc, argv);
    MainWindow w;

    w.init();
    w.show();
    a.exec();


    return 0;
}
