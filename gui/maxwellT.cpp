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

#include "runSimulation.h"

void MaxwellT::run() {
    point::Points* pts = data->pts;
    double cte = data->cte;
    double dt = data->dt;

    if (updateE) {
        if (isTMZ) {
            lane_maxwell::modules::maxwell::updateETMZ(pts,
                                                   ini,
                                                   fim,
                                                   cte,
                                                   dt);
        } else {
            lane_maxwell::modules::maxwell::updateETEZ(pts,
                                                   ini,
                                                   fim,
                                                   cte,
                                                   dt);

        }
    } else {
        if (isTMZ) {
            lane_maxwell::modules::maxwell::updateHTMZ(pts,
                                                  ini,
                                                   fim,
                                                   cte,
                                                   dt);
        } else {
            lane_maxwell::modules::maxwell::updateHTEZ(pts,
                                                  ini,
                                                   fim,
                                                   cte,
                                                   dt);
        }
    }
}


