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

/**Copyright 2012 Péricles Lopes Machado
@file tez.cc
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa

UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include "maxwell/tez.h"

namespace lane_maxwell { namespace modules { namespace maxwell {

static char fname[1000];
static char com[1000];


namespace maxwelltez {

struct CalcArgs {
    point::Points* pts;
    size_t ini;
    size_t fim;
    double cte;
    double dt;
    bool status;

    CalcArgs(point::Points* pts, size_t ini,
        size_t fim, double cte, double dt)
        : pts(pts), ini(ini), fim(fim), cte(cte),
        dt(dt), status(true) {}

    ~CalcArgs() {}

    DISALLOW_COPY_AND_ASSIGN(CalcArgs);
};

void* threadE(void* p_args) {
    CalcArgs* args = static_cast<CalcArgs*>(p_args);
    point::Points* pts = args->pts;
    size_t ini = args->ini;
    size_t fim = args->fim;
    double cte = args->cte;
    double dt = args->dt;

    updateETEZ(pts, ini, fim, cte, dt);

    pthread_exit(NULL);
}

void* threadH(void* p_args) {
    CalcArgs* args = static_cast<CalcArgs*>(p_args);
    point::Points* pts = args->pts;
    size_t ini = args->ini;
    size_t fim = args->fim;
    double cte = args->cte;
    double dt = args->dt;

    updateHTEZ(pts, ini, fim, cte, dt);

    pthread_exit(NULL);
}

bool calcEH(bool isH, point::Points* pts, size_t NTHREADS, const double& dt, const double& cte) {
    size_t N = pts->size();
    std::vector<CalcArgs*> args(NTHREADS);
    std::vector<pthread_t> threads(NTHREADS);

    size_t sT = N / NTHREADS;
    int rc;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

    for (size_t i = 0; i < NTHREADS; ++i) {

        size_t ini = i * sT, fim = (i + 1) * sT;
        fim = min(N, fim);

        args[i] = new CalcArgs(pts, ini, fim, cte, dt);
    }

    for (size_t i=0; i < NTHREADS; ++i) {
        if (isH) {
            pthread_create(&threads[i], &attr, threadH,
            static_cast<void*>(args[i]));
        } else {
            pthread_create(&threads[i], &attr, threadE,
            static_cast<void*>(args[i]));
        }
    }

    pthread_attr_destroy(&attr);
    void* status;

    for(size_t i = 0; i < NTHREADS; ++i){
        rc = pthread_join(threads[i], &status);
        if (rc) {
            fprintf(stderr, "ERROR; return code from"
            "pthread_join() is %d\n", rc);
            return false;
        }
    }

    bool ok = true;

    for (size_t i = 0; i < NTHREADS; ++i) {
        if (!args[i]->status) {
            ok = false;
        }
        delete args[i];
    }

    args.clear();
    threads.clear();

    return ok;

}

}  // maxwelltez


bool tez(point::Points* pts, const double& Ts, const double& dt,
    const double& cte, double resX, double resY,
    bool videoGen, int nf, bool isLog, int lecont,
    size_t NTHREADS, PFSource source, int sourceType,
    double duracaoPulso) {
    int nt = static_cast<int>(Ts / dt);
    Files fE[2], fH, fP[2];
    point::Points vPE, vPH;
    point::Points vSource;
    point::Points::iterator it = pts->begin(), end = pts->end();
    double Lx, Ly;
    double xmax, xmin;
    double ymax, ymin;

    while (it != end) {
        point::Point *p = *it;
        double x = p->getX(0);
        double y = p->getX(1);
        int prop = p->getProp();
        char type = p->getType();

        if (it == pts->begin()) {
            xmax = xmin = x;
            ymax = ymin = y;
        } else {
            if (x > xmax) xmax = x;
            if (y > ymax) ymax = y;

            if (x < xmin) xmin = x;
            if (y < ymin) ymin = y;
        }

        if (prop & P_MED) {
            if (type == 'H') {
                std::sprintf(fname,"Hz_%g_%g.dat", x, y);
                fH.push_back(std::fopen(fname,"w+"));

                std::sprintf(fname,"Potx_%g_%g.dat", x, y);
                fP[0].push_back(std::fopen(fname,"w+"));

                std::sprintf(fname,"Poty_%g_%g.dat", x, y);
                fP[1].push_back(std::fopen(fname,"w+"));

                vPH.push_back(p);
            } else {
                std::sprintf(fname,"Ex_%g_%g.dat", x, y);
                fE[0].push_back(std::fopen(fname,"w+"));

                std::sprintf(fname,"Ey_%g_%g.dat", x, y);
                fE[1].push_back(std::fopen(fname,"w+"));

                std::sprintf(fname,"Potx_%g_%g.dat", x, y);
                fP[0].push_back(std::fopen(fname,"w+"));

                std::sprintf(fname,"Poty_%g_%g.dat", x, y);
                fP[1].push_back(std::fopen(fname,"w+"));


                vPE.push_back(p);
            }
        }

        if (prop & P_FONTE) {
            vSource.push_back(p);
        }

        ++it;
    }

    Lx = xmax - xmin;
    Ly = ymax - ymin;

    int ni = 0;
    size_t nx = static_cast<int>(resX + 2);
    size_t ny = static_cast<int>(resY + 2);
    float** image = static_cast<float**>(std::malloc(nx * sizeof(float*)));

    for (size_t i = 0; i < nx; ++i) {
        image[i] = static_cast<float*>(std::malloc(ny * sizeof(float)));
    }

    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            image[i][j] = 0;
        }
    }

    initTEZ(pts);

    printf("L = %e x %e\n", Lx, Ly);

    for (int t = 0; t < nt; ++t) {

        printf("%d/%d\n", t+1, nt);

        //Update Hz
        maxwelltez::calcEH(true, pts, NTHREADS, dt, cte);

        point::Points::iterator it = vSource.begin(),
            end = vSource.end();

        while (it != end) {
            point::Point* p = *it;
            int prop = p->getProp();

            if (prop & P_METAL || prop & P_PEC) {
                ++it; continue;
            }

            if (sourceType == TYPE3) {
                if ((prop & P_FONTEE)) {
                    if (t * dt < duracaoPulso)
                        p->setH(0, source(t * dt, dt, p));
                } else {
                    if (prop & P_DIRFX) {
                        if (t * dt < duracaoPulso)
                            p->setE(0, source(t*dt, dt, p));
                    } else {
                        if (t * dt < duracaoPulso)
                            p->setE(1, source(t*dt, dt, p));
                    }
                }
            } else if (sourceType == HARD) {
                if ((prop & P_FONTEE)) {
                    p->setH(0, source(t * dt, dt, p));
                } else {
                    if (prop & P_DIRFX) {
                        p->setE(0, source(t*dt, dt, p));
                    } else {
                        p->setE(1, source(t*dt, dt, p));
                    }
                }
            } else {
                if ((prop & P_FONTEE)) {
                    double Hz = p->getH(0);
                    p->setH(0, Hz + source(t * dt, dt, p));
                } else {
                    if (prop & P_DIRFX) {
                        double Ex = p->getE(0);
                        p->setE(0, Ex + source(t*dt, dt, p));
                    } else {
                        double Ey = p->getE(1);
                        p->setE(1, Ey + source(t*dt, dt, p));
                    }
                }
            }


            ++it;
        }

        //Update Ex and Ey
        maxwelltez::calcEH(false, pts, NTHREADS, dt, cte);

        for (size_t i = 0; i < vPH.size(); ++i) {
            point::Point* p = vPH[i];
            std::FILE* f = fH[i];
            std::fprintf(f, "%e %e\n", t * dt, p->getH(0));

            f = fP[0][i];
            std::fprintf(f, "%e %e\n", t * dt, p->getH(0) * p->calcE(1));

            f = fP[1][i];
            std::fprintf(f, "%e %e\n", t * dt, p->getH(0) * p->calcE(0));
        }

        for (size_t i = 0; i < vPE.size(); ++i) {
            point::Point* p = vPE[i];
            std::FILE* f = fE[0][i];

            std::fprintf(f, "%e %e\n", t * dt, p->getE(0));

            f = fE[1][i];
            std::fprintf(f, "%e %e\n", t * dt, p->getE(1));

            f = fP[0][i];
            std::fprintf(f, "%e %e\n", t * dt, p->getE(1) * p->calcH(0));

            f = fP[1][i];
            std::fprintf(f, "%e %e\n", t * dt, p->getE(0) * p->calcH(0));
        }

        //Video generation
        if (videoGen) {
            if (t % nf == 0) {
                double emax = 0;
                double emin = 1e20;
                it = pts->begin();
                end = pts->end();


                for (size_t ux = 0; ux < nx; ++ux) {
                    for (size_t uy = 0; uy < ny; ++uy) {
                        image[ux][uy] = 20;
                    }
                }

                it = pts->begin();
                while (it != end) {
                    point::Point* p = *it;
                    double x = p->getX(0);
                    double y = p->getX(1);
                    int prop = p->getProp();
                    char type = p->getType();
                    if (type == 'H') {
                        size_t ux = static_cast<int>(resX * (x - xmin) / Lx);
                        size_t uy = static_cast<int>(resY * (y - ymin) / Ly);
                        double Hz = fabs(p->getH(0));

                        if (!(prop & P_UPML)) {
                            if(isLog) {
                                image[ux][uy] = Hz;
                            } else {
                                image[ux][uy] = p->getH(0);
                            }

                            if (Hz > emax) emax = Hz;
                            if (Hz < emin) emin = Hz;
                        }
                    }
                    ++it;
                }

                //interpolImage(image, nx, ny);

                std::sprintf(fname, "hz-%d", ni);
                std::sprintf(com,"t=%e |Hz|=[%e %e] eps = [%e %e]",t * dt, emin, emax);
                levelCurve(nx, ny, image, fname, com, isLog, 1, 0);
                std::sprintf(fname,"convert -level %c%d"
                " hz-%d.png hz-%d.png",'%',lecont,ni,ni);
                std::system(fname);
                std::system("*.ppm");
                ++ni;
            }
        }
    }

    if (videoGen) {
        std::system ("mencoder 'mf://hz-%d.png' -mf fps=10"
        " -o HZ.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800");
    }

    for (size_t i = 0; i < vPH.size(); ++i) {
        std::FILE* f = fH[i];
        std::fclose(f);

        f = fP[0][i];
        std::fclose(f);

        f = fP[1][i];
        std::fclose(f);
    }

    for (size_t i = 0; i < vPE.size(); ++i) {
        std::FILE* f = fE[0][i];
        std::fclose(f);

        f = fE[1][i];
        std::fclose(f);

        f = fP[0][i];
        std::fclose(f);

        f = fP[1][i];
        std::fclose(f);
    }

    for (size_t i = 0; i < nx; ++i) {
        std::free(image[i]);
    }
    std::free(image);

    return true;
}
}}}  // lane_maxwell::modules::tez

