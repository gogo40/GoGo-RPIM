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


#ifndef M_PI
const double M_PI = acos(-1);
#endif



int type;
double t_ant=-1;
double pico;
double w;
std::vector<double> fonte;
FILE* fpf;

double source(const double& t, const double& dt,
    const lane_maxwell::modules::point::Point* pts) {

    int i;

    double lamb, C, neper;
    double tau, t0;
    double f = w / (2 * M_PI);
    double sinal = 0;

    double pico1=2.0e-3;
    double tau1=400e-12;
    double tau2=2.5e-9;
    double nn=2.0e+0;

    double eta = exp(  - (tau1/tau2)*pow(  nn*tau2/tau1 ,1.0e0/nn )  );

    C = 299792458.0e0;
    neper = 2.718281e0;

    lamb = C / f;

    tau = 50.0e0*dt  , t0 = 130.0e0*dt  ;

    switch (type) {
        case SIN:
            sinal = pico * std::sin(w * t);
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case COS:
            sinal = pico * std::cos(w * t);
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case GAUSS:
            sinal = pico * std::sin(w * t) * std::exp(-(t - t0) * (t - t0) / tau / tau); //gaussiana modulada em seno
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case SIGMOID:
            sinal = pico * std::sin(w * t) / (1 + std::exp(5 - (t / tau))); //sigmoid modulada em seno
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case MONOCICLOGAUSS:
            //tau = 10.0e0*dt, t0 = 75.0e0*dt;
            tau = 40.0e0*dt, t0 = 150.0e0*dt;
            sinal = pico * std::sqrt(2*neper/(tau*tau)) * (t-t0) * std::exp(-std::pow((t-t0),2)/(tau*tau)); //geracao do monociclogaussiano.
            if(t > t_ant) { // Evita repetição do mesmo valos para fontes não pontuais
                t_ant = t;
                fprintf(fpf,"%e %e\n", t_ant, sinal);
            }
            return sinal;
        break;

        case USER:
            i = static_cast<int>(round(t / dt));
            sinal = pico * fonte[i];
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;

        case DESCARGAATM:
            sinal = pico*(((pico1/eta)/(pico1/eta)*(exp(-t/tau2))*(pow((t/tau1),1.0e0/nn)/(pow((t/tau1),nn)+1.0e0)))/0.52231);
            fprintf(fpf,"%e %e\n", t, sinal);
            return sinal;
        break;
    }

    return 0;
}

bool RunSimulation::load()
{
    this->resImgs->clear();
    this->source = ::source;
    buffer.clear();
    QTextStream s(&buffer);
    char ss[1000];
    sprintf(ss,"%s/fonte.dat", pwd.toAscii().data());
    ::fpf = fopen(ss,"w+");

    if (!::fpf) {
        fclose(::fpf);
        s << "Failed to open fonte.dat!\n";
        emit setLog(s.readAll());
        emit setEnabled();
        return false;
    }

    ::fonte = this->fonte;
    ::w = this->w;
    ::type = this->ttype;
    ::pico = this->pico;

    double eps0 = 8.854187817620e-12;
    double mu0 = 4 * acos(-1) * 1e-7;

    s << "Mode = " << ((isTMZ)?"TMZ\n":"TEZ\n");
    s << "mu0 = " << mu0 << "\n";
    s << "eps0 = " << eps0 <<"\n";
    s << "data = " << this->data << "\n";
    s << "NTHREADS = " << this->NTHREADS << "\n";
    s << "nf = " << this->nf << "\n";
    s << "dt = " << this->dt << "\n";
    s << "cte = " << this->cte << "\n";
    s << "freq = " << this->freq << "\n";
    s << "w = " << this->w << "\n";
    s << "Ts = " << this->Ts << "\n";
    s << "Source mode = ";
    switch (this->sourceType) {
    case 0: s << "HARD\n"; break;
    case 1: s << "SOFT\n"; break;
    case 2: s << "TYPE3\n"; break;
    }
    s << "Source pike = " << this->pico <<"\n";
    s << "Source duration = " << this->duracaoPulso << "\n";
    s << "Source type = ";

    switch (::type) {
    case SIN: s <<"SIN\n"; break;
    case COS: s<< "COS\n"; break;
    case GAUSS: s << "GAUSS\n"; break;
    case SIGMOID: s << "SIGMOID\n"; break;
    case MONOCICLOGAUSS: s << "MONOCICLOGAUSS\n"; break;
    case DESCARGAATM:  s << "DESCARGAATM\n"; break;
    case USER: s << "USER\n"; break;
    }

    if (this->videoGen) {
        s << "Create video = YES\n";
        s << "Res X = " << this->resX << "\n";
        s << "Res Y = " << this->resY << "\n";
        s << "Use logarithmic scale = " << this->isLog << "\n";
        s << "nf = " << this->nf << "\n";
        s << "lecont = " << this->lecont << "\n";
    } else s << "Create video = NO\n";



    threads.resize(NTHREADS);

    for (size_t i = 0; i < NTHREADS; ++i) {
        threads[i] = new MaxwellT;
        threads[i]->data = this;
    }


    setlocale(LC_ALL,"C");

    return true;

}


void RunSimulation::run()
{
    setlocale(LC_ALL,"C");
    this->online();
    QTextStream s(&buffer);
    s << "Loading points...\n";
    emit setLog(s.readAll());


    if (!lane_maxwell::modules::rpim::load(this->data, &this->pts)) {
        s << "Failed to open " << this->data <<"\n";
        fclose(::fpf);
        emit setEnabled();
        emit setLog(s.readLine());
        return;
    }
    s << "Points loaded!\n";
    emit setLog(s.readLine());

    QTime t;
    QProcess call_p;

    t.start();
    if (this->isTMZ) {
        int nt = static_cast<int>(Ts / dt);
        Files fE, fH[2];
        point::Points vPE, vPH;
        point::Points vSource;
        point::Points::iterator it = pts->begin(), end = pts->end();
        double Lx, Ly;
        double xmax, xmin;
        char fname[1000];
        char com[1000];
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
                if (type == 'E') {
                    std::sprintf(fname,"%s/Ez_%g_%g.dat",
                                 pwd.toAscii().data(),
                                 x, y);
                    fE.push_back(std::fopen(fname,"w+"));
                    vPE.push_back(p);
                    std::sprintf(fname,"%s/Ez_%g_%g.png",
                                 pwd.toAscii().data(),
                                 x, y);
                    resImgs->push_back(fname);

                } else {
                    std::sprintf(fname,"%s/Hx_%g_%g.dat",
                                 pwd.toAscii().data(),
                                 x, y);
                    fH[0].push_back(std::fopen(fname,"w+"));

                    std::sprintf(fname,"%s/Hx_%g_%g.png",
                                 pwd.toAscii().data(),
                                 x, y);
                    resImgs->push_back(fname);

                    std::sprintf(fname,"%s/Hy_%g_%g.dat",
                                 pwd.toAscii().data(),
                                 x, y);
                    fH[1].push_back(std::fopen(fname,"w+"));

                    std::sprintf(fname,"%s/Hy_%g_%g.png",
                                 pwd.toAscii().data(),
                                 x, y);
                    resImgs->push_back(fname);

                    vPH.push_back(p);
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

        lane_maxwell::modules::maxwell::initTMZ(pts);

        s << "L = " <<Lx <<" X " << Ly <<"\n";
        emit setLog(s.readLine());

        int N = pts->size();
        int fs = N / NTHREADS;

        for (size_t i = 0; i < NTHREADS; ++i) {
            threads[i]->ini = i * fs;
            threads[i]->fim = min((i+1)*fs, N);
            threads[i]->isTMZ = true;
        }

        for (int t = 0; t < nt; ++t) {
            s << t+1 <<"/" << nt <<"\n";
            emit setLog(s.readLine());

            if (this->isStoped()) {
                emit setLog("Stopped!");
                break;
            }

            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->updateE = true;
                threads[i]->start();
            }

            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->wait();
            }

            point::Points::iterator it = vSource.begin(),
                end = vSource.end();

            while (it != end) {
                point::Point* p = *it;
                int prop = p->getProp();

                if (sourceType == TYPE3) {
                    if (!(prop & P_FONTEE)) {
                        if (t * dt < duracaoPulso)
                            p->setE(0, source(t * dt, dt, p));
                    } else {
                        if (prop & P_DIRFX) {
                            if (t * dt < duracaoPulso)
                                p->setH(0, source(t*dt, dt, p));
                        } else {
                            if (t * dt < duracaoPulso)
                                p->setH(1, source(t*dt, dt, p));
                        }
                    }
                } else if (sourceType == HARD) {
                    if (!(prop & P_FONTEE)) {
                        p->setE(0, source(t * dt, dt, p));
                    } else {
                        if (prop & P_DIRFX) {
                            p->setH(0, source(t*dt, dt, p));
                        } else {
                            p->setH(1, source(t*dt, dt, p));
                        }
                    }
                } else {
                    if (!(prop & P_FONTEE)) {
                        double Ez = p->getE(0);
                        p->setE(0, Ez + source(t * dt, dt, p));
                    } else {
                        if (prop & P_DIRFX) {
                            double Hx = p->getH(0);
                            p->setH(0, Hx + source(t*dt, dt, p));
                        } else {
                            double Hy = p->getH(1);
                            p->setH(1, Hy + source(t*dt, dt, p));
                        }
                    }
                }

                ++it;
            }

            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->updateE = false;
                threads[i]->start();
            }

            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->wait();
            }

            for (size_t i = 0; i < vPE.size(); ++i) {
                point::Point* p = vPE[i];
                std::FILE* f = fE[i];
                std::fprintf(f, "%e %e\n", t * dt, p->getE(0));
            }

            for (size_t i = 0; i < vPH.size(); ++i) {
                point::Point* p = vPH[i];
                std::FILE* f = fH[0][i];

                std::fprintf(f, "%e %e\n", t * dt, p->getH(0));

                f = fH[1][i];
                std::fprintf(f, "%e %e\n", t * dt, p->getH(1));
            }

            //Video generation
            if (videoGen) {
                if (t % nf == 0) {
                    double fmax = 0, emax = 0;
                    double fmin = 1e20, emin = 1e20;
                    it = pts->begin();
                    end = pts->end();

                    for (size_t ux = 0; ux < nx; ++ux) {
                        for (size_t uy = 0; uy < ny; ++uy) {
                            image[ux][uy] = 20;
                        }
                    }

                    while (it != end) {
                        point::Point* p = *it;
                        double x = p->getX(0);
                        double y = p->getX(1);
                        int prop = p->getProp();
                        char type = p->getType();
                        if (type == 'E') {
                            size_t ux = static_cast<int>(resX * (x - xmin) / Lx);
                            size_t uy = static_cast<int>(resY * (y - ymin) / Ly);
                            double Ez = fabs(p->getE(0));

                            if (!(prop & P_UPML)) {
                                if (Ez > fmax) fmax = Ez + 1;
                                if (Ez < fmin) fmin = Ez;
                                if (p->getEps() > emax) emax = p->getEps() + 1;
                                if (p->getEps() < emin) emin = p->getEps();
                            }
                        }
                        ++it;
                    }

                    double eps0 = 8.841941e-12;
                    it = pts->begin();
                    while (it != end) {
                        point::Point* p = *it;
                        double x = p->getX(0);
                        double y = p->getX(1);
                        int prop = p->getProp();
                        char type = p->getType();
                        if (type == 'E') {
                            size_t ux = static_cast<int>(resX * (x - xmin) / Lx);
                            size_t uy = static_cast<int>(resY * (y - ymin) / Ly);
                            double Ez = fabs(p->getE(0));

                            if (!(prop & P_UPML)) {
                                if(isLog) {
                                    image[ux][uy] = Ez * (Ez - fmin) / (fmax - fmin);
                                    if (p->getEps() > 2 * eps0)
                                        image[ux][uy] += 5 * Ez + 2 * Ez * (p->getEps() - emin) / (emax - emin);
                                    //image[ux][uy] = Ez;
                                } else {
                                    image[ux][uy] = Ez * (p->getE(0) - fmin) / (fmax - fmin);
                                    if (p->getEps() > 2 * eps0)
                                        image[ux][uy] +=  5 * Ez + 2 * Ez * (p->getEps() - emin) / (emax - emin);
                                    //image[ux][uy] = p->getE(0);
                                }
                            }
                        }
                        ++it;
                    }

                    //interpolImage(image, nx, ny);

                    std::sprintf(fname, "%s/ez-%d",
                                 pwd.toAscii().data(),
                                 ni);

                    std::sprintf(com,"t=%e |Ez|=[%e %e] eps = [%e %e]",t * dt,fmin, fmax - 1, emin, emax);
                    levelCurve(nx, ny, image, fname, com, isLog, 1, 0);


                    std::sprintf(com,"convert -level %c%d"
                    " %s.png %s.png",'%',lecont,fname, fname);


#ifdef _WIN32
                    call_p.execute(com);
#else
                    system(com);
#endif
                    emit setLog(com);

                    std::sprintf(fname, "%s/ez-%d.png",
                                 pwd.toAscii().data(),
                                 ni);
                    resImgs->push_back(fname);


                    ++ni;
                }
            }
        }

#ifndef _WIN32
        if (videoGen) {
            sprintf(fname,"mencoder 'mf://%s/ez-%%d.png' -mf fps=10"
            " -o %s/EZ.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800",
                   pwd.toAscii().data(),
                   pwd.toAscii().data());
            system(fname);
        }
#endif

        for (size_t i = 0; i < vPE.size(); ++i) {
            double x = vPE.at(i)->getX(0);
            double y = vPE.at(i)->getX(1);

            std::FILE* f = fE[i];
            fclose(f);

            std::sprintf(fname,"%s/Ez_%g_%g.dat",
                         pwd.toAscii().data(),
                         x, y);

            std::FILE* fgpl = fopen("gen.gpl","w+");
            std::fprintf(fgpl, "set term png\n"
                         "set out '%s/Ez_%g_%g.png'\n"
                         "set xlabel 't'\n"
                         "set ylabel 'Ez'\n"
                         "plot '%s' notitle w l\n",
                         pwd.toAscii().data(), x, y,
                         fname);

            std::fclose(fgpl);

            call_p.execute("gnuplot gen.gpl");
            videoGen = true;
        }

        for (size_t i = 0; i < vPH.size(); ++i) {
            std::FILE* f = fH[0][i];
            std::fclose(f);

            double x = vPH.at(i)->getX(0);
            double y = vPH.at(i)->getX(1);

            std::sprintf(fname,"%s/Hx_%g_%g.dat",
                         pwd.toAscii().data(),
                         x, y);

            std::FILE* fgpl = fopen("gen.gpl","w+");
            std::fprintf(fgpl, "set term png\n"
                         "set out '%s/Hx_%g_%g.png'\n"
                         "set xlabel 't'\n"
                         "set ylabel 'Hx'\n"
                         "plot '%s' notitle w l\n",
                         pwd.toAscii().data(), x, y,
                         fname);

            std::fclose(fgpl);

            call_p.execute("gnuplot gen.gpl");

            f = fH[1][i];
            std::fclose(f);

            std::sprintf(fname,"%s/Hy_%g_%g.dat",
                         pwd.toAscii().data(),
                         x, y);

            fgpl = fopen("gen.gpl","w+");
            std::fprintf(fgpl, "set term png\n"
                         "set out '%s/Hy_%g_%g.png'\n"
                         "set xlabel 't'\n"
                         "set ylabel 'Hy'\n"
                         "plot '%s' notitle w l\n",
                         pwd.toAscii().data(), x, y,
                         fname);

            std::fclose(fgpl);

            call_p.execute("gnuplot gen.gpl");
            videoGen = true;
        }

        for (size_t i = 0; i < nx; ++i) {
            std::free(image[i]);
        }
        std::free(image);
        s << "Simulation finished!\n";
    } else {
        int nt = static_cast<int>(Ts / dt);
        Files fE[2], fH;
        point::Points vPE, vPH;
        point::Points vSource;
        char fname[1000], com[1000];
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
                    std::sprintf(fname,"%s/Hz_%g_%g.dat",
                                  pwd.toAscii().data(),
                                 x, y);
                    fH.push_back(std::fopen(fname,"w+"));
                    vPH.push_back(p);


                    std::sprintf(fname,"%s/Hz_%g_%g.png",
                                 pwd.toAscii().data(),
                                 x, y);
                    resImgs->push_back(fname);
                } else {
                    std::sprintf(fname,"%s/Ex_%g_%g.dat",
                                  pwd.toAscii().data(),
                                 x, y);
                    fE[0].push_back(std::fopen(fname,"w+"));


                    std::sprintf(fname,"%s/Ex_%g_%g.png",
                                 pwd.toAscii().data(),
                                 x, y);
                    resImgs->push_back(fname);

                    std::sprintf(fname,"%s/Ey_%g_%g.dat",
                                  pwd.toAscii().data(),
                                 x, y);
                    fE[1].push_back(std::fopen(fname,"w+"));


                    std::sprintf(fname,"%s/Ey_%g_%g.png",
                                 pwd.toAscii().data(),
                                 x, y);
                    resImgs->push_back(fname);

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


        s << "L = " <<Lx <<" X " << Ly <<"\n";
        emit setLog(s.readAll());

        int N = pts->size();
        int fs = N / NTHREADS;

        for (size_t i = 0; i < NTHREADS; ++i) {
            threads[i]->ini = i * fs;
            threads[i]->fim = min((i+1)*fs, N);
            threads[i]->isTMZ = false;
        }

        for (int t = 0; t < nt; ++t) {

            if (this->isStoped()) {
                break;
            }

            s << t+1 <<"/" << nt <<"\n";
            emit setLog(s.readAll());

            //Update Hz
            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->updateE = false;
                threads[i]->start();
            }

            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->wait();
            }

            point::Points::iterator it = vSource.begin(),
                end = vSource.end();

            while (it != end) {
                point::Point* p = *it;
                int prop = p->getProp();

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
            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->updateE = true;
                threads[i]->start();
            }

            for (size_t i = 0; i < NTHREADS; ++i) {
                threads[i]->wait();
            }

            for (size_t i = 0; i < vPH.size(); ++i) {
                point::Point* p = vPH[i];
                std::FILE* f = fH[i];
                std::fprintf(f, "%e %e\n", t * dt, p->getH(0));
            }

            for (size_t i = 0; i < vPE.size(); ++i) {
                point::Point* p = vPE[i];
                std::FILE* f = fE[0][i];

                std::fprintf(f, "%e %e\n", t * dt, p->getE(0));

                f = fE[1][i];
                std::fprintf(f, "%e %e\n", t * dt, p->getE(1));
            }

            //Video generation
            if (videoGen) {
                if (t % nf == 0) {
                    double fmax = 0, emax = 0;
                    double fmin = 1e20, emin = 1e20;
                    it = pts->begin();
                    end = pts->end();


                    for (size_t ux = 0; ux < nx; ++ux) {
                        for (size_t uy = 0; uy < ny; ++uy) {
                            image[ux][uy] = 20;
                        }
                    }

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
                                if (Hz > fmax) fmax = Hz +1;
                                if (Hz < fmin) fmin = Hz;
                                if (p->getEps() > emax) emax = p->getEps() + 1;
                                if (p->getEps() < emin) emin = p->getEps() + 1;
                            }
                        }
                        ++it;
                    }


                    double eps0 = 8.841941e-12;
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
                                    image[ux][uy] = Hz * (Hz - fmin) / (fmax - fmin);
                                    if (p->getEps() > 2 * eps0)
                                        image[ux][uy] += 5 * Hz + Hz * 2 * (p->getEps() - emin) / (emax - emin);
                                    //image[ux][uy] = Hz;
                                } else {
                                    image[ux][uy] = Hz * (p->getH(0) - fmin) / (fmax - fmin);
                                    if (p->getEps() > 2 * eps0)
                                        image[ux][uy] += 5 * Hz +  Hz * 2 * (p->getEps() - emin) / (emax - emin);
                                    //image[ux][uy] = p->getH(0);
                                }
                            }
                        }
                        ++it;
                    }

                    //interpolImage(image, nx, ny);
                    std::sprintf(fname, "%s/hz-%d",
                                 pwd.toAscii().data(),
                                 ni);


                    std::sprintf(com,"t=%e |Hz|=[%e %e] eps = [%e %e]",t * dt,fmin,fmax - 1, emin, emax);
                    levelCurve(nx, ny, image, fname, com, isLog, 1, 0);

                    std::sprintf(com,"convert -level %c%d"
                    " %s.png %s.png",'%',lecont,fname, fname);

#ifdef _WIN32
                    call_p.execute(com);
#else
                    system(com);
#endif
                    emit setLog(com);

                    std::sprintf(fname, "%s/hz-%d.png",
                                 pwd.toAscii().data(),
                                 ni);
                    resImgs->push_back(fname);

                    ++ni;
                }
            }
        }
#ifndef _WIN32
        if (videoGen) {
            sprintf(fname,"mencoder 'mf://%s/hz-%%d.png' -mf fps=10"
            " -o %s/HZ.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800",
                   pwd.toAscii().data(),
                   pwd.toAscii().data());
            system(fname);
        }
#endif

        for (size_t i = 0; i < vPH.size(); ++i) {
            std::FILE* f = fH[i];
            std::fclose(f);

            double x = vPH.at(i)->getX(0);
            double y = vPH.at(i)->getX(1);

            std::sprintf(fname,"%s/Hz_%g_%g.dat",
                         pwd.toAscii().data(),
                         x, y);

            std::FILE* fgpl = fopen("gen.gpl","w+");
            std::fprintf(fgpl, "set term png\n"
                         "set out '%s/Hz_%g_%g.png'\n"
                         "set xlabel 't'\n"
                         "set ylabel 'Hz'\n"
                         "plot '%s' notitle w l\n",
                         pwd.toAscii().data(), x, y,
                         fname);

            std::fclose(fgpl);

            call_p.execute("gnuplot gen.gpl");
            videoGen = true;

        }

        for (size_t i = 0; i < vPE.size(); ++i) {
            std::FILE* f = fE[0][i];
            std::fclose(f);

            double x = vPE.at(i)->getX(0);
            double y = vPE.at(i)->getX(1);

            std::sprintf(fname,"%s/Ex_%g_%g.dat",
                         pwd.toAscii().data(),
                         x, y);

            std::FILE* fgpl = fopen("gen.gpl","w+");
            std::fprintf(fgpl, "set term png\n"
                         "set out '%s/Ex_%g_%g.png'\n"
                         "set xlabel 't'\n"
                         "set ylabel 'Ex'\n"
                         "plot '%s' notitle w l\n",
                         pwd.toAscii().data(), x, y,
                         fname);

            std::fclose(fgpl);

            call_p.execute("gnuplot gen.gpl");

            f = fE[1][i];
            std::fclose(f);

            std::sprintf(fname,"%s/Ey_%g_%g.dat",
                         pwd.toAscii().data(),
                         x, y);

            fgpl = fopen("gen.gpl","w+");
            std::fprintf(fgpl, "set term png\n"
                         "set out '%s/Ey_%g_%g.png'\n"
                         "set xlabel 't'\n"
                         "set ylabel 'Ey'\n"
                         "plot '%s' notitle w l\n",
                         pwd.toAscii().data(), x, y,
                         fname);

            std::fclose(fgpl);

            call_p.execute("gnuplot gen.gpl");
            videoGen = true;

        }

        for (size_t i = 0; i < nx; ++i) {
            std::free(image[i]);
        }
        std::free(image);
        s << "Simulation finished!\n";
    }

    fclose(::fpf);

    emit setLog(s.readLine());

    s << "Iime elapsed = " << t.elapsed() / 1000.0 << " s\n";

    emit setLog(s.readLine());

    emit setEnabled();
    this->status = true;
}


