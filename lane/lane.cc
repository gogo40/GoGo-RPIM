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
@file lane.cc
@brief Neste arquivo são declaradas funções genéricas compartilhadas pelos diversos módulos do programa


UFPA - Universidade Federal do Pará
Laboratório de Análises Numéricas em Eletromagnetismo (LANE)

Implementação do RPIM - Radial Point Interpolation Method

Esta implementação faz parte do TCC de Engenharia da Computaçã de Péricles Lopes Machado (0608000101)
*/

#include <float.h>
#include "lane/lane.h"

#ifdef QT_GUI_LIB
#ifdef _WIN32
#define WQTL
#endif
#endif

#ifdef WQTL
#include <QProcess>
#endif

#undef RGB

namespace lane_maxwell {
Mutex printMutex;

static int dx[] = { 0, 0, 1,-1, 1,-1, 1,-1};
static int dy[] = { 1,-1, 0, 0,-1, 1, 1,-1};

float mediaInterpolImage(float **image, int i, int j, int nx, int ny) {
    float v;

    v = 0;

    for (int k = 0; k < 8; ++k) {
        int I = i + dx[k];
        int J = j + dy[k];
        if(I>=0 && I < nx && J >= 0 && J < ny) {
            if (image[I][J] != 0) {
                v += image[I][J];
            }
        }
    }
    return v/8;
}

void interpolImage(float **image, int nx, int ny) {
    int N = nx;
    while (N) {
        --N;
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                if(image[i][j] == 0) {
                    image[i][j] = mediaInterpolImage(image, i, j, nx, ny);
                }
    }
}

/*
Algoritmo para geracao de imagens PPM/PNG
a partir de uma matriz de dados.
    (Curva de Nivel)
    Rodrigo M. S. de Oliveira.
    05-06/04/2008
    atualizado em 15/04/2008

*/
void RGB(float sinal, float a0, float a1, float a2, float a3,
    float a4, float a5, float tn, float htn,
    float *R, float *G, float *B) {
    if (sinal <= a1) {
        *B = (htn / (a1 - a0))*(sinal - a0) + htn;
        *G = 0.0;
        *R = 0.0;
    } else if (a1 < sinal && sinal < a2) {
        *B = tn;
        *G = (tn / (a2 - a1)) * (sinal - a1);
        *R = 0.0;
    } else if (a2 <= sinal && sinal <= a3) {
        *B = (tn / (a2 - a3)) * (sinal - a3);
        *G = tn;
        *R = (tn / (a3 - a2)) * (sinal - a2);

    } else if ( a3 <= sinal && sinal<= a4 ) {
        *B = 0.0;
        *G = (tn / (a3 - a4)) * (sinal - a4);
        *R = tn;
    } else {
        *B = 0.0;
        *G = 0.0;
        *R = (htn / (a4 - a5)) * (sinal - a5) + htn;
    }
}

void levelCurve(
     int im,		/* Limite da Matriz (i)           */
     int jm,		/* Limite da Matriz (j)           */
     float **a,		/* Matriz                         */
     const char* fname,	/* Nome do Arquivo de saida       */
     const char* comment,	/* Texto a ser inserido na imagem */
     int LOG,		/* usar escala log    ?           */
     int PNG,		/* converter em PNG / JPG   ?     */
     int DAT		/* salvar arquivo DAT ?	          */
     ) {
    int BAR_PIXELS = 05, BLACK_PIXELS = 25;
    int i,j;
    float R,G,B;
    char mystr[1000]; // string de proposito geral

#ifdef WQTL
    QProcess call_p;
#endif
    /* Arquivo dat */
    if (DAT) {
        std::sprintf(mystr,"%s.dat",fname);
        std::FILE *ptr_dat = std::fopen(mystr,"w+");
        std::rewind(ptr_dat);

        for (j = 0; j < jm; ++j) {
            for (i = 0; i < im; ++i) {
                std::fprintf(ptr_dat,"%E\t",a[i][j]);
            }
            std::fprintf(ptr_dat,"\n");
        }
        std::fclose(ptr_dat);
    }

    /* valor absoluto */
    for (i = 0; i < im; ++i) {
        for (j = 0; j < jm; ++j) {
            a[i][j] = (a[i][j] < 0.0)? - a[i][j]: a[i][j];//abs
        }
    }

    /* (raiz quinta de a) = "LOG suave"  */
    if (LOG) {
        for (i = 0; i < im; ++i) {
            for (j = 0; j < jm; ++j) {
                a[i][j] = powf(a[i][j], 0.20);
            }
        }
    }

    /* localizacao do maximo e do minimo */
    float amax = -FLT_MAX; //-3.4028235E38; //gcc
    float amin = FLT_MAX; // 3.4028235E38;

    for (i = 0; i < im; ++i) {
        for (j = 0; j < jm; ++j) {
            amax = (a[i][j] < amax)? amax: a[i][j];
            amin = (a[i][j] > amin)? amin: a[i][j];
        }
    }

    float aw = (amax-amin)/4.0; //largura dos intervalos

    /* limites de a -> RGB */
    float	a0 = amin,
        a1 = 0.5 * aw,
        a2 = a1 + aw,
        a3 = a2 + aw,
        a4 = a3 + aw,
        a5 = amax;

    /* Total de niveis para R, G e B */
    float	tn  = 255.0,
        htn = 0.5 * tn;

/* Criacao do arquivo ppm */
    std::sprintf(mystr,"%s.ppm",fname);

    std::FILE *ptr_ppm = std::fopen(mystr, "w+");
    std::rewind(ptr_ppm);

    /* cabecalho PPM */

    std::fprintf(ptr_ppm, "P3\n%d %d\n%d\n", im,
        jm + 2 * BAR_PIXELS+BLACK_PIXELS,static_cast<int>(tn));
    std::fprintf(ptr_ppm, "#Created by LANE MAXWELL\n"); // comentario no arquivo PPM

    /*	Barra escura superior    */
    for (j = 0; j < BLACK_PIXELS; ++j) {
        for (i = 0; i < im; ++i) {
            std::fprintf(ptr_ppm,"0 0 0 ");
        }
    }

    float dcolor = (!LOG)? (amax - amin) /(static_cast<float>(im)):
            (powf(amax, 5.0) - powf(amin, 5.0))/(static_cast<float>(im));
    float color;

    /*  Barra de Cores   */
    for (j = 0; j < BAR_PIXELS; ++j) {
        for (i = 0; i < im; ++i){
            color = ( (float) i )*dcolor ;
            if (LOG) color = powf(color,0.20);

            RGB (color,a0,a1,a2,a3,a4,a5,tn,htn,&R,&G,&B);

            std::fprintf(
            ptr_ppm,"%d %d %d ",
                    static_cast<int>(nearbyintf(R)),
                    static_cast<int>(nearbyintf(G)),
                    static_cast<int>(nearbyintf(B))
            );
        }
    }

    /*	Barra escura inferior    */
    for (j = 0; j < BAR_PIXELS; ++j) { /*  dados p/ imagem PPM   */
        for (i = 0; i < im; ++i) {
            std::fprintf(ptr_ppm, "0 0 0 ");
        }
    }

    /*  conversao da matriz a[i][j] em cores (PPM)*/
    for (j = jm - 1; j >= 0; --j) {
        for (i = 0; i < im; ++i) {
            RGB (a[i][j],a0,a1,a2,a3,a4,a5,tn,htn,&R,&G,&B);
            std::fprintf(
            ptr_ppm,"%d %d %d ",
                    static_cast<int>(nearbyintf(R)),
                    static_cast<int>(nearbyintf(G)),
                    static_cast<int>(nearbyintf(B))
            );
        }
    }

    fclose(ptr_ppm);

    /*
        gera imagem png
        ( ImageMagick )
    */

    if (PNG == 1) {
        if (im < 800 || jm < 800) {
            int c;

            if (im < jm) {
                c = 800 / im;
            } else {
                c = 800 / jm;
            }

            sprintf(mystr, "convert  -resize %dx%d -fill white"
            " -pointsize 15  -draw \"text 15,17 \'%s\'\" %s.ppm %s.png",
            c * im, c * jm, comment, fname, fname);
        } else {
            sprintf(mystr, "convert  -fill white -pointsize 15"
            "  -draw \"text 15,17 \'%s\'\" %s.ppm %s.png",
            comment, fname, fname);
        }
#ifdef WQTL
        call_p.execute(mystr);
#else
        system(mystr);
#endif
    } else if (PNG == 2) {
        if (im < 800 || jm < 800) {
            int c;

            if (im < jm) {
                c = 800 / im;
            } else {
                c = 800 / jm;
            }

            sprintf(mystr,"convert  -resize %dx%d -fill black"
            " -pointsize 15  -draw \"text 15,17 \'%s\'\" %s.ppm %s.png",
            c * im, c * jm, comment, fname, fname);
        }else{
            sprintf(mystr,"convert  -fill black -pointsize 15"
            "  -draw \"text 15,17 \'%s\'\" %s.ppm %s.png",
            comment, fname, fname);
        }
#ifdef WQTL
        call_p.execute(mystr);
#else
        system(mystr);
#endif
    }
}
}  // lane_maxwell
