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

/*
GoGoNeuro - A RPIM implementation and 2-D Eletromagnetic simulator
    Copyright (C) 2011  Péricles Lopes Machado (LANE-UFPA)

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

#include <cmath>
#include <cstdio>
#include <vector>

using namespace std;


/*
sinTAXE:
I ou R = IFouri ou Rfouri (funcao,2.*pi*f,dt,nt*dt,nt);
(ponto da transformada na freq. f Hz)
*/

/* Parte Imaginaria */

double Ifouri(const vector<double>& f,double t,double a,double b,int m){
	double n;
	double h,th,sum,f1,f2,s1,s2,c1,c2,a1,a2,sinth2,costh2,al,be,ga,fouri,su1,su;
	int i;
	
	n=m-1;
	h=(b-a)/(n);
	th=t*h;
	sum=0.;
	f1=f[1];
	f2=f[m];
	s1=sin(a*t);
	s2=sin(b*t);
	c1=cos(a*t);
	c2=cos(b*t);
	a1=a+h;
	sinth2=(sin(th))*(sin(th));
	costh2=(cos(th))*(cos(th));
	al=(th*th+th*sin(th)*cos(th)-2.*sinth2)/(th*th*th);
	be=2.*(th*(1.+costh2)-2.*sin(th)*cos(th))/(th*th*th);
	ga=4.*(sin(th)-th*cos(th))/(th*th*th);

	su =f2*s2-f1*s1;
	su1=-0.5*(f2*c2-f1*c1);
	
	for (i=2;i<=n;i+=2){
		sum=sum+f[i]*cos(a1*t);
		a1 =a1+h;
		su1=su1+f[i+1]*cos(a1*t);
		a1 =a1+h;
	}

	fouri=h*(al*su+be*su1+ga*sum)/(b-a);
	return (fouri);
}

/* Parte Real */

double Rfouri(const vector<double>& f,double t,double a,double b,int m){
	double n;
	double h,th,sum,f1,f2,s1,s2,c1,c2,a1,a2,sinth2,costh2,al,be,ga,fouri,su1,su;
	int i;
	
	n=m-1;
	h=(b-a)/(n);
	th=t*h;
	sum=0.;
	f1=f[1];
	f2=f[m];
	s1=sin(a*t);
	s2=sin(b*t);
	c1=cos(a*t);
	c2=cos(b*t);
	a1=a+h;
	sinth2=(sin(th))*(sin(th));
	costh2=(cos(th))*(cos(th));
	al=(th*th+th*sin(th)*cos(th)-2.*sinth2)/(th*th*th);
	be=2.*(th*(1.+costh2)-2.*sin(th)*cos(th))/(th*th*th);
	ga=4.*(sin(th)-th*cos(th))/(th*th*th);

	su =f1*c1-f2*c2;
	su1=-0.5*(f2*s2-f1*s1);
	
	for (i=2;i<=n;i+=2){
		sum=sum+f[i]*sin(a1*t);
		a1 =a1+h;
		su1=su1+f[i+1]*sin(a1*t);
		a1 =a1+h;
	}

	fouri=h*(al*su+be*su1+ga*sum)/(b-a);
	return (fouri);
}

int main(int argc, char** argv)
{
	if (argc != 4) {
		printf("%s <frequencia> <ang> <data>\n", argv[0]);
		return 0;
	}

	FILE* f = fopen(argv[3], "r");

	double freq, ang;
	vector<double> v;

	sscanf(argv[1], "%lf", &freq);
    sscanf(argv[2], "%lf", &ang);

	double t1, t2;
	int c = 0;
	while (!feof(f)) {
		double t, x;
		fscanf(f, "%lf%lf", &t, &x);

		if (fabs(x) > 1e-7) {
			if (c < 2){
				t1 = t2;
				t2 = t;
			}
		
			c++;
			v.push_back(x);
		}
	}

	int nt = v.size()-1;
	double dt = t2 - t1;
	double df = (2.5 * freq) / nt;


	double iF = Ifouri(v, 2 * M_PI * freq, dt, nt * dt, nt);
	double rF = Rfouri(v, 2 * M_PI * freq, dt, nt * dt, nt);

	double z = sqrt(iF*iF + rF*rF);
	iF = Ifouri(v, -2 * M_PI * freq, dt, nt * dt, nt);
	rF = Rfouri(v, -2 * M_PI * freq, dt, nt * dt, nt);
	z+= sqrt(iF*iF + rF*rF);
	printf("%e %e\n", ang, z);
	
	return 0;
}


