#GoGoRPIM - A RPIM implementation and 2-D Eletromagnetic simulator
#    Copyright (C) 2012  Péricles Lopes Machado (LANE-UFPA)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

function test_rpim_1(tipo)
    fmax = 1e9;
    lmb = 3e8 / fmax;
    dx = dy = lmb / 20;

    if tipo == 1
        f = fopen("pontos_igual.dat");
        fp = fopen("dist_1.dat","w+");
    end

    if tipo == 2
        f = fopen("pontos_igual.dat");
        fp = fopen("dist_2.dat","w+");    
    end

    if tipo == 3
        f = fopen("pontos_igual.dat");
        fp = fopen("dist_3.dat","w+");
    end

    m = 1;
    NR = 3;

    while (!feof(f))    
        [px, py, count] = fscanf(f,"%lf%lf", '2');

        if count == 0
            break;
        end

        x(m) = px * dx;
        y(m) = py * dy;

        
        if tipo == 3 && m > 1       
            x(m) = (1 + 0.8 * rand() - 0.4) * x(m);
            y(m) = (1 + 0.7 * rand() - 0.37) * y(m); 
        end

        if tipo == 2 && m > 1  && NR > 0     
            x(m) = (1 + 3 * rand() - 1.5) * x(m);
            y(m) = (1 + 3 * rand() - 1.5) * y(m);
            NR = NR - 1; 
        end


        fprintf(fp,"%e %e\n", x(m), y(m));       

        m = m+1;

    end

    fclose(fp);
    fclose(f);


    pts = [x' y'];

    N = size(pts)(1);
    M = size(pts)(2);
    pi = acos(-1);

    K = 2 * pi * fmax / 3e8;

    F = [];
    dF = [];

    for n = 2:N
        F(n-1) = 0;
        F(n-1) = F(n-1) + sin(K*pts(n,1)) + cos(K*pts(n, 2));        
        %F(n-1) = sin(K * F(n-1));
    end

    Fr = 1;


    dFr(1) = K; 
    dFr(2) = 0;

    if tipo == 1
        f = fopen("err_fx_1.dat", "w+");
        f2 = fopen("err_fy_1.dat", "w+");
        f3 = fopen("err_f_1.dat", "w+");
    end

    if tipo == 2
        f = fopen("err_fx_2.dat", "w+");
        f2 = fopen("err_fy_2.dat", "w+");
        f3 = fopen("err_f_2.dat", "w+");
    end

    if tipo == 3
        f = fopen("err_fx_3.dat", "w+");
        f2 = fopen("err_fy_3.dat", "w+");
        f3 = fopen("err_f_3.dat", "w+");
    end
    
    errant1 = -100;
    errant2 = -100;
    errant3 = -100;

    nio = 0;
    for c=0.01:0.5:40
        
        if c > -0.1 && c < 0.1
            continue;        
        end

        [phi, dphi] = rpim(pts, c);
        Fi = irpim(dphi(1,:), F);
        if (abs(dFr(1)) > 1e-5)
            fprintf(f, "%e %e\n", c, err=100*(Fi-dFr(1))/dFr(1));
        else
            fprintf(f, "%e %e\n", c, err=Fi-dFr(1));
        end

        
        if( (err < 0 && errant1 > 0) || 
            (err > 0 && errant1 < 0)
            )
            %printf("c1 = %e err %e\n", c, err);
        end

        errant1 = err;

        Fi = irpim(dphi(2,:), F);
        if (abs(dFr(2)) > 1e-5)
            fprintf(f2, "%e %e\n", c, err=100*(Fi-dFr(2))/dFr(2));
        else
            fprintf(f2, "%e %e\n", c, err=Fi-dFr(2));
        end

        
        if( (err < 0 && errant2 > 0) || 
            (err > 0 && errant2 < 0)
            )
           % printf("c2 = %e err %e\n", c, err);
        end

        errant2 = err;

        Fi = irpim(phi, F);
        if (abs(Fr) > 1e-5)
            fprintf(f3, "%e %e\n", c, err=100*(Fi-Fr)/Fr);
        else
            fprintf(f3, "%e %e\n", c, err=Fi-Fr);
        end

        if( (err < 0 && errant3 > 0) || 
            (err > 0 && errant3 < 0)
            )
            printf("c3 = %e err %e\n", c, err);
        end

        errant3 = err;
        nio = nio + 1;
    end



    fclose(f);
    fclose(f2);
    fclose(f3);

    errmax = 1e-10;
    NP = 25;

    cmin = 1;
    cmax = 100;

    nio

    [c, err, ni] = lsfcm(pts, fmax, cmin, cmax, NP, errmax);

    c
    err
    ni


    K


    [phi, dphi] = rpim(pts,c(1));


    dFr
    Fr

    irpim(phi,F)
    irpim(dphi(1,:),F)
    irpim(dphi(2,:),F)
end


