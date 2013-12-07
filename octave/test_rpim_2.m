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

function test_rpim_2(tipo)
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

        xo(m) = x(m) = px * dx;
        yo(m) = y(m) = py * dy;

        
        if tipo == 3 && m > 1       
            x(m) = (1 + 0.2 * rand() - 0.1) * x(m);
            y(m) = (1 + 0.2 * rand() - 0.1) * y(m); 
        end

        if tipo == 2 && m > 1  && NR > 0     
            x(m) = (1 + 0.2 * rand() - 0.1) * x(m);
            y(m) = (1 + 0.2 * rand() - 0.1) * y(m);
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


    Fr = 1;


    dFr(1) = K; 
    dFr(2) = 0;

    

    errmax = 1e-10;
    NP = 24;

    cmin = 1;
    cmax = 101;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tipo == 3
        M = 7;
    else
        M = 1;
    end
    
    


    if tipo == 1
        f11 = fopen("err_1lsfcm_1_1.dat","w+");
        f12 = fopen("err_1lsfcm_2_1.dat","w+");

        f21 = fopen("err_2lsfcm_1_1.dat","w+");
        f22 = fopen("err_2lsfcm_2_1.dat","w+");

        f31 = fopen("err_3lsfcm_1_1.dat","w+");
        f32 = fopen("err_3lsfcm_2_1.dat","w+");
    end

    if tipo == 2
        f11 = fopen("err_1lsfcm_1_2.dat","w+");
        f12 = fopen("err_1lsfcm_2_2.dat","w+");

        f21 = fopen("err_2lsfcm_1_2.dat","w+");
        f22 = fopen("err_2lsfcm_2_2.dat","w+");

        f31 = fopen("err_3lsfcm_1_2.dat","w+");
        f32 = fopen("err_3lsfcm_2_2.dat","w+");
    end

    if tipo == 3
        f11 = fopen("err_1lsfcm_1_3.dat","w+");
        f12 = fopen("err_1lsfcm_2_3.dat","w+");

        f21 = fopen("err_2lsfcm_1_3.dat","w+");
        f22 = fopen("err_2lsfcm_2_3.dat","w+");

        f31 = fopen("err_3lsfcm_1_3.dat","w+");
        f32 = fopen("err_3lsfcm_2_3.dat","w+");
    end

    printf("lsfcm\n");

    for NP = 2:2:16
        errm(1, NP) = 0;
        errm2(1, NP) = 0;

        errm(2, NP) = 0;
        errm2(2, NP) = 0;

        errm(3, NP) = 0;
        errm2(3, NP) = 0;
    end

    for NI = 1:M

        if tipo == 3
            for m = 2:N       
                x(m) = xo(m)+(0.02 * rand() - 0.01) * dx;
                y(m) = yo(m)+(0.02 * rand() - 0.01) * dy; 
            end
        end



        pts = [x' y'];

        F = [];
        dF = [];

        for n = 2:N
            F(n-1) = 0;
            F(n-1) = F(n-1) + sin(K*pts(n,1)) + cos(K*pts(n, 2));        
        end

        for  NP = 2:2:16

            printf("NP = %d\n", NP);

            cmin = 1;
            cmax = 10;
            [c, err, ni] = lsfcm(pts, fmax, cmin, cmax, NP, errmax);

            %[phi, dphi] = rpim(pts, c(1));
            %Fi = irpim(phi, F);
            %err(1) = Fi - Fr;

            
            %[phi, dphi] = rpim(pts, c(2));
            %Fi = irpim(dphi(1,:), F);
            %err(2) = Fi - dFr(1);

            
            %[phi, dphi] = rpim(pts, c(3));
            %Fi = irpim(dphi(2,:), F);
            %err(3) = Fi - dFr(2);

            errm(1, NP) = errm(1, NP) + abs(err(1));
            errm(2, NP) = errm(2, NP) + abs(err(2));
            errm(3, NP) = errm(3, NP) + abs(err(3));

            cmin = 1;
            cmax = 10;
            [c, err, ni] = lsfcm2(pts, fmax, cmin, cmax, NP, errmax);
           
            %[phi, dphi] = rpim(pts, c(1));
            %Fi = irpim(phi, F);
            %err(1) = Fi - Fr;

            
            %[phi, dphi] = rpim(pts, c(2));
            %Fi = irpim(dphi(1,:), F);
            %err(2) = Fi - dFr(1);

            
            %[phi, dphi] = rpim(pts, c(3));
            %Fi = irpim(dphi(2,:), F);
            %err(3) = Fi - dFr(2);

            errm2(1, NP) = errm2(1, NP) + abs(err(1));
            errm2(2, NP) = errm2(2, NP) + abs(err(2));
            errm2(3, NP) = errm2(3, NP) + abs(err(3));

        end

         printf("M = %d\n", M);

    end

    for NP = 2:2:16
        fprintf(f11, "%d %e\n", NP, errm(1, NP)/M);
        fprintf(f12, "%d %e\n", NP, errm2(1, NP)/M);

        fprintf(f21, "%d %e\n", NP, errm(2, NP)/M);
        fprintf(f22, "%d %e\n", NP, errm2(2, NP)/M);

        fprintf(f31, "%d %e\n", NP, errm(3, NP)/M);
        fprintf(f32, "%d %e\n", NP, errm2(3, NP)/M);
    end

    fclose(f11); fclose(f12);
    fclose(f21); fclose(f22);
    fclose(f31); fclose(f32);
end


