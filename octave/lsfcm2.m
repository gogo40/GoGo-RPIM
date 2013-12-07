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

function [c, errmin, ni] = lsfcm2(pts, fmax, cmin, cmax, NP, errmax)
    
    N = size(pts)(1);
    M = size(pts)(2);

    pi = acos(-1);    
    vo = 3e8;
    K = 2 * pi * fmax / vo;
K
    dbeta = (cmax - cmin) / NP;

    ni(1) = 0;

    for m = 1:M
        ni(m + 1) = 0;
    end

    F = [];
    dFr = [];

    for n = 2:N
        F(n-1) = 0;
        Fr = 0;
    
        acc = 0;
        for m = 1:M
            if mod(m,2) == 1
                F(n-1) = F(n-1) + sin(K*pts(n,m));
            else
                F(n-1) = F(n-1) + cos(K*pts(n,m));
                Fr = Fr + 1;
            end
            acc = acc + K*pts(n,m);
        end

        F(n-1) = cos(acc);
    end

    Fr = 1;


    for m = 1:M
        if mod(m,2) == 1
            dFr(m) = K;
        else
            dFr(m) = 0;
        end
         dFr(m) = 0;
    end

%
%Otimização de PHI
%
    ok = false;
    bestc = cmin * cmin;
    errmin(1) = 100;

    a = 0;
    b = 0;
    fa = 0;
    fb = 0;

    fib = false;

    for beta=cmin:dbeta:cmax
        xa = beta;
        xb = beta + 0.5 * dbeta;

        [phi, dphi] = rpim(pts, xa * xa);

        Fi = irpim(phi, F);

        dfa = Fi - Fr;

        if !ok || abs(dfa) < errmin(1)
            ok = true;
            errmin(1) = abs(dfa);
            bestc = xa * xa;        
        end

        [phi, dphi] = rpim(pts, xb * xb);

        Fi = irpim(phi, F);

        dfb = Fi - Fr;

        if !ok || abs(dfb) < errmin(1)
            ok = true;
            errmin(1) = abs(dfb);
            bestc = xb * xb;        
        end

        if dfa * dfb < 0 
            fa = dfa;
            fb = dfb;
            a = xa;
            b = xb;
            fib = true;
        end

        ni(1) = ni(1) + 1;
    end

    k = 0;
    fx = fa;
    fib;
   
    while k < NP && errmin(1) > errmax
        x = (a*fb - b*fa)/(fb-fa);

        [phi, dphi] = rpim(pts, x * x);

        Fi = irpim(phi, F);
        fxo = fx;
        fx = Fi - Fr;
        
        if  fx*fa > 0
            a =  x;
            fa = fx;
            if fx*fxo > 0
                fb = 0.5*fb;
            end
        else
            b = x;
            fb = fx;
            if fx*fxo > 0
                fa = 0.5*fa;            
            end
        end

        if abs(fx) < errmin(1)
            errmin(1) = abs(fx);
            bestc = x * x;
        end

        ni(1) = ni(1) + 1;
        k = k + 1;                
    end

    c(1) = bestc;

%
%Otimização de DPHI
%
    for m=1:M
        ok = false;
        bestc = cmin;
        errmin(m+1) = 100;

        a = 0;
        b = 0;
        fa = 0;
        fb = 0;
        ni(m+1) = 0;
        fib = false;
        for beta=cmin:dbeta:cmax
            xa = beta;
            xb = beta + 0.5 * dbeta;

            [phi, dphi] = rpim(pts, xa * xa);

            Fi = irpim(dphi(m,:), F);

            dfa = Fi - dFr(m);

            if !ok || abs(dfa) < errmin(m+1)
                ok = true;
                errmin(m+1) = abs(dfa);
                bestc = xa * xa;        
            end

            [phi, dphi] = rpim(pts, xb * xb);

            Fi = irpim(dphi(m,:), F);

            dfb = Fi - dFr(m);

            if !ok || abs(dfb) < errmin(m+1)
                ok = true;
                errmin(m+1) = abs(dfb);
                bestc = xb * xb;        
            end

            if dfa * dfb < 0 
                fa = dfa;
                fb = dfb;
                a = xa;
                b = xb;
                fib = true;
            end

            ni(m+1) = ni(m+1) + 1;
        end

        k = 0;
        fx = fa;
        fib;

        while k < NP && errmin(m+1) > errmax
            x = (a*fb - b*fa)/(fb-fa);

            [phi, dphi] = rpim(pts, x * x);

            Fi = irpim(dphi(m,:), F);
            fxo = fx;
            fx = Fi - dFr(m);
            
            if  fx*fa > 0
                a =  x;
                fa = fx;
                if fx*fxo > 0
                    fb = 0.5*fb;
                end
            else
                b = x;
                fb = fx;
                if fx*fxo > 0
                    fa = 0.5*fa;            
                end
            end

            if abs(fx) < errmin(m+1)
                errmin(m+1) = abs(fx);
                bestc = x * x;
            end

            ni(m+1) = ni(m+1) + 1;
            k = k + 1;                
        end

        c(m+1) = bestc;

    end
end



