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

function [phi, dphi] = rpim(pts, c)
N = size(pts)(1);
M = size(pts)(2);

PxT(1) = 1;
for n = 1:M
    PxT(n+1) = pts(1, n);
end

rmax = 0;
for n = 2:N
    r = distE(pts(n,:), pts(1,:));
    if r > rmax
        rmax = r;    
    end
end

rmax = rmax * rmax;

for n = 2:N
    r = distE(pts(n,:), pts(1,:));
    RxT(n-1) = exp(-c * r * r / rmax);
end


for n = 2:N
    for m = 2:N
        r = distE(pts(n,:), pts(m,:));
        Ro(n-1, m-1) = exp(-c * r * r / rmax);
    end
end


for n = 2:N
    Po(n-1,1) = 1;
    for m = 1:M
        Po(n-1, m+1) = pts(n, m);
    end
end

Sb = inv(Po' * inv(Ro) * Po) * Po' * inv(Ro);
Sa = (inv(Ro) - inv(Ro) * Po * Sb);

phi = RxT * Sa + PxT * Sb;

for m = 1:M
    for n = 2:N    
        dphi(m, n-1) = 0;

        for k = 2:N
            r = distE(pts(k,:),pts(1,:));
            dRxT(k-1) = -2 * c * (pts(k, m) - pts(1, m)) * exp(-c * r * r / rmax);
            dphi(m, n-1) = dphi(m, n-1)  * Sa(k-1, n-1);
        end

        dphi(m, n-1) = dphi(m, n-1) + Sb(m+1, n-1);
        
    end
end


endfunction 


