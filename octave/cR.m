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

function R=cR(x, c)

N=size(x)(1);
rmax = 0;

for n=2:N
    for m=2:N
        r = distE(x(n,:), x(m,:));
        if r > rmax
            rmax = r;
        end
    end    
end


rmax = rmax * rmax;

for n=2:N
    for m=2:N
        r = distE(x(n,:), x(m,:));    
        R(n-1,m-1) = exp(-c*r*r/rmax);
    end
end

endfunction
