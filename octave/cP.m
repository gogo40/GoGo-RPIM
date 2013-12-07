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

function P=cP(x)

N = size(x)(1);
M = size(x)(2);
N 
M
for n = 2:N
    P(n-1, 1) = 1;
    for m = 1:M
        P(n-1, m+1) = x(n,m);
    end
end

endfunction
