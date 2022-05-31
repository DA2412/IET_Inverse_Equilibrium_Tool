% IET (Inverse Equilibrium Tool) is a MATLAB computational tool constituted by three main modules,
% each with a dedicated GUI. It allows to compute the coil currents needed to obtain a predetermined 
% plasma shape with defined plasma global parameters (i.e. total plasma current and total poloidal 
% magnetic flux at the boundary) by solving a constrained minimization problem. 
% 
%     Copyright (C) 2019 Domenico Abate - domenico.abate@igi.cnr.it
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
function [xx,ww] = funGaussLeg(nn)
switch nn
    case 1
        xx(1)=+0.000000000000000e+00;
        ww(1)=+2.000000000000000e+00;
    case 2
        xx(1)=-5.773502691896260e-01;
        xx(2)=+5.773502691896260e-01;
        ww(1)=+1.000000000000000e+00;
        ww(2)=+1.000000000000000e+00;
    case 3
        xx(1)=-7.745966692414834e-01;
        xx(2)=+0.000000000000000e+00;
        xx(3)=+7.745966692414834e-01;
        ww(1)=+5.555555555555530e-01;
        ww(2)=+8.888888888888890e-01;
        ww(3)=+5.555555555555530e-01;
    case 4
        xx(1)=-8.611363115940530e-01;
        xx(2)=-3.399810435848563e-01;
        xx(3)=+3.399810435848563e-01;
        xx(4)=+8.611363115940530e-01;
        ww(1)=+3.478548451374480e-01;
        ww(2)=+6.521451548625464e-01;
        ww(3)=+6.521451548625464e-01;
        ww(4)=+3.478548451374480e-01;
    case 5
        xx(1)=-9.061798459386640e-01;
        xx(2)=-5.384693101056830e-01;
        xx(3)=+0.000000000000000e+00;
        xx(4)=+5.384693101056830e-01;
        xx(5)=+9.061798459386640e-01;
        ww(1)=+2.369268850561820e-01;
        ww(2)=+4.786286704993664e-01;
        ww(3)=+5.688888888888890e-01;
        ww(4)=+4.786286704993664e-01;
        ww(5)=+2.369268850561820e-01;
    case 6
        xx(1)=-9.324695142031520e-01;
        xx(2)=-6.612093864662650e-01;
        xx(3)=-2.386191860831970e-01;
        xx(4)=+2.386191860831970e-01;
        xx(5)=+6.612093864662650e-01;
        xx(6)=+9.324695142031520e-01;
        ww(1)=+1.713244923791623e-01;
        ww(2)=+3.607615730481390e-01;
        ww(3)=+4.679139345726900e-01;
        ww(4)=+4.679139345726900e-01;
        ww(5)=+3.607615730481390e-01;
        ww(6)=+1.713244923791623e-01;
    case 7
        xx(1)=-9.491079123427590e-01;
        xx(2)=-7.415311855993950e-01;
        xx(3)=-4.058451513773971e-01;
        xx(4)=+0.000000000000000e+00;
        xx(5)=+4.058451513773971e-01;
        xx(6)=+7.415311855993950e-01;
        xx(7)=+9.491079123427590e-01;
        ww(1)=+1.294849661688624e-01;
        ww(2)=+2.797053914892770e-01;
        ww(3)=+3.818300505051190e-01;
        ww(4)=+4.179591836734694e-01;
        ww(5)=+3.818300505051190e-01;
        ww(6)=+2.797053914892770e-01;
        ww(7)=+1.294849661688624e-01;
    case 8
        xx(1)=-9.602898564975362e-01;
        xx(2)=-7.966664774136270e-01;
        xx(3)=-5.255324099163290e-01;
        xx(4)=-1.834346424956500e-01;
        xx(5)=+1.834346424956500e-01;
        xx(6)=+5.255324099163290e-01;
        xx(7)=+7.966664774136270e-01;
        xx(8)=+9.602898564975362e-01;
        ww(1)=+1.012285362903700e-01;
        ww(2)=+2.223810344533744e-01;
        ww(3)=+3.137066458778873e-01;
        ww(4)=+3.626837833783620e-01;
        ww(5)=+3.626837833783620e-01;
        ww(6)=+3.137066458778873e-01;
        ww(7)=+2.223810344533744e-01;
        ww(8)=+1.012285362903700e-01;
    case 9
        xx(1)=-9.681602395076261e-01;
        xx(2)=-8.360311073266360e-01;
        xx(3)=-6.133714327005904e-01;
        xx(4)=-3.242534234038090e-01;
        xx(5)=+0.000000000000000e+00;
        xx(6)=+3.242534234038090e-01;
        xx(7)=+6.133714327005904e-01;
        xx(8)=+8.360311073266360e-01;
        xx(9)=+9.681602395076261e-01;
        ww(1)=+8.127438836156873e-02;
        ww(2)=+1.806481606948580e-01;
        ww(3)=+2.606106964029360e-01;
        ww(4)=+3.123470770400020e-01;
        ww(5)=+3.302393550012600e-01;
        ww(6)=+3.123470770400020e-01;
        ww(7)=+2.606106964029360e-01;
        ww(8)=+1.806481606948580e-01;
        ww(9)=+8.127438836156873e-02;
    case 10
        xx(1)=-9.739065285171720e-01;
        xx(2)=-8.650633666889850e-01;
        xx(3)=-6.794095682990244e-01;
        xx(4)=-4.333953941292471e-01;
        xx(5)=-1.488743389816311e-01;
        xx(6)=+1.488743389816311e-01;
        xx(7)=+4.333953941292471e-01;
        xx(8)=+6.794095682990244e-01;
        xx(9)=+8.650633666889850e-01;
        xx(10)=+9.739065285171720e-01;
        ww(1)=+6.667134430868290e-02;
        ww(2)=+1.494513491505810e-01;
        ww(3)=+2.190863625159820e-01;
        ww(4)=+2.692667193099920e-01;
        ww(5)=+2.955242247147530e-01;
        ww(6)=+2.955242247147530e-01;
        ww(7)=+2.692667193099920e-01;
        ww(8)=+2.190863625159820e-01;
        ww(9)=+1.494513491505810e-01;
        ww(10)=+6.667134430868290e-02;
    otherwise
        error('max order: 10')
end
end
