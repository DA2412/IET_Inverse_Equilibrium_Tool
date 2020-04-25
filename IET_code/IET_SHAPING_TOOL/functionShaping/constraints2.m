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
function F=constraints2(a,delta_u,k_second,csi,s)
    c=csi(2); %csi upper inner 2quadrant
    F=zeros(3,1);
    F(1)=(1-a*(1-delta_u)/s(1))^(s(3))+(k_second*a/s(2))^s(3)-1;
    F(2)=(1+a*(1-delta_u)*(1-1/(2^0.5))*(c-1)/s(1))^s(3)+(k_second*a*(1/(2^0.5) + c*(1-1/(2^0.5)))/s(2) )^s(3)-1;
    if delta_u>=0
        F(3)=(s(2)* ((s(1)-a*(1-delta_u))/(k_second*a))^(s(3)-1) ) * (s(2)/s(1))^(s(3)-1)- (1-delta_u);
    else
        F(3)=(s(2)* ((s(1)-a*(1-delta_u))/(k_second*a))^(s(3)-1) ) * (s(2)/s(1))^(s(3)-1)- 1/(1+delta_u);
    end
end