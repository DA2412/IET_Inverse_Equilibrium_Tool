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
function [ res ] = fun_Field_Loop( source, point ) %#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Computation of Aphi, Br, Bz, Psi for axisymmetric loop current (no
%   thickness of the loop)
%
%   source (sources' geometry) - structure including:  
%     - R: radial distance form the axis of the coil's centre [m] 
%     - Z: vertical distance from z=0 plane of the coil's centre [m]
%
%   point (evaluation points) - structure including:  
%     - RR: array of the radial coordinate of the evaluatin points [m]
%     - ZZ: array of the vertical coordinate of the evaluatin points [m]
%
%   res (results) - structure including:
%    - a (npt x ncoil)
%    - br (npt x ncoil)
%    - bz (npt x ncoil)
%    - psi (npt x ncoil)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
npt=numel(point.RR);
ncoil=numel(source.R); 
a=zeros(npt,ncoil);
br=zeros(npt,ncoil);
bz=zeros(npt,ncoil);
psi=zeros(npt,ncoil);
RR=point.RR;
ZZ=point.ZZ;

%% Computation (vectorization of quantities)
for  jj=1:ncoil
    r0=source.R(jj);
    z0=source.Z(jj);
    kk=sqrt(4*r0*RR./((r0+RR).^2+(ZZ-z0).^2));
    [J1,J2] = ellipke(kk.^2);
    res_a=4.d-7./kk.*sqrt(r0./RR).*((1-kk.^2/2).*J1-J2);
    res_psi=4.d-7./kk.*sqrt(r0./RR).*((1-kk.^2/2).*J1-J2).*(2*pi*RR);
    res_br=+1.d-7.*kk.*(ZZ-z0)./(RR.*sqrt(r0.*RR)).*...
        (-J1+(r0^2+RR.^2+(ZZ-z0).^2)./((r0-RR).^2+(ZZ-z0).^2).*J2);
    res_bz=+1.d-7.*kk./sqrt(r0.*RR).*(J1+(r0^2-RR.^2-(ZZ-z0).^2)./...
        ((r0-RR).^2+(ZZ-z0).^2).*J2);
    a(:,jj)=res_a;
    psi(:,jj)=res_psi;
    br(:,jj)=res_br;
    bz(:,jj)=res_bz;
end

%% Find points on axis (r=0)
ind_axis=find(point.RR==0);
a(ind_axis,:)=0;
psi(ind_axis,:)=0;
br(ind_axis,:)=0;
if isempty(ind_axis)==0
    for  jj=1:ncoil
        r0=source.R(jj);
        z0=source.Z(jj);
        bz(ind_axis,jj)=4.d-7*pi*r0/(2*(r0^2+(ZZ(ind_axis)-z0).^2).^1.5d0);
    end
end

%% Output
res.a=a;
res.br=br;
res.bz=bz;
res.psi=psi;

end

