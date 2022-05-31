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
function [res]=fun_ThinSolenoid(r1,z1,z2,points)

%% parameters
eps=1.0d-15;
mu0=4*pi*1.e-7;
%%
ier=0;

PP.RR=points(:,1);
PP.ZZ=points(:,2);
aphi=zeros(size(PP.RR));
br=aphi;
bz=aphi;

IK=aphi;
IE=aphi;
IP=aphi;
IPmIK=aphi;
IKmIE=aphi;

J=1/(z2-z1); %linear current density [A/m] for unit current I [A]
rr=PP.RR;

rho1=rr/r1; % normalized r
sqrt_rho1=sqrt(rho1);
sqn=4.d0*rho1./((rho1+1).*(rho1+1));

zeta(:,1)=(z1-PP.ZZ)/r1; % normalized deltaz1
zeta(:,2)=(z2-PP.ZZ)/r1; % normalized deltaz2

%% find points on the axis
ind_axis=find(rho1<eps);
aphi(ind_axis)=0.d0;
br(ind_axis)=0.d0;
bz(ind_axis)=mu0*J/2.d0*(zeta(ind_axis,2)./sqrt(1.d0+zeta(ind_axis,2).*zeta(ind_axis,2))-zeta(ind_axis,1)./...
    sqrt(1.d0+zeta(ind_axis,1).*zeta(ind_axis,1)));
ind_all=1:numel(rho1);
ind_not=setdiff(ind_all,ind_axis);

%% All the points off the axis
ind_off=setdiff(ind_all(ind_not),ind_axis)';
ind_off_1=find(sqn(ind_not)>1.d0-eps);
ind_off_2=setdiff(ind_off,ind_off_1)';

% If the calculation points belongs to the thin solenoid    
if ~isempty(ind_off_1)
    for ii=1:2
        fk=(rho1+1).*(rho1+1)+zeta(:,ii).*zeta(:,ii);
        sqk=4d0*rho1./fk;
        kk=sqrt(sqk);
        for hh=1:numel(ind_off_1)
            hh_ok=ind_off_1(hh);
            [ik,ie,ikmie,~]=ellip_ke(sqk(hh_ok));
            IK(hh_ok,1)=ik;
            IE(hh_ok,1)=ie;
            IKmIE(hh_ok,1)=ikmie;
        end       
        aphi(ind_off_1)=aphi(ind_off_1)+(-1.d0)^ii.*zeta(ind_off_1,ii)./...
            kk(ind_off_1).*IKmIE(ind_off_1);
        br(ind_off_1)=br(ind_off_1)+(-1.d0)^ii.*1.d0./kk(ind_off_1).*...
            ((1-0.5d0.*sqk(ind_off_1)).*IK(ind_off_1)-IE(ind_off_1));
        bz(ind_off_1)=bz(ind_off_1)+(-1.d0)^ii.*zeta(ind_off_1,ii).*...
            kk(ind_off_1).*IK(ind_off_1);
    end
    aphi(ind_off_1)=aphi(ind_off_1)*2.d-7*J*r1;
    br(ind_off_1)=br(ind_off_1)*4.d-7*J;
    bz(ind_off_1)=bz(ind_off_1)*1.d-7*J./sqrt_rho1(ind_off_1);
end

% If the calculation points does not belong to the thin solenoid     
for ii=1:2
    fk=(rho1+1).*(rho1+1)+zeta(:,ii).*zeta(:,ii);
    sqk=4d0*rho1./fk;
    kk=sqrt(sqk);
    for hh=1:numel(ind_off_2)
        hh_ok=ind_off_2(hh);
        [ik,ie,ip,ipmik,ipmie,~]=ellip_kep(sqn(hh_ok),sqk(hh_ok));
        IK(hh_ok,1)=ik;
        IE(hh_ok,1)=ie;
        IP(hh_ok,1)=ip;
        IPmIK(hh_ok,1)=ipmik;
        IKmIE(hh_ok,1)=ipmie;
    end
    aphi(ind_off_2)=aphi(ind_off_2)+(-1.d0)^ii*zeta(ind_off_2,ii)./kk(ind_off_2).*(sqk(ind_off_2).*...
        (sqn(ind_off_2)-1.d0)./sqn(ind_off_2).*IPmIK(ind_off_2)+IKmIE(ind_off_2));
    br(ind_off_2)=br(ind_off_2)+(-1.d0)^ii./kk(ind_off_2).*((1-0.5d0.*kk(ind_off_2).*kk(ind_off_2)).*...
        IK(ind_off_2)-IE(ind_off_2));
    bz(ind_off_2)=bz(ind_off_2)+(-1.d0)^ii.*zeta(ind_off_2,ii).*kk(ind_off_2).*(-(rho1(ind_off_2)-1.d0)./...
        (rho1(ind_off_2)+1.d0).*IP(ind_off_2)+IK(ind_off_2));
end
aphi(ind_off_2)=aphi(ind_off_2).*2.d-7.*J.*r1./sqrt_rho1(ind_off_2);
br(ind_off_2)=br(ind_off_2).*4.d-7*J./sqrt_rho1(ind_off_2);
bz(ind_off_2)=bz(ind_off_2).*1.d-7*J./sqrt_rho1(ind_off_2);


%% Output
res.a=aphi;
res.br=br;
res.bz=bz;
res.psi=aphi*2*pi.*rr;
end 
