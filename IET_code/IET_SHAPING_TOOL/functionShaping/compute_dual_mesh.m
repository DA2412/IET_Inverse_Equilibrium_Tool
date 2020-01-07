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
function [Area_duale, kk_nu] = compute_dual_mesh (app,node, tri, mu0)
nn=size(node,1);
ntri=size(tri,1);
mur=1;
mur=ones(1,ntri)*mur;

%% costruzione delle matrice con assemblaggio locale
kk_nu=zeros(nn,nn);

Area=zeros(nn,1);
Psi=zeros(nn,1);
Iphi=zeros(nn,1);

Ctilde = [ 0 -1 -1; 
          -1  0  1;
           1  1  0];
        
C = [ 0 -1 1; 
     -1  0 1;
     -1  1 0];

disegnaMesh=1;

%% Assemblaggio matrice globale
tic
%disp('computing global matri')
for ii=1:ntri
    % riordino i nodi, dal pi? basso al pi? alto
    nodi_loc = sort(tri(ii,:));
    
    % nodi prmali
    n1 = node(nodi_loc(1),1:2);
    n2 = node(nodi_loc(2),1:2);
    n3 = node(nodi_loc(3),1:2);  
    
    % baricentro
    nbar = (n1+n2+n3)/3;
    nb(ii,:)=nbar;
    
    % lati primali
    e1 = n3 - n2;
    e2 = n3 - n1;
    e3 = n2 - n1;
    
    % nodi duali
    n1d = 0.5*(n2+n3);
    n2d = 0.5*(n1+n3);
    n3d = 0.5*(n1+n2);

    % lati duali
    e1d = nbar - n1d;
    e2d = n2d- nbar;
    e3d = nbar - n3d;
      
% Area cella primale
Area2=abs(e2(1)*e3(2)-e3(1)*e2(2));
indn_glob(1,1:3)=tri(ii,1:3);

%% calcolo delle matrici Me ed Mf:
M_mu=[1 0;0 1]/(mur(ii)*mu0);
Me=[e1d;e2d;e3d];

if (n3(1)~=0)
    Mf=[-e2(1)/n1d(1), +e1(1)/n2d(1), 0;
        -e2(2)/n1d(1), +e1(2)/n2d(1), 0]/(Area2*2*pi);
else  if (n1(1)~=0)
        Mf=[0, -e3(1)/n2d(1), +e2(1)/n3d(1);
            0, -e3(2)/n2d(1), +e2(2)/n3d(1)]/(Area2*2*pi);
       else
        Mf=[-e3(1)/n1d(1), 0, +e1(1)/n3d(1);
            -e3(2)/n1d(1), 0, +e1(2)/n3d(2)]/(Area2*2*pi);
       end
end

%% calcolo della matrice locale: kk_loc=Ctilde*(Me*M_mu*Mf)*C-i*omega*Msigma
 kk_loc=Ctilde*(Me*M_mu*Mf)*C;

%% assegnazione valori alla matrice globale 
 kk_nu(nodi_loc(1), nodi_loc(1)) = kk_nu(nodi_loc(1), nodi_loc(1)) + kk_loc(1,1);
 kk_nu(nodi_loc(1), nodi_loc(2)) = kk_nu(nodi_loc(1), nodi_loc(2)) + kk_loc(1,2);
 kk_nu(nodi_loc(1), nodi_loc(3)) = kk_nu(nodi_loc(1), nodi_loc(3)) + kk_loc(1,3);

 kk_nu(nodi_loc(2), nodi_loc(1)) = kk_nu(nodi_loc(2), nodi_loc(1)) + kk_loc(2,1);
 kk_nu(nodi_loc(2), nodi_loc(2)) = kk_nu(nodi_loc(2), nodi_loc(2)) + kk_loc(2,2);
 kk_nu(nodi_loc(2), nodi_loc(3)) = kk_nu(nodi_loc(2), nodi_loc(3)) + kk_loc(2,3);

 kk_nu(nodi_loc(3), nodi_loc(1)) = kk_nu(nodi_loc(3), nodi_loc(1)) + kk_loc(3,1);
 kk_nu(nodi_loc(3), nodi_loc(2)) = kk_nu(nodi_loc(3), nodi_loc(2)) + kk_loc(3,2);
 kk_nu(nodi_loc(3), nodi_loc(3)) = kk_nu(nodi_loc(3), nodi_loc(3)) + kk_loc(3,3);

%% reticolo duale 
% Area duale 1 (affacciata al nodo 1)
Atmp1 = abs(det([e3d; e3/2]));
Atmp2 = abs(det([e2d; e2/2]));
A1d=0.5*(Atmp1+Atmp2);
%nb_A1d(nodi_loc(1)) = (tmp1 * Atmp1 + tmp2 * Atmp2)/A1d*2;
    
% Area duale 2 (affacciata al nodo 2)    
Atmp1 = abs(det([e3d; e3/2]));
Atmp2 = abs(det([e1d; e1/2]));
A2d=0.5*(Atmp1+Atmp2);
%nb_A2d(nodi_loc(2)) = (tmp1 * Atmp1 + tmp2 * Atmp2)/A2d*2;
    
% Area duale 3 (affacciata al nodo 3)   
Atmp1 = abs(det([e1d; e1/2]));
Atmp2 = abs(det([e2d; e2/2]));
A3d=0.5*(Atmp1+Atmp2);
%nb_A3d(nodi_loc(3)) = (tmp1 * Atmp1 + tmp2 * Atmp2)/A3d*2;
            
Area(nodi_loc(1))=Area(nodi_loc(1))+A1d;
Area(nodi_loc(2))=Area(nodi_loc(2))+A2d;
Area(nodi_loc(3))=Area(nodi_loc(3))+A3d;
 
%if disegnaduale 
%      plot([nbar(1) n1d(1)],[nbar(2) n1d(2)],'-r')
%      plot([nbar(1) n2d(1)],[nbar(2) n2d(2)],'-r')
%      plot([nbar(1) n3d(1)],[nbar(2) n3d(2)],'-r')
%end  
% 
% % end
% %% reticolo duale 
% % Area duale 1 (affacciata al nodo 1)
% tmp1 = (n1 + n3d + nbar)/3;
% tmp2 = (n1 + n2d + nbar)/3;
% Atmp1 = abs(det([e3d; e3/2]));
% Atmp2 = abs(det([e2d; e2/2]));
% A1d=0.5*(Atmp1+Atmp2);
% %nb_A1d(nodi_loc(1)) = (tmp1 * Atmp1 + tmp2 * Atmp2)/A1d*2;
%     
% % Area duale 2 (affacciata al nodo 2)    
% tmp1 = (n2 + n3d + nbar)/3;
% tmp2 = (n2 + n1d + nbar)/3;
% Atmp1 = abs(det([e3d; e3/2]));
% Atmp2 = abs(det([e1d; e1/2]));
% A2d=0.5*(Atmp1+Atmp2);
% %nb_A2d(nodi_loc(2)) = (tmp1 * Atmp1 + tmp2 * Atmp2)/A2d*2;
%     
% % Area duale 3 (affacciata al nodo 3)   
% tmp1 = (n3 + n1d + nbar)/3;
% tmp2 = (n3 + n2d + nbar)/3;
% Atmp1 = abs(det([e1d; e1/2]));
% Atmp2 = abs(det([e2d; e2/2]));
% A3d=0.5*(Atmp1+Atmp2);
% %nb_A3d(nodi_loc(3)) = (tmp1 * Atmp1 + tmp2 * Atmp2)/A3d*2;
%             
% Area(nodi_loc(1))=Area(nodi_loc(1))+A1d;
% Area(nodi_loc(2))=Area(nodi_loc(2))+A2d;
% Area(nodi_loc(3))=Area(nodi_loc(3))+A3d;
%  
% if disegnaduale 
%      plot([nbar(1) n1d(1)],[nbar(2) n1d(2)],'-r')
%      plot([nbar(1) n2d(1)],[nbar(2) n2d(2)],'-r')
%      plot([nbar(1) n3d(1)],[nbar(2) n3d(2)],'-r')
% end  
if disegnaMesh

            plot(app.PlasmaShape,[nbar(1) n1d(1)],[nbar(2) n1d(2)],'-r')
     plot(app.PlasmaShape,[nbar(1) n2d(1)],[nbar(2) n2d(2)],'-r')
     plot(app.PlasmaShape,[nbar(1) n3d(1)],[nbar(2) n3d(2)],'-r')

end

Area_duale = Area;
kk_nu=sparse(kk_nu);
%save cell Area_duale kk_nu
end

