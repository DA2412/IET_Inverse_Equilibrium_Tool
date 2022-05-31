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
function [IK,IE,IP,IPmIK,IKmIE,ier]=ellip_kep(sqn_hh,sqk_hh)

% complete elliptic integrals of first, second, third type: K(sqk), E(sqk), P(sqn,sqk)
%
% accuracy: eps
%
% ier:      errore index:     =0: OK
%                             =1: singularity for sqk
%                             =2: singularity for sqn
%                             =3: max iteration exceeded
%
% Iterative procedure, based on Landen's transformation, described in 
% M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions", 
% Dover Publications", 1965, Chapter 17.6, pag. 598-599

%% parameters
      imax=100;
      eps=1.0d-15;
      pihalf=1.570796326794897d0;
%% main code
      
      ier=0;
      
      if(sqk_hh>1.d0-eps) 
        disp('*** FATAL ERROR in ellip_kep: sqk>1-eps ***');
        ier=1;
      end

      
      if(sqn_hh>1.d0-eps) 
        disp('*** FATAL ERROR in ellip_kep: sqn>1-eps ***');
        ier=2;
      end

      if(sqk_hh <= eps) 
	    IK=pihalf;
	    IE=pihalf;
	    IP=pihalf;	  
      end

      aa0=1.d0;
      bb0=sqrt(1.d0-sqk_hh);
      cc0=sqrt(sqk_hh);
      dd0=(1.d0-sqn_hh)/bb0;
      ee0=sqn_hh/(1.d0-sqn_hh);
      ff0=0;
      
      sumc=cc0*cc0;
      ii=0;
      
     while ((cc0>eps) || (dd0-1.0d0>eps))
        ii=ii+1;
      
        aa=0.5d0*(aa0+bb0);
        bb=sqrt(aa0*bb0);
        cc=0.5d0*(aa0-bb0);
        
        dd=bb/(4.d0*aa)*(2.d0+dd0+1.d0/dd0);
        ee=(dd0*ee0+ff0)/(1.d0+dd0);
        ff=0.5d0*(ee0+ff0);
      
        aa0=aa;
        bb0=bb;
        cc0=cc;
        dd0=dd;
        ee0=ee;
        ff0=ff;
      
        sumc=sumc+2^ii*cc*cc;
         
        if(ii>imax) 
            disp(['*** ii>imax=',num2str(imax),' in ellipkep ***']);
            disp(['*** Too many iterations ***']);
            ier=3;
            return
        end
      
     end
      
% 17.6.3 kkk=pihalf/aa   
      IK=pihalf/(aa);
% 17.6.4 eee=kkk*(1-sumc/2)
      IE=IK*(1.d0-sumc*0.5d0);
      IKmIE=IK*sumc*0.5d0;
%  ppp=kkk*(1+ff)
      IP=IK*(1+ff);
%  pmk=kkk*ff)     
      IPmIK=IK*ff;
end