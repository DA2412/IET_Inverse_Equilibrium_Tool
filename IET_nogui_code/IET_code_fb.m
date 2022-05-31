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

% mesh data
tri = OutputShape.t;
node = OutputShape.p;
ntri=size(tri,1);
nfaces=size(tri,2);
nodi_BCs=reshape(OutputShape.b,numel(OutputShape.b),1);
nn=size(node,1);
Iphi = zeros(nn,1);
Area_duale = OutputShape.Area_duale;
Areatot = sum(Area_duale);
%% plasma eq data

R_boundary_finale = OutputShape.R_boundary_finale;
Z_boundary_finale = OutputShape.Z_boundary_finale;

Jphi = Ipla/Areatot;
Iphi = Jphi * Area_duale;
lambda =Ipla/sum(Iphi);    %defining lambda to scale plasma current
Iphi=Iphi*lambda;           %total plasma current

%% Imposing boundary conditions BCs
Psi=zeros(nn,1);
disp('imposizione BCs e soluzione')
kk_nu = OutputShape.kk_nu;
kk_bc = OutputShape.kk_nu;
nodi_BCs = OutputShape.b;
nBCs = length(nodi_BCs);
Psi(nodi_BCs) = psiBoundary;

kk_bc(nodi_BCs,:)=0;
kk_bc(:,nodi_BCs)=0;
for ii=1:nBCs
    kk_bc(nodi_BCs(ii),nodi_BCs(ii))=1;
end
Iphi=Iphi-kk_nu(:,nodi_BCs)*Psi(nodi_BCs);
Iphi(nodi_BCs)=Psi(nodi_BCs);

% first solution
tic
sol_0=kk_bc\Iphi;
toc
solk=sol_0;

[psiAxis0] = max(solk); %finding magnetic axis

%% post processing: magnetic axis, contours of initial solution
rr = node(:,1);
zz = node(:,2);
npoint = ngrid;

rgrid=linspace(min(rr),max(rr),npoint);
zgrid=linspace(min(zz),max(zz),npoint);
[RR,ZZ] = meshgrid(rgrid,zgrid);
PSI = griddata(rr,zz,solk,RR,ZZ);
[psi_ak] = max(solk); %finding magnetic axis
errorAxis(1,1)=(psi_ak - psiAxis0)/abs(psiAxis0)*100; %per cent relative magnetic axis error

%% iterative procedure solution for the magnetic poloidal flux
iter=0;
err = 1;%
scarto = 1;%

while err > tol
    
    psin = (solk-psi_ak)/(psiBoundary-psi_ak); %normalized poloidal magnetic flux
    
    
    switch JpType
        case 'Shafranov'
            psibar = psibar;
            FdF = FdF;
            dP = dP;
            fdfn= interp1(psibar,FdF,psin,'pchip');%interp1(psibar,FdF,psin,'spline');   %cubic interpolation
            dpn= interp1(psibar,dP,psin,'pchip');%interp1(psibar,dP,psin,'spline');     %cubic interpolation
            Jphi=(rr(:,1).*dpn+fdfn./(mu0*rr(:,1))); %computation of nodes plasma current density

        case 'Blum'     
            alpha_M = AlphaM;
           alpha_N = AlphaN;
            beta_0 = Beta0;
            R_0 = R0;
            g_psi=((1-psin.^alpha_M).^alpha_N);
            h_f = (rr(:,1)*beta_0/R_0+R_0*(1-beta_0)./rr(:,1));
             lambda=1;
            Jphi=lambda*h_f.*g_psi;
          
        case 'Arbitrary'
         Jphi=interp1(psibar,Jpoints,psin,'pchip');   %cubic interpolation
                
    end
            
                        
            Iphi=zeros(nn,1);
            Iphi=Area_duale.*Jphi;
            lambda = Ipla/sum(Iphi);                    %calculation of lambda and total current
            Iphi = Iphi*lambda;
            Iphi_tot = sum(Iphi);
            Iphi_plasma=Iphi;
  
            Psi(nodi_BCs) = psiBoundary;
            Iphi=Iphi-kk_nu(:,nodi_BCs)*Psi(nodi_BCs);
            Iphi(nodi_BCs)=Psi(nodi_BCs);
            
            
            solk1=kk_bc\Iphi;                           %solution
            
            
            [psi_aknew,ind_psiak] = max(solk1);            %finding magnetic axis
            errorAxis(iter+1,1) = abs((psi_aknew-psi_ak)/psi_ak)*100;  %calculation of error at axis
            
            psi_ak = psi_aknew;

            %calculation of iterative error (Picard)
            err = norm(solk1-solk);
            scarto(iter+1,:) = norm(solk1-solk)/numel(solk1);
         
            
            disp(['iteration n°:',num2str(iter+1)]);
            disp(['lambda: ',num2str(lambda)])
            disp(['I tot: ',num2str(Iphi_tot)])
            disp(['psi_ak magnetic axis: ',num2str(psi_ak)]);
            disp(['ra magnetic axis: ',num2str(rr(ind_psiak))]);
            disp(['za magnetic axis: ',num2str(zz(ind_psiak))]);
            disp(['magnetic flux error at axis: ',num2str(errorAxis(iter+1))])
            disp(['Picard error: ',num2str(err)])
            disp(['---------------------------'])
            
             solk = solk1;                       % re-assignement of iterative solution
            Residual_Picard(iter+1,:)=err;
            iter = iter+1;    
            
            if iter>50
                disp(['------EXCEEDED 50 ITERATIONS--------------------'])
               break
            end
            
       
end
    
EqSolution.Jphi_final = (Iphi_plasma./Area_duale);
EqSolution.iter = iter;
EqSolution.psi_tot_direct = solk;
EqSolution.errorAxis = errorAxis;
EqSolution.Raxis = rr(ind_psiak);
EqSolution.Zaxis = zz(ind_psiak);
EqSolution.Psi_axis = psi_ak;
EqSolution.psiBoundary = psiBoundary;
EqSolution.Iphi_plasma = Iphi_plasma;
EqSolution.major_radius = (max(R_boundary_finale) + min(R_boundary_finale))/2;
EqSolution.minor_radius = (max(R_boundary_finale) - min(R_boundary_finale))/2;
EqSolution.RR = RR;
EqSolution.ZZ = ZZ;
EqSolution.PSI_final = griddata(node(:,1),node(:,2),EqSolution.psi_tot_direct ,RR,ZZ,'cubic');

            
 switch JpType
     
    case 'Shafranov'
        EqSolution.psin = psin;
        EqSolution.fdfn = fdfn;
        EqSolution.dpn = dpn;
        
    case 'Blum'
        EqSolution.alpha_M = alpha_M;
        EqSolution.alpha_N = alpha_N;
        EqSolution.beta_0 = beta_0;
        EqSolution.R_0 = R_0;
        EqSolution.Jphi = Jphi;
        EqSolution.psin = psin;
        EqSolution.h_f = h_f;
        EqSolution.g_psi = g_psi;
        EqSolution.lambda = lambda;
        
    case 'Arbitrary'
        Jphi=interp1(psibar,Jpoints,psin,'pchip');   %cubic interpolation
        
    case 'RFP'
        EqSolution.psin = psin;
        EqSolution.fn = fn;
        EqSolution.dfn = dfn;
        
 end   

OutputEquil.Equil = EqSolution;
OutputEquil.OutputShape = OutputShape;