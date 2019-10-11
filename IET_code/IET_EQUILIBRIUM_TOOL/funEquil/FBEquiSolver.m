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
function [EqSolution] = FBEquiSolver(app)
% imposizione correnti sui nodi del plasma
nn=size(app.node,1);
Iphi = zeros(nn,1);
Areatot = sum(app.Area_duale);
Jphi = app.Ipla/Areatot;
Iphi = Jphi * app.Area_duale;
lambda =app.Ipla/sum(Iphi);    %defining lambda to scale plasma current
Iphi=Iphi*lambda;           %total plasma current

%% Imposing boundary conditions BCs
Psi=zeros(nn,1);
disp('imposizione BCs e soluzione')
kk_bc = app.kk_nu;
nBCs = length(app.nodi_BCs);
Psi(app.nodi_BCs) = app.psiBoundary;

kk_bc(app.nodi_BCs,:)=0;
kk_bc(:,app.nodi_BCs)=0;
for ii=1:nBCs
    kk_bc(app.nodi_BCs(ii),app.nodi_BCs(ii))=1;
end
Iphi=Iphi-app.kk_nu(:,app.nodi_BCs)*Psi(app.nodi_BCs);
Iphi(app.nodi_BCs)=Psi(app.nodi_BCs);

% first solution
tic
sol_0=kk_bc\Iphi;
toc
solk=sol_0;

%% post processing: magnetic axis, contours of initial solution
rr = app.node(:,1);
zz = app.node(:,2);
npoint = 500;

rgrid=linspace(min(rr),max(rr),npoint);
zgrid=linspace(min(zz),max(zz),npoint);
[RR,ZZ] = meshgrid(rgrid,zgrid);
PSI = griddata(rr,zz,solk,RR,ZZ);

contour(app.PlasmaSurfacesPlot,RR,ZZ,PSI,100)
colormap(app.PlasmaSurfacesPlot,'jet')
colorbar(app.PlasmaSurfacesPlot,'vert')
plot(app.PlasmaSurfacesPlot,app.R_boundary_finale,app.Z_boundary_finale,'-r','linewidth',2);
title(app.PlasmaSurfacesPlot,'1st solution plasma ring')
% axis(app.PlasmaSurfacesPlot,[min(rr) max(rr) min(zz) max(zz)])

JP.faces=app.tri;
JP.vertices=app.node;
JP.facevertexcdata=Iphi;
patch(app.JprofileSolutionPlot,JP,'facecolor','interp','edgecolor','none');
colormap(app.JprofileSolutionPlot,'jet')
 colorbar(app.JprofileSolutionPlot,'vert')
axis(app.JprofileSolutionPlot,'equal');
grid(app.JprofileSolutionPlot,'on');
title(app.JprofileSolutionPlot,'1st Jpla uniform')
xlabel(app.ConvergencePlot,'R[m]')
ylabel(app.ConvergencePlot,'Z[m]')
% axis(app.PlasmaSurfacesPlot,[min(rr) max(rr) min(zz) max(zz)])

[psi_ak] = max(solk); %finding magnetic axis
errorAxis(1,1)=abs((100-psi_ak)/100)*100; %per cent relative magnetic axis error

%% iterative procedure solution for the magnetic poloidal flux
iter=0;
err = max(abs(solk)/(max(abs(solk))))*100; %Initial Picard error
tol = 0.01;

while err > tol
    
    psin = (solk-psi_ak)/(app.psiBoundary-psi_ak); %normalized poloidal magnetic flux
    
    
    switch app.JpType
        case 'Shafranov'
                fdfn= spline(app.psibar,app.FdF,psin);%interp1(psibar,FdF,psin,'spline');   %cubic interpolation
                dpn= spline(app.psibar,app.dP,psin);%interp1(psibar,dP,psin,'spline');     %cubic interpolation
                Jphi=(app.node(:,1).*dpn+fdfn./(app.mu0*app.node(:,1))); %computation of nodes plasma current density

        case 'Blum'            
            g_psi=abs((1-psin.^app.AlfaM).^app.AlfaN);
            Jphi=lambda*(app.node(:,1)*app.Beta0/app.R0+app.R0*(1-app.Beta0)./app.node(:,1)).*g_psi;
          
        case 'Arbitrary'
         Jphi=interp1(app.psibar,app.Jpoints,psin,'pchip');   %cubic interpolation
    end
            
                        
            Iphi=zeros(nn,1);
            Iphi=app.Area_duale.*Jphi;
            lambda = app.Ipla/sum(Iphi);                    %calculation of lambda and total current
            Iphi = Iphi*lambda;
            Iphi_tot = sum(Iphi);
            Iphi_plasma=Iphi;
  
            Psi(app.nodi_BCs) = app.psiBoundary;
            Iphi=Iphi-app.kk_nu(:,app.nodi_BCs)*Psi(app.nodi_BCs);
            Iphi(app.nodi_BCs)=Psi(app.nodi_BCs);
            
            
            solk1=kk_bc\Iphi;                           %solution
            
            
            [psi_aknew,ind_psiak] = max(solk1);            %finding magnetic axis
            errorAxis(iter+1,1) = abs((psi_aknew-psi_ak)/psi_ak)*100;  %calculation of error at axis
            
            psi_ak = psi_aknew;

            %calculation of iterative error (Picard)
            err = norm(solk1-solk)/(norm(abs(solk1)))*100;
            
            
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
            
end
    
            EqSolution.iter = iter;
            EqSolution.psi_tot_direct = solk;
            EqSolution.errorAxis = errorAxis;
            EqSolution.Raxis = rr(ind_psiak);
            EqSolution.Zaxis = zz(ind_psiak);
            EqSolution.Psi_axis = psi_ak;
            EqSolution.Iphi_plasma = Iphi_plasma;
    
            app.ConvergencePlot.cla;
plot(app.ConvergencePlot,1:iter,Residual_Picard,'r-x','linewidth',2)
            hold(app.ConvergencePlot,'on')
            grid(app.ConvergencePlot,'on')
xlabel(app.ConvergencePlot,'Iterations')
axis(app.ConvergencePlot,[0 iter 0 max(Residual_Picard)])
ylabel(app.ConvergencePlot,'\epsilon_{ \psi_a }')

app.PlasmaSurfacesPlot.cla;
PSI_final = griddata(app.node(:,1),app.node(:,2),EqSolution.psi_tot_direct ,RR,ZZ);
contour(app.PlasmaSurfacesPlot,RR,ZZ,PSI_final,30)
contour(app.PlasmaSurfacesPlot,RR,ZZ,PSI_final,[app.psiBoundary,app.psiBoundary],'r','linewidth',2) 
plot(app.PlasmaSurfacesPlot,app.R_boundary_finale,app.Z_boundary_finale,'-r','linewidth',2);
% axis(app.PlasmaSurfacesPlot,[min(rr) max(rr) min(zz) max(zz)])
title(app.PlasmaSurfacesPlot,'\Psi final equilibrium')

app.JprofileSolutionPlot.cla;
JP.facevertexcdata=Iphi_plasma./app.Area_duale;
patch(app.JprofileSolutionPlot,JP,'facecolor','interp','edgecolor','none');
hold(app.JprofileSolutionPlot,'on')
plot(app.JprofileSolutionPlot,app.R_boundary_finale,app.Z_boundary_finale,'-r','linewidth',2);
            grid(app.JprofileSolutionPlot,'on')
            axis(app.JprofileSolutionPlot,'equal')
title(app.JprofileSolutionPlot,'Jp final')

            EqSolution.RR = RR;
            EqSolution.ZZ = ZZ;
            EqSolution.PSI_final = PSI_final;
            EqSolution.JP = JP;
            EqSolution.psin = psin;
             EqSolution.Jphi = Jphi;

 switch app.JpType
        case 'Shafranov'
                EqSolution.fdfn= fdfn;
                EqSolution.dpn= dpn;
        case 'Blum'            
                EqSolution.g_psi= g_psi;
        case 'Arbitrary'
    end
   
end

