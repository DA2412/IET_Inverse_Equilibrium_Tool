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
Area_duale = app.Area_duale;
Areatot = sum(Area_duale);
psiBoundary = app.psiBoundary;
R_boundary_finale = app.R_boundary_finale;
Z_boundary_finale = app.Z_boundary_finale;

JpType = app.JpType;

Jphi = app.Ipla/Areatot;
Iphi = Jphi * Area_duale;
lambda =app.Ipla/sum(Iphi);    %defining lambda to scale plasma current
Iphi=Iphi*lambda;           %total plasma current
tol = app.tol;

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

[psiAxis0] = max(solk); %finding magnetic axis

%% post processing: magnetic axis, contours of initial solution
rr = app.node(:,1);
zz = app.node(:,2);
npoint = app.ngrid;

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
errorAxis(1,1)=(psi_ak - psiAxis0)/abs(psiAxis0)*100; %per cent relative magnetic axis error

%% iterative procedure solution for the magnetic poloidal flux
iter=0;
err = 1;%
scarto = 1;%

while err > tol
    
    psin = (solk-psi_ak)/(app.psiBoundary-psi_ak); %normalized poloidal magnetic flux
    
    
    switch app.JpType
        case 'Shafranov'
            fdfn= interp1(psibar,FdF,psin,'pchip');%interp1(psibar,FdF,psin,'spline');   %cubic interpolation
            dpn= interp1(psibar,dP,psin,'pchip');%interp1(psibar,dP,psin,'spline');     %cubic interpolation
            Jphi=(rr(:,1).*dpn+fdfn./(mu0*rr(:,1))); %computation of nodes plasma current density

        case 'Blum'     
            alpha_M = app.AlphaM;
           alpha_N = app.AlphaN;
            beta_0 = app.Beta0;
            R_0 = app.R0;
            g_psi=((1-psin.^alpha_M).^alpha_N);
            h_f = (rr(:,1)*beta_0/R_0+R_0*(1-beta_0)./rr(:,1));
             lambda=1;
            Jphi=lambda*h_f.*g_psi;
          
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
EqSolution.Psi_boundary = psiBoundary;
EqSolution.Iphi_plasma = Iphi_plasma;
EqSolution.major_radius = (max(R_boundary_finale) + min(R_boundary_finale))/2;
EqSolution.minor_radius = (max(R_boundary_finale) - min(R_boundary_finale))/2;
EqSolution.RR = RR;
EqSolution.ZZ = ZZ;
EqSolution.PSI_final = griddata(app.node(:,1),app.node(:,2),EqSolution.psi_tot_direct ,RR,ZZ,'cubic');

            
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
        
%     case 'RFP'
%         EqSolution.psin = psin;
%         EqSolution.fn = fn;
%         EqSolution.dfn = dfn;
        
 end   

 
 app.ConvergencePlot.cla;
 semilogy(app.ConvergencePlot,Residual_Picard,'b-x','linewidth',2)
 hold(app.ConvergencePlot,'on')
  semilogy(app.ConvergencePlot,scarto,'r-o','linewidth',2)
 grid(app.ConvergencePlot,'on')
 xlabel(app.ConvergencePlot,'Iterations')
 axis(app.ConvergencePlot,[0 iter 0 max(Residual_Picard)])
 ylabel(app.ConvergencePlot,'Picard convergence')
 legend(app.ConvergencePlot,'error','residual')

app.PlasmaSurfacesPlot.cla;
contour(app.PlasmaSurfacesPlot,RR,ZZ,EqSolution.PSI_final,30)
contour(app.PlasmaSurfacesPlot,RR,ZZ,EqSolution.PSI_final,[app.psiBoundary,app.psiBoundary],'r','linewidth',2) 
plot(app.PlasmaSurfacesPlot,app.R_boundary_finale,app.Z_boundary_finale,'-r','linewidth',2);
title(app.PlasmaSurfacesPlot,'\Psi final equilibrium')


app.JprofileSolutionPlot.cla;
JP.faces=app.tri;
JP.vertices=app.node;
JP.facevertexcdata=EqSolution.Jphi_final./(max(EqSolution.Jphi_final));
patch(app.JprofileSolutionPlot,JP,'facecolor','interp','edgecolor','none');
hold(app.JprofileSolutionPlot,'on')
plot(app.JprofileSolutionPlot,app.R_boundary_finale,app.Z_boundary_finale,'-r','linewidth',2);
grid(app.JprofileSolutionPlot,'on')
axis(app.JprofileSolutionPlot,'equal')
title(app.JprofileSolutionPlot,'Jp final')


end

