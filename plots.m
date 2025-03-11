%%

time=[0:time_saved_resolution:stoptime-time_saved_resolution];
depths=[0:0.002:0.4];
%%
F_Cai=(D_dCa*phi(1)*(dCaf(1,:)-dCaw)./dbl)';
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
F_H2Si=(D_dtH2S*phi(1)*(dtH2Sf(1,:)-dtH2Sw)./dbl)';

%in 1 m2, over 10 years
%Bud_acc=Foc_s*3/12; 
%Bud_PIC=(pcalcitef(:,10000)-pcalcitef(:,1))./phi';


CO2SYS_data = CO2SYS(dtalkf(1,:)*1e6/rho_sw,dtCO2f(1,:)*1e6/rho_sw,1,2,S,T,T,P,P,dtSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
rHCO3=CO2SYS_data(:,21)'./(dtCO2f(1,:)*1e6/rho_sw);
rCO3=CO2SYS_data(:,22)'./(dtCO2f(1,:)*1e6/rho_sw);
rCO2=CO2SYS_data(:,23)'./(dtCO2f(1,:)*1e6/rho_sw);
pCO2=CO2SYS_data(:,19);

Bud_acc(1,1)=0;
for a=2:10000
    Bud_DICflux(1,a)=sum((F_DICi(1:a,1)-F_DICi(1,1)).*time_saved_resolution);
    Bud_CO2flux(1,a)=sum((F_DICi(1:a,1)-F_DICi(1,1)).*time_saved_resolution.*rCO2(1,a)');
    Bud_HCO3flux(1,a)=sum((F_DICi(1:a,1)-F_DICi(1,1)).*time_saved_resolution.*rHCO3(1,a)');
    Bud_CO3flux(1,a)=sum((F_DICi(1:a,1)-F_DICi(1,1)).*time_saved_resolution.*rCO3(1,a)');
    Bud_TAflux(1,a)=sum((F_TAi(1:a,1)-F_TAi(1,1)).*time_saved_resolution);
    Bud_O2flux(1,a)=sum((F_O2i(1:a,1)-F_O2i(1,1)).*time_saved_resolution);
    Bud_Caflux(1,a)=sum((F_Cai(1:a,1)-F_Cai(1,1)).*time_saved_resolution);

    if a<250 %250   83   8 
        Bud_acc(1,a)=Foc_s.*time_saved_resolution+Bud_acc(1,a-1);
    else
        Bud_acc(1,a)=Bud_acc(1,a-1);
    end
end

CO2SYS_data = CO2SYS(dtalkf(1,1)*1e6/rho_sw,dtCO2f(1,1)*1e6/rho_sw,1,2,S,T,T,P,P,dtSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
rHCO3_ini=CO2SYS_data(1,21);
rCO3_ini=CO2SYS_data(1,22);
rCO2_ini=CO2SYS_data(1,23);
CO2SYS_data = CO2SYS(dtalkf(1,10000)*1e6/rho_sw,dtCO2f(1,10000)*1e6/rho_sw,1,2,S,T,T,P,P,dtSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
rHCO3_fin=CO2SYS_data(1,21);
rCO3_fin=CO2SYS_data(1,22);
rCO2_fin=CO2SYS_data(1,23);

a(1,1)=Bud_acc(1,10000);
a(2,1)=Bud_TAflux(1,10000);
a(3,1)=Bud_DICflux(1,10000);
a(4,1)=Bud_O2flux(1,10000);
a(5,1)=Bud_CO2flux(1,10000);
a(6,1)=Bud_HCO3flux(1,10000);
a(7,1)=Bud_CO3flux(1,10000);
a(8,1)=Bud_Caflux(1,10000);
a(9,1)=sum(procsf(:,10000).*z_res'.*phiS')+sum(psocsf(:,10000).*z_res'.*phiS')+sum(pfocsf(:,10000).*z_res'.*phiS');


%%
CO2SYS_data = CO2SYS(dtalkf(1,:)*1e6/rho_sw,dtCO2f(1,:)*1e6/rho_sw,1,2,S,T,T,P,P,dtSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
rHCO3=CO2SYS_data(:,21);
rCO3=CO2SYS_data(:,22);
rCO2=CO2SYS_data(:,23);
rtot=rHCO3+rCO3+rCO2;
clf
plot(time,rHCO3./rtot)
hold on
plot(time,rCO3./rtot)
hold on
plot(time,rCO2./rtot)

%% Fate of POC
% not sure about x axis: POC acc flux during deployment? global acc flux
% after 10 years?
%faire ça pour 3 Fpoc et 2 sed types
load('results.mat')
figure (1)
clf

subplot 321
plot(a(3,:),a(4,:).*a(16,:),'color',[82/255, 55/255, 156/255],'linewidth',2)
hold on
plot(a(3,:),a(4,:).*a(16,:)+a(4,:).*a(15,:),'color',[26/255, 174/255, 185/255],'linewidth',2)
plot(a(3,:),a(4,:),'color',[226/255, 97/255, 0/255],'linewidth',2)
xlim([0 15])
ylim([0 4])

subplot 322
plot(d(3,:),d(4,:).*d(16,:),'color',[82/255, 55/255, 156/255],'linewidth',2)
hold on
plot(d(3,:),d(4,:).*d(16,:)+d(4,:).*d(15,:),'color',[26/255, 174/255, 185/255],'linewidth',2)
plot(d(3,:),d(4,:),'color',[226/255, 97/255, 0/255],'linewidth',2)
xlim([0 15])
ylim([0 4])

subplot 323
plot(b(3,:),b(4,:).*b(16,:),'color',[82/255, 55/255, 156/255],'linewidth',2)
hold on
plot(b(3,:),b(4,:).*b(16,:)+b(4,:).*b(15,:),'color',[26/255, 174/255, 185/255],'linewidth',2)
plot(b(3,:),b(4,:),'color',[226/255, 97/255, 0/255],'linewidth',2)
xlim([0 15])
ylim([0 1.5])

subplot 324
plot(e(3,:),e(4,:).*e(16,:),'color',[82/255, 55/255, 156/255],'linewidth',2)
hold on
plot(e(3,:),e(4,:).*e(16,:)+e(4,:).*e(15,:),'color',[26/255, 174/255, 185/255],'linewidth',2)
plot(e(3,:),e(4,:),'color',[226/255, 97/255, 0/255],'linewidth',2)
xlim([0 15])
ylim([0 1.5])

subplot 325
plot(c(3,:),c(4,:).*c(16,:),'color',[82/255, 55/255, 156/255],'linewidth',2)
hold on
plot(c(3,:),c(4,:).*c(16,:)+c(4,:).*c(15,:),'color',[26/255, 174/255, 185/255],'linewidth',2)
plot(c(3,:),c(4,:),'color',[226/255, 97/255, 0/255],'linewidth',2)
xlim([0 15])
ylim([0 0.3])

subplot 326
plot(f(3,:),f(4,:).*f(16,:),'color',[82/255, 55/255, 156/255],'linewidth',2)
hold on
plot(f(3,:),f(4,:).*f(16,:)+f(4,:).*f(15,:),'color',[26/255, 174/255, 185/255],'linewidth',2)
plot(f(3,:),f(4,:),'color',[226/255, 97/255, 0/255],'linewidth',2)
xlim([0 15])
ylim([0 0.3])

%% Timeline over 10 years
figure (2)
clf

%plot([0.25 0.25],[-6.5 6.5],'color',[55/255, 53/255, 129/255],'linewidth',1)
hold on
%plot([1/12 1/12],[-6.5 6.5],'color',[251/255, 72/255, 37/255],'linewidth',1)
%plot([3/365 3/365],[-6.5 6.5],'color',[241/255, 192/255, 11/255],'linewidth',1)
plot([0 0.5],[0 0],'color',[0 0 0],'linewidth',1)

load('4.mat')
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
plot(time,F_TAi,'-','color',[55/255, 53/255, 129/255],'linewidth',2)
plot(time,F_DICi,'--','color',[55/255, 53/255, 129/255],'linewidth',2)
plot(time,F_O2i,':','color',[55/255, 53/255, 129/255],'linewidth',2)
%plot(time,F_H2Si)

load('13.mat')
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
plot(time,F_TAi,'-','color',[251/255, 72/255, 37/255],'linewidth',2)
plot(time,F_DICi,'--','color',[251/255, 72/255, 37/255],'linewidth',2)
plot(time,F_O2i,':','color',[251/255, 72/255, 37/255],'linewidth',2)
%plot(time,F_H2Si)

load('22.mat')
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
plot(time,F_TAi,'-','color',[241/255, 192/255, 11/255],'linewidth',2)
plot(time,F_DICi,'--','color',[241/255, 192/255, 11/255],'linewidth',2)
plot(time,F_O2i,':','color',[241/255, 192/255, 11/255],'linewidth',2)
%plot(time,F_H2Si)

xlim([0 0.5])
ylim([-6.5 6.5])

figure (3)
clf

%plot([0.25 0.25],[-6.5 6.5],'color',[55/255, 53/255, 129/255],'linewidth',1)
hold on
%plot([1/12 1/12],[-6.5 6.5],'color',[251/255, 72/255, 37/255],'linewidth',1)
%plot([3/365 3/365],[-6.5 6.5],'color',[241/255, 192/255, 11/255],'linewidth',1)
plot([0 0.5],[0 0],'color',[0 0 0],'linewidth',1)

load('31.mat')
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
plot(time,F_TAi,'-','color',[55/255, 53/255, 129/255],'linewidth',2)
plot(time,F_DICi,'--','color',[55/255, 53/255, 129/255],'linewidth',2)
plot(time,F_O2i,':','color',[55/255, 53/255, 129/255],'linewidth',2)
%plot(time,F_H2Si)

load('40.mat')
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
plot(time,F_TAi,'-','color',[251/255, 72/255, 37/255],'linewidth',2)
plot(time,F_DICi,'--','color',[251/255, 72/255, 37/255],'linewidth',2)
plot(time,F_O2i,':','color',[251/255, 72/255, 37/255],'linewidth',2)
%plot(time,F_H2Si)

load('49.mat')
F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./dbl)';
F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./dbl)';
F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./dbl)';
plot(time,F_TAi,'-','color',[241/255, 192/255, 11/255],'linewidth',2)
plot(time,F_DICi,'--','color',[241/255, 192/255, 11/255],'linewidth',2)
plot(time,F_O2i,':','color',[241/255, 192/255, 11/255],'linewidth',2)
%plot(time,F_H2Si)

xlim([0 0.5])
ylim([-6.5 6.5])

%% Fpoc=8, depth/time profiles

figure (4)

subplot 331
pcolor(time,-depths,dtalkf)
shading flat 
colorbar
xlim([0 1])
title('TA (mol/m3)')

subplot 332
pcolor(time,-depths,dFef)
shading flat 
colorbar
xlim([0 1])
title('Fe (mol/m3)')

subplot 333
pcolor(time,-depths,dCaf)
shading flat 
colorbar
xlim([0 1])
title('Ca (mol/m3)')

subplot 334
pcolor(time,-depths,dMnf)
shading flat 
colorbar
xlim([0 1])
title('Mn (mol/m3)')

subplot 335
pcolor(time,-depths,dO2f)
shading flat 
colorbar
xlim([0 1])
title('O2 (mol/m3)')

subplot 336
pcolor(time,-depths,dtCO2f)
shading flat 
colorbar
xlim([0 1])
title('DIC (mol/m3)')

subplot 337
pcolor(time,-depths,dtH2Sf)
shading flat 
colorbar
xlim([0 1])
title('H2S (mol/m3)')

subplot 338
pcolor(time,-depths,dtNH4f)
shading flat 
colorbar
xlim([0 1])
title('NH4 (mol/m3)')

subplot 339
pcolor(time,-depths,dtNO3f)
shading flat 
colorbar
xlim([0 1])
title('NO3 (mol/m3)')

%%
figure (5)
subplot 331
pcolor(time,-depths,dtPO4f)
shading flat 
colorbar
xlim([0 1])
title('PO4 (mol/m3)')

subplot 332
pcolor(time,-depths,dtSO4f)
shading flat 
colorbar
xlim([0 1])
title('SO4 (mol/m3)')

subplot 333
pcolor(time,-depths,fdCH4f)
shading flat 
colorbar
xlim([0 1])
title('fraction POC degraded into CH4')

subplot 334
pcolor(time,-depths,fdO2f)
shading flat 
colorbar
xlim([0 1])
title('fraction POC degraded by O2')

subplot 335
pcolor(time,-depths,fdtNO3f)
shading flat 
colorbar
xlim([0 1])
title('fraction POC degraded by NO3')

subplot 336
pcolor(time,-depths,fdtSO4f)
shading flat 
colorbar
xlim([0 1])
title('fraction POC degraded by SO4')

subplot 337
pcolor(time,-depths,fpFeOH3f)
shading flat 
colorbar
xlim([0 1])
title('fraction POC degraded by FeOH3')

subplot 338
pcolor(time,-depths,fpMnO2f)
shading flat 
colorbar
xlim([0 1])
title('fraction POC degraded by MnO2')

subplot 339
pcolor(time,-depths,OmegaCf)
shading flat 
colorbar
xlim([0 1])
title('Omega calcite')

%%
figure (6)

subplot 331
pcolor(time,-depths,pcalcitef)
shading flat 
colorbar
xlim([0 1])
title('calcite (mol/m3)')

subplot 332
pcolor(time,-depths,pclayf)
shading flat 
colorbar
xlim([0 1])
title('clay (mol/m3)')

subplot 333
pcolor(time,-depths,pFeOH3f)
shading flat 
colorbar
xlim([0 1])
title('FeOH3 (mol/m3)')

subplot 334
pcolor(time,-depths,pfocf)
shading flat 
colorbar
xlim([0 1])
title('fast-decay natural POC (mol/m3)')

subplot 335
pcolor(time,-depths,pfocsf)
shading flat 
colorbar
xlim([0 1])
title('fast-decay seaweed POC (mol/m3)')

subplot 336
pcolor(time,-depths,psocf)
shading flat 
colorbar
xlim([0 1])
title('slow-decay natural POC (mol/m3)')

subplot 337
pcolor(time,-depths,psocsf)
shading flat 
colorbar
xlim([0 1])
title('slow-decay seaweed POC (mol/m3)')

subplot 338
pcolor(time,-depths,procf)
shading flat 
colorbar
xlim([0 1])
title('refractory natural POC (mol/m3)')

subplot 339
pcolor(time,-depths,procsf)
shading flat 
colorbar
xlim([0 1])
title('refractory seaweed POC (mol/m3)')

%%
figure (7)
subplot 331
pcolor(time,-depths,pMnO2f)
shading flat 
colorbar
xlim([0 1])
title('MnO2 (mol/m3)')

subplot 332
pcolor(time,-depths,Rd_calcitef)
shading flat 
colorbar
xlim([0 1])
ylim([-0.2 0])
title('rate diss calcite (mol/m3/a)')

subplot 333
pcolor(time,-depths,Rf_totf)
shading flat 
colorbar
xlim([0 1])
title('respiration rate natural POC (mol/m3/a)')

subplot 334
pcolor(time,-depths,Rf_totsf)
shading flat 
colorbar
xlim([0 1])
title('respiration rate seaweed POC (mol/m3/a)')

subplot 335
pcolor(time,-depths,Rp_calcitef)
shading flat 
colorbar
xlim([0 1])
title('rate prec calcite (mol/m3/a)')

subplot 336
pcolor(time,-depths,uf)
shading flat 
colorbar
xlim([0 1])
title('u (m/a)')

subplot 337
pcolor(time,-depths,wf)
shading flat 
colorbar
xlim([0 1])
title('w (m/a)')

%%

kfastsf=NaN(201,10000);
kslowsf=NaN(201,10000);
kfastf=NaN(201,10000);
kslowf=NaN(201,10000);
for i=1:10000
    kfastsf(:,i)=kfasts;
    kfastf(:,i)=kfast;
    kslowsf(:,i)=kslows;
    kslowf(:,i)=kslow;
end

%%

figure (8)

subplot 221
pcolor(time,-depths,log10(Rd_calcitef))
hold on
contour(time,-depths,OmegaCf,[1.0001 1.0001],'color',[1 1 1],'linewidth',2)
shading flat 
colorbar
xlim([0 10])
ylim([-0.4 0])
cmocean('haline')
caxis([-1 3])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

subplot 222
pcolor(time,-depths,log10(Rp_calcitef))
hold on
contour(time,-depths,OmegaCf,[1.0001 1.0001],'color',[1 1 1],'linewidth',2)
shading flat 
colorbar
xlim([0 1])
cmocean('haline')
xlim([0 10])
ylim([-0.4 0])
caxis([-1 2])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

subplot 223
pcolor(time,-depths,log10(fdO2f.*(psocf.*kslowf+pfocf.*kfastf+psocsf.*kslowsf+pfocsf.*kfastsf)))
pcolor(time,-depths,log10(mouloud))
shading flat 
colorbar
xlim([0 1])
cmocean('thermal')
xlim([0 4])
caxis([-2 4])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

subplot 224
pcolor(time,-depths,log10((fpMnO2f+fpFeOH3f+fdCH4f+fdtNO3f+fdtSO4f).*(psocf.*kslowf+pfocf.*kfastf+psocsf.*kslowsf+pfocsf.*kfastsf)))
shading flat 
colorbar
xlim([0 1])
cmocean('thermal')
xlim([0 4])
caxis([-2 4])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% subplot 313
% pcolor(time,-depths,R(fdtNO3f+fdCH4f+fdtSO4f+fpMnO2f+fpFeOH3f))
% shading flat 
% colorbar
% xlim([0 1])
% %ylim([-0.05 0])

%%

figure (9)
clf

subplot 321
plot(dtalkic,-depths,'color',[55/255, 53/255, 129/255],'linewidth',2)
hold on
plot(dtalkf(:,251),-depths,'color',[251/255, 72/255, 37/255],'linewidth',2)
plot(dtalkf(:,500),-depths,'color',[241/255, 192/255, 11/255],'linewidth',2)
ylabel('depth (m)')
xlabel('TA (mol/m3)')
xlim([2 7])

subplot 322
plot(dtCO2ic,-depths,'color',[55/255, 53/255, 129/255],'linewidth',2)
hold on
plot(dtCO2f(:,251),-depths,'color',[251/255, 72/255, 37/255],'linewidth',2)
plot(dtCO2f(:,500),-depths,'color',[241/255, 192/255, 11/255],'linewidth',2)
ylabel('depth (m)')
xlabel('DIC (mol/m3)')
xlim([2 6])

subplot 323
plot(dO2ic,-depths,'color',[55/255, 53/255, 129/255],'linewidth',2)
hold on
plot(dO2f(:,251),-depths,'color',[251/255, 72/255, 37/255],'linewidth',2)
plot(dO2f(:,500),-depths,'color',[241/255, 192/255, 11/255],'linewidth',2)
ylabel('depth (m)')
xlabel('O2 (mol/m3)')
xlim([0 0.3])

subplot 324
plot(pcalciteic,-depths,'color',[55/255, 53/255, 129/255],'linewidth',2)
hold on
plot(pcalcitef(:,251),-depths,'color',[251/255, 72/255, 37/255],'linewidth',2)
plot(pcalcitef(:,500),-depths,'color',[241/255, 192/255, 11/255],'linewidth',2)
ylabel('depth (m)')
xlabel('calcite (mol/m3)')
xlim([0 5500])

subplot 325
plot(dtNO3ic,-depths,'color',[55/255, 53/255, 129/255],'linewidth',2)
hold on
plot(dtNO3f(:,251),-depths,'color',[251/255, 72/255, 37/255],'linewidth',2)
plot(dtNO3f(:,500),-depths,'color',[241/255, 192/255, 11/255],'linewidth',2)
ylabel('depth (m)')
xlabel('NO3 (mol/m3)')
xlim([0 0.04])