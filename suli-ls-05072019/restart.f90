!!
SUBROUTINE restart_ini()
!!so deve rodar caso seja restart
use velpre
use obst
use tempo
use ls_param
use mms_m
use smag
implicit none
open (unit=6662, action= 'read', file= 'dados//arquivorestart', status= 'old', form='formatted')
 read(6662,*)  u,v,w,ku,kv,kw,prd0,prd1,prd,rho,ls_nu,ls,mod_ls,bxx0,bxx1,bxy0,bxz0,bxxf,bxxf1,bxyf,bxzf,byx0,byy0,&
&byy1,byz0,byxf,byyf,byyf1,byzf,bzx0,bzy0,bzz0,bzz1,bzxf,bzyf,bzzf,bzzf1,ub,vb,wb,d_max,d_min,ls_m,&
hsx,hsy,hsz,hs,kurv,xnut,ynut,znut,tf_u,tf_v,tf_w,vol_ini,vol_ins
 close (unit=6662)

!inicialização do código
it = 0
t=0.
cont = 10

! Método da integração temporal
write(*,*) "Esquema temporal"
if ((t_tempo == 0)  .or. (t_tempo == 3)) then
write(*,*) "Euler"
ntt = 1

a_dt = 1.0*dt

elseif (t_tempo == 1) then
write(*,*) "RK2"
ntt = 2

a_dt(1) = 0.5*dt
a_dt(2) = 1.0*dt
a_dt(3) = 1.0*dt

elseif (t_tempo == 2) then
write(*,*) "RK3"
ntt = 3

a_dt(1) = 0.5*dt
a_dt(2) = 1.0*dt
a_dt(3) = 1.0*dt

endif
if (mms_t == 0) then
	   tf_p = 0.
	   tf_u = 0.
	   tf_v = 0.
	   tf_w = 0.
endif
if (wave_t > 0) then
	write(*,*) "Com entrada de onda"
	CALL waves_coef()
else
	write(*,*) "Sem entrada de onda"
endif

alpha1 = 1.5 ! número de células que varão parte da espessura, pois é uma função suave

rho_f1 = 1.204!kg/m³ ar (ls negativo) 20°C
mi_f1 = 0.000018253 !Pa/s  !0.00001516!m²/s ar (ls negativo) 20°C !!!$$$$$$$$$ no incompact3d tá como NI!
rho_f2 =998.0!kg/m³ água saturada (ls positivo) 20°C.
mi_f2 = 0.00100798  !Pa/s  !0.00000101!m²/s água saturada (ls positivo) 20°C !!!$$$$$$$$$$ no incompact3d tá como NI!

sigma = 0.0728!N/m tensão superficial da água 20°C
	
rho_m = abs(rho_f2-rho_f1)*0.5

! coeficientes de integração RK3 TVD

 adtl(1)=1.
 bdtl(1)=0.
 gdtl(1)=1.

 adtl(2)=3./4.
 bdtl(2)=1./4.
 gdtl(2)=1./4.

 adtl(3)=1./3.
 bdtl(3)=2./3.
 gdtl(3)=2./3.





END SUBROUTINE restart_ini

SUBROUTINE restart_salva()
use velpre
use obst
use ls_param
use smag
use mms_m
implicit none
!u,v,w
!ku,kv,kw
!prd0,prd1,prd,rho,ls_nu,ls,mod_ls
open (unit=6661, action= 'write', file= 'dados//arquivorestart', status= 'unknown', form='formatted')
 write(6661,*) u,v,w,ku,kv,kw,prd0,prd1,prd,rho,ls_nu,ls,mod_ls,bxx0,bxx1,bxy0,bxz0,bxxf,bxxf1,bxyf,bxzf,byx0,byy0,&
&byy1,byz0,byxf,byyf,byyf1,byzf,bzx0,bzy0,bzz0,bzz1,bzxf,bzyf,bzzf,bzzf1,ub,vb,wb,d_max,d_min,ls_m,&
hsx,hsy,hsz,hs,kurv,xnut,ynut,znut,tf_u,tf_v,tf_w,vol_ini,vol_ins
 close (unit=6661)
print *, 'salvou restart'


END SUBROUTINE restart_salva
