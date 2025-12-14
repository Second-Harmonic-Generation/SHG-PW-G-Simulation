! Note: 
!     We have aligned the code with VS Code's formatting standards.
!     But it may appears cluttered on GitHub due to formatting differences.


!            ************************************************************************************
!            *                                                                                  *
!            * File name:                                                                       *
!            *     Code_SHG-PW-G-Simulation.F90                                                 *
!            *                                                                                  *
!            * This Fortran code models pulsed wave Gaussian Second Harmonic Generation         *
!            * using five coupled equations (Heat-equation & Phase-equation & Fields-equations) *
!            * solved by Finite Difference Method.                                              *
!            *                                                                                  *
!            ************************************************************************************

program Coupled_G_PW

implicit none

!*********************************************************************************************************************
!                                       Variables Definition
!*********************************************************************************************************************

!--------------------------------- common variables
integer       i            ,j          ,k          ,l                                                                &
             ,nt           ,nr         ,nz         ,Np         ,inn        ,kn                                       &
	          ,nt1          ,run                                                                                      &
	          ,nomegaf                                                          

real*8        E            ,t          ,z          ,r          ,x          ,y                                        &                                 
             ,pi           ,tp         ,y1         ,y2                                                               &
             ,freq                                                                                                   &
	          ,timet        ,alpha      ,gama1      ,gama2      ,gama3      ,power      ,nnrom                        &   
	          ,no1T0        ,ne1T0      ,ne2T0      ,no1rT      ,ne2rT      ,ne1rT                                    &
	          ,omegaf       ,length     ,deltar     ,deltaz     ,deltat     ,radius                                   &  
	          ,lambda1      ,lambda2    ,deltar1    ,deltar2                                                          & 
             ,tbetween                                                             

complex*16    Ii                                                                                                       

character*35  freqf      ,Npf        ,tpf        ,EE

!--------------------------------- temperature variables
real*8        h                                                                                                      &  
             ,T0         ,Cp                                                                                         &
	          ,roh        ,aa1        ,aa2        ,aa3        ,aa4        ,aa5        ,KT0                            &
             ,Tinf       ,Tamb       ,Temp                                                                           & 
	          ,sigma                                                                                                  &
	          ,Tempmax                                                                                                &
             ,epsilong                                                                                               &
	          ,stability                                                                                              &
		 	  
	          ,temperature[allocatable](:,:,:)    ,KT[allocatable](:,:)                                              

character*35  filenameTt   ,filenameTr   ,filenameTz                                                                 &
             ,filenameTempmaxl
!--------------------------------- phase variables
real*8        phi                                                                                                    &
             ,B1T0         ,B2T0         ,C1T0          ,C2T0         ,B1rT         ,B2rT         ,C1rT              &
	          ,C2rT                                                                                                   &		 
             ,Phase                                                                                                  &
	          ,aa1T0        ,bb1T0        ,cc1T0         ,nx1T0        ,ny1T0        ,nz1T0        ,Term1             &
	          ,aa2T0        ,bb2T0        ,cc2T0         ,nx2T0        ,ny2T0        ,nz2T0        ,Term2             &
	          ,aa1rT        ,bb1rT        ,cc1rT         ,nx1rT        ,ny1rT        ,nz1rT        ,Term3             &
	          ,aa2rT        ,bb2rT        ,cc2rT         ,nx2rT        ,ny2rT        ,nz2rT        ,theta             &
	          ,B1r0T        ,B2r0T        ,C1r0T         ,C2r0T                                                       &                                

	          ,aa1r0T       ,bb1r0T       ,cc1r0T        ,nx1r0T       ,ny1r0T       ,nz1r0T       ,dnx1dT            &
	          ,aa2r0T       ,bb2r0T       ,cc2r0T        ,nx2r0T       ,ny2r0T       ,nz2r0T       ,dnx2dT            &
             ,dny1dT       ,dny2dT       ,dnz1dT        ,dnz2dT       ,no1r0T       ,ne1r0T       ,ne2r0T            &   
             
             ,Phasemin                                                                                               &

	          ,deltano1rT   ,deltane1rT   ,deltane2rT                                                                 &   
	          ,deltano1r0T  ,deltane1r0T  ,deltane2r0T                                                     

complex*8     deltaphase[allocatable](:,:)                                                                           &
             ,phasechange[allocatable](:,:,:)

character*35  filenamePt   ,filenamePr   ,filenamePz                                                                 &
             ,filenamePhaseminl

!--------------------------------- fields variables
integer       f            ,ibest

real*8        c                                                                                                      &
             ,fi                                                                                                     &
	          ,deff                                                                                                   &
             ,omega        ,Psi22         ,Psi32                                                                     &
             ,Lscale       ,Elec12        ,Elec22         ,Elec32                                                    &
	         ,epsilon0     ,Psi2max       ,Psi3max                                                                                          

complex*16    cc1          ,cc2          ,cc3           ,cc4          ,cc5                                           &
             ,dd1          ,dd2          ,dd3           ,dd4          ,dd5                                           &
             ,ee1          ,ee2          ,ee3           ,ee4          ,ee5                                           &

             ,Psi1[allocatable](:,:,:)   ,Elec1[allocatable](:,:,:)                                                  &
	          ,Psi2[allocatable](:,:,:)   ,Elec2[allocatable](:,:,:)                                                  &  
	         ,Psi3[allocatable](:,:,:)   ,Elec3[allocatable](:,:,:) 
                     

character*35  filenameibestl                                                                                         &
             ,filenameElec12t      ,filenameElec12r       ,filenameElec12z                                           &
             ,filenameElec22t      ,filenameElec22r       ,filenameElec22z                                           &         
	          ,filenameElec32t      ,filenameElec32r       ,filenameElec32z                                           &
             
	          ,filenamePsi3picksl                                                                                     & 
             ,filenamePsi2picksl 

!**********************************************************************************************************************
!                                    Giving Zero to variables
!**********************************************************************************************************************

!--------------------------------- Giving Zero to common variables
              i = 0            ;j = 0            ;k = 0        ;l = 0               
             nt = 0           ;nr = 0           ;nz = 0       ;Np = 0       ;inn= 0        ;kn = 0                                                  
	         nt1 = 0           ;run = 0                                                 
        nomegaf = 0

              E = 0.            ;t = 0.            ;z = 0.        ;r = 0.        ;x = 0.         ;y = 0.         
             pi = 0.           ;tp = 0.           ;y1 = 0.       ;y2 = 0.           
           freq = 0.                                                                                                
	       timet = 0.        ;alpha = 0.        ;gama1 = 0.    ;gama2 = 0.    ;gama3 = 0.     ;power = 0.    ;nnrom = 0.
          no1T0 = 0.        ;ne2T0 = 0.        ;ne1T0 = 0.    ;no1rT = 0.    ;ne2rT = 0.     ;ne1rT = 0.      
	      omegaf = 0.       ;length = 0.       ;deltar = 0.   ;deltaz = 0.   ;deltat = 0.    ;radius = 0.                                               
      	lambda1 = 0.      ;lambda2 = 0.      ;deltar1 = 0.  ;deltar2 = 0.
       tbetween = 0.                                              

             Ii = (0.,0.)  

!-------------------------------- Giving Zero to temperature variable
              h = 0.                                                                                                                       
             T0 = 0.           ;Cp = 0.                                                    
	         roh = 0.          ;aa1 = 0.          ;aa2 = 0.      ;aa3 = 0.      ;aa4 = 0.      ;aa5 = 0.      ;KT0 = 0.      
           Tinf = 0.         ;Tamb = 0.         ;Temp = 0.
	       sigma = 0.                                                                       
        Tempmax = 0. 
       epsilong = 0.                                      
      stability = 0.   

!------------------------------- Giving Zero to phase variable
            phi = 0.                                                                                                    
           B1T0 = 0.         ;B2T0 = 0.         ;C1T0 = 0.     ;C2T0 = 0.     ;B1rT = 0.     ;B2rT = 0.     ;C1rT = 0.        
           C2rT = 0.    
          Phase = 0.			 
          aa1T0 = 0.        ;bb1T0 = 0.        ;cc1T0 = 0.    ;nx1T0 = 0.    ;ny1T0 = 0.    ;nz1T0 = 0.    ;Term1 = 0.         
	       aa2T0 = 0.        ;bb2T0 = 0.        ;cc2T0 = 0.    ;nx2T0 = 0.    ;ny2T0 = 0.    ;nz2T0 = 0.    ;Term2 = 0.        
	       aa1rT = 0.        ;bb1rT = 0.        ;cc1rT = 0.    ;nx1rT = 0.    ;ny1rT = 0.    ;nz1rT = 0.    ;Term3 = 0.              
          aa2rT = 0.        ;bb2rT = 0.        ;cc2rT = 0.    ;nx2rT = 0.    ;ny2rT = 0.    ;nz2rT = 0.    ;theta = 0.       
	       B1r0T = 0.        ;B2r0T = 0.        ;C1r0T = 0.    ;C2r0T = 0.                                                                     
                          
            
	      aa1r0T = 0.       ;bb1r0T = 0.       ;cc1r0T = 0.   ;nx1r0T = 0.   ;ny1r0T = 0.   ;nz1r0T = 0.   ;dnx1dT = 0.       
	      aa2r0T = 0.       ;bb2r0T = 0.       ;cc2r0T = 0.   ;nx2r0T = 0.   ;ny2r0T = 0.   ;nz2r0T = 0.   ;dnx2dT = 0.     
         dny1dT = 0.       ;dny2dT = 0.       ;dnz1dT = 0.   ;dnz2dT = 0.   ;no1r0T = 0.   ;ne1r0T = 0.   ;ne2r0T = 0.    
     deltano1rT = 0.   ;deltane1rT = 0.   ;deltane2rT = 0.                                                                   
    deltano1r0T = 0.  ;deltane1r0T = 0.  ;deltane2r0T = 0.                                                    

!-------------------------------- Giving Zero to field equation variable
              f = 0.        ;ibest = 0.
              c = 0.                                                    
             fi = 0.         

            cc1 = 0.          ;cc2 = 0.          ;cc3 = 0.      ;cc4 = 0.      ;cc5 = 0.                                           
            dd1 = 0.          ;dd2 = 0.          ;dd3 = 0.      ;dd4 = 0.      ;dd5 = 0.                                          
            ee1 = 0.          ;ee2 = 0.          ;ee3 = 0.      ;ee4 = 0.      ;ee5 = 0.                                          
 
           deff = 0.                                        
	       omega = 0.        ;Psi22 = 0.        ;Psi32 = 0.
	      Lscale = 0.       ;Elec12 = 0.       ;Elec22 = 0.   ;Elec32 = 0.  
       epsilon0 = 0.      ;Psi2max = 0.      ;Psi3max = 0.

!**********************************************************************************************************************
!                                             Inputs		  
!**********************************************************************************************************************

! Note: 
!     This code lets the user enter values twice: once numerically (for calculations) 
!     and once as a string (for filenames or labels).  
!     For example, `E` is a number, while `EE` stores the same value as a string.  
!     This dual input ensures accurate calculations and meaningful file naming.

!write(*,'(/,2x,a,\)') '            Enter the Energy value  : '
!read(*,*) E
!write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) EE

!write(*,'(/,2x,a,\)') '         Enter the frequency value  : '
!read(*,*) freq
!write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) freqf

!write(*,'(/,2x,a,\)') '        Enter the Number of Pulses  : '
!read(*,*) Np
!write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) Npf

!write(*,'(/,2x,a,\)') '                      Enter the tp  : '
!read(*,*) tp
!write(*,'(/,2x,a,\)') '                             Again  : '
!read(*,*) tpf

! For Calculation
E = 0.45       
freq = 4000
Np = 1
tp = 50e-6

! For Generating Filenames based on the values above
EE = '045'            
freqf = '4000'
Npf = '1'
tpf = '50'

!**********************************************************************************************************************
!                          Determination of Filenames and Opening files
!**********************************************************************************************************************

! Note:
!      To achieve both efficiency and clarity in managing output data,
!      below, we generate filenames based on input information.

!------------------------------------  temperature 

!------------------------------------------------ Heat Equation Files
filenameTt = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Tt.plt'
open(1,file=filenameTt)
!write(1,'(/,a,/)')    ' variables=         "t"                             "temperature"'

filenameTr = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Tr.plt'
open(2,file=filenameTr)
!write(2,'(/,a,/)')    ' variables=         "r"                             "temperature"'

filenameTz = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Tz.plt'
open(3,file=filenameTz)
!write(3,'(/,a,/)')    ' variables=         "z"                             "temperature"' 

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameTt  &
                                                                                    ,filenameTr  &
																				                     	,filenameTz
 read(*,*)

!------------------------------------------------ Phase Equation Files
filenamePt = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Pt.plt'
open(4,file=filenamePt)
!write(4,'(/,a,/)')    ' variables=         "t"          "deltaphase_real"      "deltaphase_imaginary"'

filenamePr = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Pr.plt'
open(5,file=filenamePr)
!write(5,'(/,a,/)')    ' variables=         "r"          "deltaphase_real"      "deltaphase_imaginary"'

filenamePz = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Pz.plt'
open(6,file=filenamePz)
!write(6,'(/,a,/)')    ' variables=         "z"           "deltaphase_real"      "deltaphase_imaginary"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenamePt  &
                                                                                    ,filenamePr  &
	             								                                             ,filenamePz
 read(*,*)

!------------------------------------------------ Field Equations Files
filenameElec12t = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec12t.plt'
open(7,file=filenameElec12t)
!write(7,'(/,a,/)')    ' variables =     "t"                              "Elec1 ** 2"'

filenameElec12r = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec12r.plt'
open(8,file=filenameElec12r)
!write(8,'(/,a,/)')    ' variables =     "r"                              "Elec1 ** 2"'

filenameElec12z = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec12z.plt'
open(9,file=filenameElec12z)
!write(9,'(/,a,/)')    ' variables =     "z"                              "Elec1 ** 2"'


write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameElec12t  &
                                                                                    ,filenameElec12r  &
                                                                                    ,filenameElec12z
 read(*,*)

!------------------
filenameElec22t = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec22t.plt'
open(10,file=filenameElec22t)
!write(10,'(/,a,/)')    ' variables =     "t"                              "Elec2 ** 2"'

filenameElec22r = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec22r.plt'
open(11,file=filenameElec22r)
!write(11,'(/,a,/)')    ' variables =     "r"                              "Elec2 ** 2"'

filenameElec22z = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec22z.plt'
open(12,file=filenameElec22z)
!write(12,'(/,a,/)')    ' variables =     "z"                              "Elec2 ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameElec22t  &
                                                                                    ,filenameElec22r  &
										                                                      ,filenameElec22z
 read(*,*)

!------------------
filenameElec32t = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec32t.plt'
open(13,file=filenameElec32t)
!write(13,'(/,a,/)')    ' variables =     "t"                              "Elec3 ** 2"'

filenameElec32r = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec32r.plt'
open(14,file=filenameElec32r)
!write(14,'(/,a,/)')    ' variables =     "r"                              "Elec3 ** 2"'

filenameElec32z = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Elec32z.plt'
open(15,file=filenameElec32z)
!write(15,'(/,a,/)')    ' variables =     "z"                              "Elec3 ** 2"'

write(*,'(2/,a,/,40x,a,/,40x,a,/,40x,a,/)')' Results will be saved in these files :',filenameElec32t  &
                                                                                    ,filenameElec32r  &
										                                                      ,filenameElec32z
 read(*,*)

!------------------
filenamePsi3picksl = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Psi3picks_l.plt'
open(16,file=filenamePsi3picksl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "Psi3 picks ** 2"'

!------------------
filenamePsi2picksl = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Psi2picks_l.plt'
open(17,file=filenamePsi2picksl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "Psi2 picks ** 2"'

!------------------
filenameTempmaxl = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Tempmax_l.plt'
open(18,file=filenameTempmaxl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "Tempmax"'

!------------------
filenamePhaseminl= 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_Phasemin_l.plt'
open(19,file=filenamePhaseminl)
!write(18,'(/,a,/)')    ' variables =     "l"                              "Phasemin"'

!------------------
filenameibestl = 'E'//trim(EE)//'_f'//trim(freqf)//'_Np'//trim(Npf)//'_tp'//trim(tpf)//'_ibest_l.plt'
open(20,file=filenameibestl)
!write(16,'(/,a,/)')    ' variables =     "l"                              "ibest ** 2"'


!**********************************************************************************************************************
!                                           Constants
!**********************************************************************************************************************

!----------------- The constants of Common
       pi = 4*atan(1.)                                                              !dimensionless
       Ii = (0.,1.)

 tbetween = 1./freq              !4.*tp                                             !s
    timet = Np*tbetween                                                             !s
      !nt = 1600                 !int(1600./20.) !80   int(tbetween/deltat)         !dimensionless 
  !deltat = tbetween / nt                                                           !s     
      !inn= 20

       nz = 12000                                                                   !dimensionless
   length = 0.02                 !length of crystal                                 !m 
   deltaz = length/nz                                                               !m
       kn = 5

       nr = 120
   radius = 0.002                !radius of crystal                                 !m
   omegaf = 80.e-6              !spot size                                          !m
  nomegaf = 5 
    nnrom = 4./5.
  !deltar = omegaf/10.                                                              !m
  deltar1 = (nomegaf*omegaf)/(int(nnrom*nr))
  deltar2 = (radius-(nomegaf*omegaf))/(int((1-nnrom)*nr))
 !deltar1 = deltar
 !deltar2 = deltar1
      !nr = int(radius/deltar)                                                      !dimensionless 
 
       Cp = 728.016              !heat capacity at constant pressure                !J/(kg.K)
 
      KT0 = 13.                  !thermal conductivity of KTP crystal               !W/(m.K)
      roh = 2945.                !mass density                                      !kg/m^3

stability = 0.5
   deltat = stability * ( (roh*Cp)/(2.*KT0) ) * ( (deltar1**2.*deltaz**2.)/(deltar1**2.+deltaz**2.) ) !s   
      nt1 = int(tbetween/deltat)                                                    !dimensionless 
       inn= int(nt1/80)
       nt = nt1 - mod(nt1,inn)

    power = E/(sqrt(pi)*tp)                                                         !wat 

    gama1 = 0.5                  !the absorption coefficient of fundomental wave    !1/m
    gama2 = 0.5                  !the absorption coefficient of fundomental wave    !1/m
    gama3 = 4.                   !the absorption coefficient of fundomental SHW     !1/m        

   !no1T0 = 1.8296
   !ne1T0 = 1.7466   
   !ne2T0 = 1.7881   

  lambda1 = 1064.e-9             !wavelength fundamental                            !m
  lambda2 =  532.e-9             !wavelength second harmonic                        !m

!----------------- The constants of Temperature
        h = 10.                  !heat transfer coefficient (convection - cylinder) !W/(m^2.K)
       T0 = 300.                 !initial temperature                               !K

     Tinf = 300.
     Tamb = 300.
    sigma = 5.669e-8             !Stephan-Bultzman constant                         !W/(m^2.K^4) 

 epsilong = 0.9                  !surface emissivity                                !dimensionless

!------------------ The constants of phase equation
      phi = 24.77*pi/180.
    theta =   90.*pi/180.

   dnx1dT = (0.1323*(lambda1*1.e6)**(-3.) - 0.4385*(lambda1*1.e6)**(-2.) + 1.2307*(lambda1*1.e6)**(-1.) + 0.7709)*1.e-5
   dny1dT = (0.5014*(lambda1*1.e6)**(-3.) - 2.0030*(lambda1*1.e6)**(-2.) + 3.3016*(lambda1*1.e6)**(-1.) + 0.7498)*1.e-5
   dnz1dT = (0.3896*(lambda1*1.e6)**(-3.) - 1.3332*(lambda1*1.e6)**(-2.) + 2.2762*(lambda1*1.e6)**(-1.) + 2.1151)*1.e-5

   dnx2dT = (0.1323*(lambda2*1.e6)**(-3.) - 0.4385*(lambda2*1.e6)**(-2.) + 1.2307*(lambda2*1.e6)**(-1.) + 0.7709)*1.e-5
   dny2dT = (0.5014*(lambda2*1.e6)**(-3.) - 2.0030*(lambda2*1.e6)**(-2.) + 3.3016*(lambda2*1.e6)**(-1.) + 0.7498)*1.e-5
   dnz2dT = (0.3896*(lambda2*1.e6)**(-3.) - 1.3332*(lambda2*1.e6)**(-2.) + 2.2762*(lambda2*1.e6)**(-1.) + 2.1151)*1.e-5

    nx1T0 = sqrt(3.0065+0.03901/((lambda1*1.e6)**2.-0.04251)-0.01327*(lambda1*1.e6)**2.) 
    ny1T0 = sqrt(3.0333+0.04154/((lambda1*1.e6)**2.-0.04547)-0.01408*(lambda1*1.e6)**2.) 
        nz1T0 = sqrt(3.3134+0.05694/((lambda1*1.e6)**2.-0.05658)-0.01682*(lambda1*1.e6)**2.) 

        nx2T0 = sqrt(3.0065+0.03901/((lambda2*1.e6)**2.-0.04251)-0.01327*(lambda2*1.e6)**2.) 
        ny2T0 = sqrt(3.0333+0.04154/((lambda2*1.e6)**2.-0.04547)-0.01408*(lambda2*1.e6)**2.) 
        nz2T0 = sqrt(3.3134+0.05694/((lambda2*1.e6)**2.-0.05658)-0.01682*(lambda2*1.e6)**2.) 

        aa1T0 = 1. / nx1T0 ** 2.
    	bb1T0 = 1. / ny1T0 ** 2.
        cc1T0 = 1. / nz1T0 ** 2.

        aa2T0 = 1. / nx2T0 ** 2.
        bb2T0 = 1. / ny2T0 ** 2. 
	cc2T0 = 1. / nz2T0 ** 2.

        Term1 = sin(theta)**2. * cos(phi)**2.
        Term2 = sin(theta)**2. * sin(phi)**2.
        Term3 = cos(theta)**2.

         B1T0 = -Term1 * ( bb1T0 + cc1T0 )  &
	        -Term2 * ( aa1T0 + cc1T0 )       &
                -Term3 * ( aa1T0 + bb1T0 ) 
     		 
	 C1T0 = +Term1 * bb1T0 * cc1T0           &
	        +Term2 * aa1T0 * cc1T0           &
	        +Term3 * aa1T0 * bb1T0 


	 B2T0 = -Term1 * ( bb2T0 + cc2T0 )       &
	        -Term2 * ( aa2T0 + cc2T0 )       &
	        -Term3 * ( aa2T0 + bb2T0 )
             
	 C2T0 = +Term1 * bb2T0 * cc2T0           &
	        +Term2 * aa2T0 * cc2T0           &
	        +Term3 * aa2T0 * bb2T0 


 	no1T0 = (2.**0.5) / sqrt( -B1T0 - sqrt( B1T0**2. - 4.*C1T0 ) )  
        ne1T0 = (2.**0.5) / sqrt( -B1T0 + sqrt( B1T0**2. - 4.*C1T0 ) ) 
        ne2T0 = (2.**0.5) / sqrt( -B2T0 + sqrt( B2T0**2. - 4.*C2T0 ) ) 
   
!----------------- The field equations constants
           c = 3.e8                                                                 !m/s

       omega = 2.*pi*c/lambda1                                                      !rad/s 

    epsilon0 = 8.85e-12                                                             !C**2/N*m**2  or F(farad)/m

        deff = 7.3e-12           !nonlinear  effective coefficient                  !m/v       
	           
!**********************************************************************************************************************
!                                        Arrays Allocattion 
!**********************************************************************************************************************
 
!------------------- Allocate Array temperature
allocate(temperature(1:2,0:nr,0:nz/kn))     
allocate( KT(0:nr,0:nz/kn) )

!------------------- Allocate Array phase
allocate(deltaphase (0:nr,0:nz))
allocate(phasechange(0:nt/inn,0:nr,0:nz))                     

!------------------- Allocate Array field equations
allocate(Psi1(0:nt/inn,0:nr,1:2))                               
allocate(Psi2(0:nt/inn,0:nr,1:2))             
allocate(Psi3(0:nt/inn,0:nr,1:2))             

allocate(Elec1(0:nt/inn,0:nr,0:nz/kn))
allocate(Elec2(0:nt/inn,0:nr,0:nz/kn)) 
allocate(Elec3(0:nt/inn,0:nr,0:nz/kn))  

!**********************************************************************************************************************
!                                     Giving Zero to Arrays
!********************************************************************************************************************** 
  
!------------------- Zero to Array temperature
forall (i=1:2,j=0:nr,k=0:nz/kn)
                    temperature(i,j,k) = 0.        
end forall !i

!--------
forall (j=0:nr,k=0:nz/kn)
                               KT(j,k) = 0.
end forall

!-------------------- Zero to Array phase equation
forall (j=0:nr,k=0:nz)
                       deltaphase(j,k) = (0.,0.)        
end forall !i

!--------
forall (i=0:nt/inn,j=0:nr,k=0:nz)                            
                    phasechange(i,j,k) = (0.,0.)        
end forall !i

!--------------------- Zero to Array field equations
forall (i=0:nt/inn,j=0:nr,k=1:2)
                           Psi1(i,j,k) = (0.,0.)       
	  	           Psi2(i,j,k) = (0.,0.)       
                           Psi3(i,j,k) = (0.,0.)       
end forall !i

forall (i=0:nt/inn,j=0:nr,k=0:nz/kn)
                          Elec1(i,j,k) = (0.,0.)     
		          Elec2(i,j,k) = (0.,0.)      
                          Elec3(i,j,k) = (0.,0.)      
end forall !i

!**********************************************************************************************************************
!                                       Printing Constants     
!**********************************************************************************************************************

!------------------------------------------------ Common 
write(*,*)
write(*,*)'------- Common Constants ---------------------------------------------------'
write(*,*) 
write(*,'(A15,F15.10,//)')  '            E = ',E                

write(*,'(A15,2F15.10,/)')  '           Ii = ',Ii                  
write(*,'(A15,I5    ,/ )')  '           Nt = ',Nt
write(*,'(A15,I5    ,/ )')  '           inn= ',inn             
write(*,'(A15,I5    ,/ )')  '           Nr = ',Nr                
write(*,'(A15,I5    ,/ )')  '           Nz = ',Nz
write(*,'(A15,I5    ,/ )')  '           kn = ',kn 
write(*,'(A15,F15.10,//)')  '           pi = ',pi

write(*,'(A15,f15.10,/ )')  '        gama1 = ',gama1               
write(*,'(A15,f15.10,/ )')  '        gama2 = ',gama2               
write(*,'(A15,f15.10,//)')  '        gama3 = ',gama3              

write(*,'(A15,F15.10,/ )')  '        timet = ',timet           
write(*,'(A15,f15.10,//)')  '        power = ',power              

write(*,'(A15,F15.10,/ )')  '       omegaf = ',omegaf          
write(*,'(A15,F15.10,/ )')  '       length = ',length          
write(*,'(A15,F15.10,/ )')  '       deltat = ',deltat          
write(*,'(A15,F15.10,/ )')  '       deltaz = ',deltaz          
write(*,'(A15,F15.10,/ )')  '      deltar1 = ',deltar1          
write(*,'(A15,F15.10,//)')  '      deltar2 = ',deltar2          


write(*,'(A15,F15.10,//)')  '       radius = ',radius           

write(*,'(A15,F15.10,/ )')  '      lambda1 = ',lambda1         
write(*,'(A15,f15.10,//)')  '      lambda2 = ',lambda2        

write(*,'(A15,F15.10,//)')  '     tbetween = ',tbetween                                        

!------------------------------------------------ For Heat Equation 
write(*,*)
write(*,*)'------- Heat Equation Constants --------------------------------------------'
write(*,*)
write(*,'(A15,F15.10,/ )')  '            h = ',h               

write(*,'(A15,F15.10,/ )')  '           T0 = ',T0              
write(*,'(A15,F15.10,//)')  '           Cp = ',Cp              

write(*,'(A15,F15.10,/ )')  '          KT0 = ',KT0             
write(*,'(A15,F15.10,//)')  '          roh = ',roh              

write(*,'(A15,F15.10,//)')  '        sigma = ',sigma           

write(*,'(A15,F15.10,//)')  '     epsilong = ',epsilong        

write(*,'(A15,F15.10,//)')  '    stability = ',stability                                                                

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)

!------------------------------------------------ For Phase Equation 
write(*,*)
write(*,*)'------- Phase Equation Constants -------------------------------------------'
write(*,*)
write(*,'(A15,F15.10,//)')  '          phi = ',phi             

write(*,'(A15,f15.10,/ )')  '         B1T0 = ',B1T0                  
write(*,'(A15,F15.10,/ )')  '         B2T0 = ',B2T0            
write(*,'(A15,F15.10,/ )')  '         C1T0 = ',C1T0            
write(*,'(A15,F15.10,//)')  '         C2T0 = ',C2T0            

write(*,'(A15,F15.10,//)')  '        theta = ',theta           

write(*,'(A15,F15.10,/ )')  '        aa1T0 = ',aa1T0           
write(*,'(A15,F15.10,/ )')  '        bb1T0 = ',bb1T0           
write(*,'(A15,F15.10,/ )')  '        cc1T0 = ',cc1T0           
write(*,'(A15,F15.10,/ )')  '        aa2T0 = ',aa2T0           
write(*,'(A15,F15.10,/ )')  '        bb2T0 = ',bb2T0           
write(*,'(A15,F15.10,/ )')  '        cc2T0 = ',cc2T0           
write(*,'(A15,F15.10,/ )')  '        nx1T0 = ',nx1T0           
write(*,'(A15,F15.10,/ )')  '        ny1T0 = ',ny1T0           
write(*,'(A15,F15.10,/ )')  '        nz1T0 = ',nz1T0           
write(*,'(A15,F15.10,/ )')  '        nx2T0 = ',nx2T0           
write(*,'(A15,F15.10,/ )')  '        ny2T0 = ',ny2T0           
write(*,'(A15,F15.10,//)')  '        nz2T0 = ',nz2T0           

write(*,'(A15,F15.10,/ )')  '        no1T0 = ',no1T0           
write(*,'(A15,F15.10,/ )')  '        ne1T0 = ',ne1T0           
write(*,'(A15,F15.10,//)')  '        ne2T0 = ',ne2T0           

write(*,'(A15,F15.10,/ )')  '        Term1 = ',Term1           
write(*,'(A15,F15.10,/ )')  '        Term2 = ',Term2           
write(*,'(A15,F15.10,//)')  '        Term3 = ',Term3           

write(*,'(A15,F15.10,/ )')  '       dnx1dT = ',dnx1dT          
write(*,'(A15,F15.10,/ )')  '       dny1dT = ',dny1dT          
write(*,'(A15,F15.10,/ )')  '       dnz1dT = ',dnz1dT          
write(*,'(A15,F15.10,/ )')  '       dnx2dT = ',dnx2dT          
write(*,'(A15,F15.10,/ )')  '       dny2dT = ',dny2dT          
write(*,'(A15,F15.10,//)')  '       dnz2dT = ',dnz2dT          

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)

!------------------------------------------------ For fields Equation 
write(*,*)
write(*,*)'------- Field Equations Constants ------------------------------------------'
write(*,*)
write(*,'(A15,f15.3 ,/ )')  '            c = ',c                  

write(*,'(A15,f35.30,/ )')  '         deff = ',deff                

write(*,'(A15,f25.5 ,/ )')  '        omega = ',omega               

write(*,'(A15,f25.20,//)')  '     epsilon0 = ',epsilon0           

write(*,*)'----------------------------------------------------------------------------'
write(*,'(A,\)')' Please press any key to continue '
read(*,*)
  
!**********************************************************************************************************************
!                                   Main Block of the Program     
!**********************************************************************************************************************
!----------- Optimization
!write(*,*)

!run=run+1

!write(*,*) 'Run = ',run

!--------------------------------------------------------
do j=0,nr
   do k=0,nz
	  
      if (mod(k,kn)==0) then	  
         
		 temperature(1,j,k/kn) = T0
	                    KT(j,k/kn) = KT0
      end if
		   
   end do !k
end do !j	     

!-------------
no1rT = no1T0 
ne1rT = ne1T0 
ne2rT = ne2T0

!----------------------------------------------------------------- Run program for NP pulse 
do l=1,Np 
         
   !----------------------------------------------- Field
   !--------
   do k=0,nz
      z=k*deltaz 
       
      !--------
      do i=1,nt              
         t=(i-1)*deltat

         if (mod(i,inn)==0) then 

	    !--------	        
   	    do j=1,nr-1
                  
               if (j<=int(nnrom*nr)) then
                  r=j*deltar1
                  deltar=deltar1
                 else
                  r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
                  deltar=deltar2
               end if  
            
               !--------------------------------- Constants
               Lscale =  sqrt(no1rT*ne1rT*ne2rT)                                                &
			   
		       * sqrt( (epsilon0*c**3.*pi*omegaf**2.) / (4.*omega**2.*deff**2.*power) ) 

	       cc1 = deltaz * no1rT / c
               cc2 = deltaz *  Ii*c / (2.*no1rT*omega)
               cc3 = deltaz * gama1 / 2.
	       cc4 = deltaz *    Ii / Lscale

	       dd1 = deltaz * ne1rT / c
               dd2 = deltaz *  Ii*c / (2.*ne1rT*omega)
	       dd3 = deltaz * gama2 / 2.
	       dd4 = deltaz *    Ii / Lscale

               ee1 = deltaz * ne2rT / c
               ee2 = deltaz *  Ii*c / (4.*ne2rT*omega)
	       ee3 = deltaz * gama3 / 2.
	       ee4 = deltaz *    Ii / Lscale
	       !--------------------------------- Bounday conditions For Field Equations			

               !------------- Psi1
               if (k==0) Psi1(i/inn,j,1) = exp( (-(t-2.*tp)**2.)/(tp**2.) ) * exp(-r**2./omegaf**2.)!for input surface

	       Psi1(i/inn,0 ,1 ) = Psi1(i/inn,1,1 )                                                 !for crystal axis
           
               Psi1(i/inn,nr,1 ) = (0.,0.)                                                     !for lateral surface 
			

	       if (k==0 ) Psi1(i/inn,0 ,1 ) = Psi1(i/inn,1,1 )                                      !for (0 ,0 )
           
   	       if (k==0 ) Psi1(i/inn,nr,1 ) = (0.,0.)                                            !for (nr,0 ) 

	       if (k==nz) Psi1(i/inn,nr,1 ) = (0.,0.)                                               !for (nr,nz) 

	       if (k==nz) Psi1(i/inn,0 ,1 ) = Psi1(i/inn,1,1)                                       !for (0 ,nz) 

	       !------------- Psi2
               if (k==0) Psi2(i/inn,j,1) = exp( (-(t-2.*tp)**2.)/(tp**2.) ) * exp(-r**2./omegaf**2.)!for innput  surface

	       Psi2(i/inn,0 ,1 ) = Psi2(i/inn,1,1 )                                                 !for crystal axis
           
               Psi2(i/inn,nr,1 ) = (0.,0.)                                                     !for lateral surface 
          
	       if (k==0 ) Psi2(i/inn,0 ,1 ) = Psi2(i/inn,1,1 )                                      !for (0 ,0 )
           
	       if (k==0 ) Psi2(i/inn,nr,1 ) = (0.,0.)                                               !for (nr,0 ) 

               if (k==nz) Psi2(i/inn,nr,1) = (0.,0.)                                           !for (nr,nz) 

	       if (k==nz) Psi2(i/inn,0 ,1) = Psi2(i/inn,1 ,1)                                       !for (0 ,nz) 

 	       !------------- Psi3 
               if (k==0)Psi3(i/inn,j,1) = (0.,0.)                                              !for innput  surface
            
	       Psi3(i/inn,0  ,1) = Psi3(i/inn,1,1)                                                  !for crystal axis
					      
	       Psi3(i/inn,nr ,1) = (0.,0.)                                                          !for lateral surface 

	       if (k==0 ) Psi3(i/inn,0 ,1 ) = (0.,0.)                                               !for (0 ,0 )
           
	       if (k==0 ) Psi3(i/inn,nr,1 ) = (0.,0.)                                               !for (nr,0 ) 

	       if (k==nz) Psi3(i/inn,nr,1 ) = (0.,0.)                                               !for (nr,nz) 

	       if (k==nz) Psi3(i/inn,0 ,1 ) = Psi3(i/inn,1,1)                                       !for (0 ,nz) 
  
	       !--------------------------------- End of Bounday conditions
           
               !--------------------------------- Field Equations		   
	       !-------------
	       Psi1(i/inn,j,2) =  Psi1(i/inn,j,1)                                                                &
			
		               - cc1  * ( Psi1(i/inn,j,1) - Psi1(i/inn-1,j,1) ) / deltat                              &
			
		               + cc2  * ( Psi1(i/inn,j+1,1) - Psi1(i/inn,j-1,1) ) / (2*r*deltar)                      &
														
		    	       + cc2  * ( Psi1(i/inn,j+1,1) - 2*Psi1(i/inn,j,1) + Psi1(i/inn,j-1,1) ) / deltar**2       &
														
			       - cc3  *   Psi1(i/inn,j,1)                                                                  &

			       + cc4  * conjg(Psi2(i/inn,j,1)) * Psi3(i/inn,j,1) * exp(-Ii*phasechange(i/inn,j,k) )   

 	       !-------------
	       Psi2(i/inn,j,2) =  Psi2(i/inn,j,1)                                                                &
		
		               - dd1 * ( Psi2(i/inn,j,1) - Psi2(i/inn-1,j,1) ) / deltat                               &
			
		               + dd2 * ( Psi2(i/inn,j+1,1) - Psi2(i/inn,j-1,1) ) / (2*r*deltar)                       &
													
		               + dd2 * ( Psi2(i/inn,j+1,1) - 2*Psi2(i/inn,j,1) + Psi2(i/inn,j-1,1) ) / deltar**2      &
														
		               - dd3 *   Psi2(i/inn,j,1)                                                              &

		               + dd4 * conjg(Psi1(i/inn,j,1)) * Psi3(i/inn,j,1) * exp(-Ii*phasechange(i/inn,j,k) )     
 
               !-------------			
	       Psi3(i/inn,j,2) =  Psi3(i/inn,j,1)                                                                &
		
		               - ee1 * ( Psi3(i/inn,j,1) - Psi3(i/inn-1,j,1) ) / deltat                               &
			
		               + ee2 * ( Psi3(i/inn,j+1,1) -   Psi3(i/inn,j-1,1) ) / (2*r*deltar)                     &
														
		               + ee2 * ( Psi3(i/inn,j+1,1) - 2*Psi3(i/inn,j,1) + Psi3(i/inn,j-1,1) ) / deltar**2      &
														
			       - ee3 *   Psi3(i/inn,j,1)                                                                   &
														
			       + ee4 *   Psi1(i/inn,j,1) * Psi2(i/inn,j,1) * exp(Ii*phasechange(i/inn,j,k) )          

               !------------- Maximom 
	       if (k==nz) then
			      
	          Psi22 = Psi2(i/inn,0,1) * conjg(Psi2(i/inn,0,1)) * 100
                  Psi32 = Psi3(i/inn,0,1) * conjg(Psi3(i/inn,0,1)) * 100
		  
                  !------
                  if (Psi22 >= Psi2max) then

                     Psi2max = Psi22 
                     
                  end if

                  !------
                  if (Psi32 >= Psi3max) then

                     Psi3max = Psi32 
                     ibest = i

                  end if
				    
               end if

               !------------------------------                 
	       if (mod(k,kn)==0) then
                     
	          !------------- Elec1
                  if (k==0) Elec1(i/inn,j,0) = Psi1(i/inn,j,1)       !for input  surface
 
	          Elec1(i/inn,0 ,k/kn ) = Psi1(i/inn,0,1)                 !for crystal axis
           
	          Elec1(i/inn,nr,k/kn ) = (0.,0.)                         !for lateral surface 

                  if (k==0) Elec1(i/inn,0,0) = Psi1(i/inn,0,1)       !for (0 ,0 )
           
                  Elec1(i/inn,nr,0) = (0.,0.)                        !for (nr,0 ) 

                  Elec1(i/inn,nr,nz/kn) = (0.,0.)                    !for (nr,nz) 

                  if (k==nz)Elec1(i/inn,0 ,nz/kn) = Psi1(i/inn,0,1)  !for (0 ,nz) 

	          !------------- Elec2
                  if (k==0) Elec2(i/inn,j,0) = Psi2(i/inn,j,1)       !for input  surface
 
	          Elec2(i/inn,0 ,k/kn ) = Psi2(i/inn,0,1)                 !for crystal axis
           
	          Elec2(i/inn,nr,k/kn ) = (0.,0.)                         !for lateral surface 

                  if (k==0) Elec2(i/inn,0,0) = Psi2(i/inn,0,1)       !for (0 ,0 )
           
                  Elec2(i/inn,nr,0) = (0.,0.)                        !for (nr,0 ) 

                  Elec2(i/inn,nr,nz/kn) = (0.,0.)                    !for (nr,nz) 

                  if (k==nz) Elec2(i/inn,0 ,nz/kn) = Psi2(i/inn,0,1) !for (0 ,nz) 

	          !------------- Elec3 
                  Elec3(i/inn,j,0) = (0.,0.)                         !for input  surface
 
	          Elec3(i/inn,0 ,k/kn ) = Psi3(i/inn,1,1)                 !for crystal axis
           
	          Elec3(i/inn,nr,k/kn ) = (0.,0.)                         !for lateral surface 

                  Elec3(i/inn,0,0) = (0.,0.)                         !for (0 ,0 )
           
                  Elec3(i/inn,nr,0) = (0.,0.)                        !for (nr,0 ) 

                  Elec3(i/inn,nr,nz/kn) = (0.,0.)                    !for (nr,nz) 

                  if (k==nz) Elec3(i/inn,0,nz/kn) = Psi3(i/inn,1,2)  !for (0 ,nz) 

                  !---------------------
	          Elec1(i/inn,j,k/kn) = Psi1(i/inn,j,2)
		  Elec2(i/inn,j,k/kn) = Psi2(i/inn,j,2)
		  Elec3(i/inn,j,k/kn) = Psi3(i/inn,j,2)
              
	       end if

            end do !j
	 end if
      end do !i

      !-------------- End-Psi of each deltaz  ==> Initial Psi for next deltaz
      do i=0,nt   
         do j=1,nr-1
      	
            Psi1(i/inn,j,1) = Psi1(i/inn,j,2)
            Psi2(i/inn,j,1) = Psi2(i/inn,j,2)
            Psi3(i/inn,j,1) = Psi3(i/inn,j,2)

	     end do !j
       end do !i	     
       !-------------

   end do !k
   !----------------------------------------------- End of Field

   !----------------------------------------------- Heat_Phase
   !--------
   do i=0,nt

       !--------
      do j=1,nr-1
         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

	 !--------
         do k=1,nz-1
	    z=k*deltaz 
	       
	       if (mod(k,kn)==0) then
               
	          !-------- stability coefficient
                  !stability = ( (2.*KT0*deltat)/(roh*Cp) ) * ( (deltar**2.+deltaz**2.)/(deltar**2.*deltaz**2.) ) 

	          !--------------------------------- Constants
                  aa1 = (h*deltaz)/(KT(j,k/kn))
            
	          aa2 = (epsilong*sigma*deltaz)/(KT(j,k/kn))
            
	          aa3 = ( deltat/(roh*Cp) ) * KT(j,k/kn)          
              
	          aa4 = ( deltat/(roh*Cp) ) * power/(pi*omegaf**2)

	          aa5 = ( deltat/(roh*Cp) ) * (1./4.)
                  !------------------------------------ Boundary conditions
			
	          !-------------- For Heat Equation
	           
	          temperature(1,0 ,k/kn ) = temperature(1,1,k/kn)      !Thermal insulation condition for crystal axis

                  temperature(1,nr,k/kn ) = T0                         !Temperature-fixed condition for lateral surface
            
	          temperature(1,j ,0    ) = temperature(1,j,1)                        &  
			                       
					   - aa1*( temperature(1,j,1)    - Tinf )     &         !Convection & Radiation 
			                       
				           - aa2*( temperature(1,j,1)**4 - Tamb**4 )       !condition for input surface
            
	           
	          temperature(1,j ,nz/kn) = temperature(1,j,nz/kn-1)                  &
			   
			                   - aa1*( temperature(1,j,nz/kn-1) - Tinf )  &         !Convection & Radiation 

			                   - aa2*( temperature(1,j,nz/kn-1)**4 - Tamb**4 ) !condition for outputsurface
               
                  !---------------------
	          temperature(1,0 ,0    ) = temperature(1,0,1)                        &
			   
			                   - aa1*( temperature(1,0,1) - Tinf )        &         !Convection & Radiation 
			                                           
				           - aa2*( temperature(1,0,1)**4 - Tamb**4 )               !condition for (0,0)
		    
	           
	          temperature(1,0 ,nz/kn) = temperature(1,0,nz/kn-1)                  &
			   
			                   - aa1*( temperature(1,0,nz/kn-1) - Tinf )  &         !Convection & Radiation 
			                                                     
					   - aa2*( temperature(1,0,nz/kn-1)**4 - Tamb**4 )        !condition for (0,nz)

              
	          temperature(1,nr,0    ) = T0                                 !Temperature-fixed condition for (nr,0 )

	          temperature(1,nr,nz/kn) = T0                                 !Temperature-fixed condition for (nr,nz)
					
	          !------------------------------------ End of Boundary conditions
           
                  !------------------------------------ Heat Equation
	          temperature(2,j,k/kn) =                                                                        &
			   
		         + temperature(1,j,k/kn)                                                                      &
		                        
		         + aa3 * ( ( temperature(1,j+1,k/kn) -  temperature(1,j-1,k/kn) ) / (2*r*deltar)              &
								 
			          +( temperature(1,j+1,k/kn) -2*temperature(1,j,k/kn) + temperature(1,j-1,k/kn) )          &

                                     / (deltar**2)                                                                    & 

                                  +( temperature(1,j,k/kn-1) -2*temperature(1,j,k/kn) + temperature(1,j,k/kn+1) )     &

                                     / (deltaz**2) )                                                                  & 

                         + aa4 * (    gama1 * Elec1(i/inn,j,k/kn)*conjg(Elec1(i/inn,j,k/kn))                          &
									         
		                  +   gama2 * Elec2(i/inn,j,k/kn)*conjg(Elec2(i/inn,j,k/kn))                                    &
											 
			          + 2*gama3 * Elec3(i/inn,j,k/kn)*conjg(Elec3(i/inn,j,k/kn)) )                                       &

                         + aa5 * ( ( kT(j+1,k/kn)-kT(j-1,k/kn) ) * (temperature(1,j+1,k/kn)-temperature(1,j-1,k/kn) ) &
			                         
			             / (4*deltar**2)                                                                                 &

			          +( kT(j,k/kn+1)-kT(j,k/kn-1) ) * (temperature(1,j,k/kn+1)-temperature(1,j,k/kn-1))                 &
						   
			             / (4*deltaz**2) )

                  !------------- Maximom 
                  Temp = temperature(1,j,nz/kn)

                  if (Temp > Tempmax) then

                     Tempmax = Temp
                     
                  end if

	          !----------------------------------- Phase Equation constants

	          !----------------------------------- Phase Equation constants
  		      nx1r0T = nx1T0 + dnx1dT * ( temperature(1,0,k/kn) - T0 )
		      ny1r0T = ny1T0 + dny1dT * ( temperature(1,0,k/kn) - T0 )
		      nz1r0T = nz1T0 + dnz1dT * ( temperature(1,0,k/kn) - T0 )
			 
		      nx2r0T = nx2T0 + dnx2dT * ( temperature(1,0,k/kn) - T0 )
		      ny2r0T = ny2T0 + dny2dT * ( temperature(1,0,k/kn) - T0 )
		      nz2r0T = nz2T0 + dnz2dT * ( temperature(1,0,k/kn) - T0 )
          		   
	              aa1r0T = 1. / ( nx1r0T )**2 
                      bb1r0T = 1. / ( ny1r0T )**2 
                      cc1r0T = 1. / ( nz1r0T )**2 

                      aa2r0T = 1. / ( nx2r0T )**2 
                      bb2r0T = 1. / ( ny2r0T )**2 
                      cc2r0T = 1. / ( nz2r0T )**2       

		       B1r0T = -Term1 * ( bb1r0T + cc1r0T )                      &
			       -Term2 * ( aa1r0T + cc1r0T )                           &
			       -Term3 * ( aa1r0T + bb1r0T ) 
     		 
		       C1r0T =  Term1 * bb1r0T * cc1r0T                          &
			       +Term2 * aa1r0T * cc1r0T                               &
		               +Term3 * aa1r0T * bb1r0T 

                       B2r0T = -Term1 * ( bb2r0T + cc2r0T )            &
			       -Term2 * ( aa2r0T + cc2r0T )                           &
			       -Term3 * ( aa2r0T + bb2r0T )
             
		       C2r0T =  Term1 * bb2r0T * cc2r0T                          &
			       +Term2 * aa2r0T * cc2r0T                               &
			       +Term3 * aa2r0T * bb2r0T 

             
 	              no1r0T = (2**0.5) / sqrt( -B1r0T  - sqrt( B1r0T**2 - 4*C1r0T ) )  
                      ne1r0T = (2**0.5) / sqrt( -B1r0T  + sqrt( B1r0T**2 - 4*C1r0T ) )
                      ne2r0T = (2**0.5) / sqrt( -B2r0T  + sqrt( B2r0T**2 - 4*C2r0T ) ) 
			 
	         deltano1r0T = no1r0T - no1T0
                 deltane1r0T = ne1r0T - ne1T0
                 deltane2r0T = ne2r0T - ne2T0

	          !----------------------------------- Phase Equation constants
  		       nx1rT = nx1T0 + dnx1dT * ( temperature(1,j,k/kn) - T0 )
		       ny1rT = ny1T0 + dny1dT * ( temperature(1,j,k/kn) - T0 )
		       nz1rT = nz1T0 + dnz1dT * ( temperature(1,j,k/kn) - T0 )
			 
		       nx2rT = nx2T0 + dnx2dT * ( temperature(1,j,k/kn) - T0 )
	               ny2rT = ny2T0 + dny2dT * ( temperature(1,j,k/kn) - T0 )
	               nz2rT = nz2T0 + dnz2dT * ( temperature(1,j,k/kn) - T0 )
          		   
		       aa1rT = 1 / ( nx1rT )**2 
                       bb1rT = 1. / ( ny1rT )**2 
                       cc1rT = 1. / ( nz1rT )**2 

                       aa2rT = 1. / ( nx2rT )**2 
                       bb2rT = 1. / ( ny2rT )**2 
                       cc2rT = 1. / ( nz2rT )**2       

	                B1rT = -Term1 * ( bb1rT + cc1rT )                   &
			       -Term2 * ( aa1rT + cc1rT )                             &
			       -Term3 * ( aa1rT + bb1rT ) 
     		 
		        C1rT =  Term1 * bb1rT * cc1rT                            &
			       +Term2 * aa1rT * cc1rT                                 &
			       +Term3 * aa1rT * bb1rT 

                        B2rT = -Term1 * ( bb2rT + cc2rT )              &
			       -Term2 * ( aa2rT + cc2rT )                             &
		               -Term3 * ( aa2rT + bb2rT )
             
		        C2rT =  Term1 * bb2rT * cc2rT                            &
			       +Term2 * aa2rT * cc2rT                                 &
			       +Term3 * aa2rT * bb2rT 

             
 	               no1rT = (2**0.5) / sqrt( -B1rT  - sqrt( B1rT**2 - 4*C1rT ) )  
                       ne1rT = (2**0.5) / sqrt( -B1rT  + sqrt( B1rT**2 - 4*C1rT ) )
                       ne2rT = (2**0.5) / sqrt( -B2rT  + sqrt( B2rT**2 - 4*C2rT ) ) 
			 
	          deltano1rT = no1rT - no1T0
                  deltane1rT = ne1rT - ne1T0
                  deltane2rT = ne2rT - ne2T0
           
               end if			
               
               !------------------------------------ For Phase Equation
               deltaphase(j ,0 ) = (0.,0.)                                               !for input surface
            
               deltaphase(nr,k ) = (0.,0.)                                               !for lateral surface

	       deltaphase(j ,nz) = deltaphase(j,nz-1)                                 &       !for output surface
				  + ( 2*pi*deltaz / lambda1 )                                        &
	                          * ( deltano1rT + deltane1rT - 2*deltane2rT )                      

               !------
               deltaphase(nr,0 ) = (0.,0.)                                               !for (nr,0 )

	       deltaphase(nr,nz) = (0.,0.)                                                    !for (nr,nz)

	       deltaphase(0 ,0 ) = (0.,0.)                                                    !for ( 0,0 )

	       !------
	       deltaphase(0 ,k ) = deltaphase(0,k-1)                                  &       !for crystal axis
				  + ( 2*pi*deltaz / lambda1 )                                        &
				  * ( deltano1r0T  + deltane1r0T  - 2*deltane2r0T  )                                            


	       deltaphase(0 ,nz) = deltaphase(0,nz-1)                                 &       !for ( 0,nz)
				  + ( 2*pi*deltaz / lambda1 )                                        &
			          * ( deltano1rT + deltane1rT - 2*deltane2rT )                     
        
               !----------------------------------- Phase Equation
	       deltaphase(j,k ) = deltaphase(j,k-1)                                   &
				 + ( 2*pi*deltaz / lambda1 )                                         &  
				 * ( deltano1rT + deltane1rT  - 2*deltane2rT  )                                 

               !------------- Minimom 
               Phase = deltaphase(j,nz)  
 
               if (Phase .LE. Phasemin) then

                  Phasemin = Phase
                     
               end if

         end do !k
      end do !j
      !-------------------------------------- End of Heat-Phase
      
      !------------ End-temprature of each deltat  ==> Initial temperature for next deltat
      do j=1,nr-1
         do k=1,nz-1
 	  
            if (mod(k,kn)==0) temperature(1,j,k/kn) = temperature(2,j,k/kn)

         end do !k
      end do !j
   
      !------------
      do j=0,nr
         do k=0,nz
 	  
	     if (mod(k,kn)==0)    KT(j,k/kn) = KT0 * T0 / temperature(1,j,k/kn)

         end do !k
      end do !j	     
       
      !------------ deltaphaze(j,k)  ==>  phase-changes(i/inn,j,k)
      do j=0,nr
         do k=0,nz
 	  
            phasechange(i/inn,j,k) = deltaphase(j,k)      

         end do !k
      end do !j

      !============================================= Print Results for each deltat
      t=(l-1)*nt*deltat + i*deltat 

      !--------- For Heat Equation
      write(1,'(2x,f25.10,5x,f25.10)')  t , temperature(1,0,nz/kn)
      
      if (mod(i,inn)==0) then
	  
	 !------ For Phase Equation
         write(4,'(2x,f25.10,5x,2f25.10)') t , phasechange(i/inn,0,nz) 
 
         !----- For Field Equations
 	 Elec12 = Elec1(i/inn,0,nz/kn) * conjg(Elec1(i/inn,0,0    )) * 100
	 write( 7,'(2x,f25.10,5x,f25.10)') t , Elec12 

	 Elec22 = Elec2(i/inn,0,nz/kn) * conjg(Elec2(i/inn,0,nz/kn)) * 100      
	 write(10,'(2x,f25.10,5x,f25.10)') t , Elec22 

	 Elec32 = Elec3(i/inn,0,nz/kn) * conjg(Elec3(i/inn,0,nz/kn)) * 100      
	 write(13,'(2x,f25.10,5x,f25.10)') t , Elec32
      
      end if
      !--------------------------------------------- End of run for each deltat     
 
   end do !i 
   !-----

   !-----
   write(16,'(2x,I5,2x,f8.2)') l   ,Psi3max
   psi3max = 0.

   !-----
   write(17,'(2x,I5,2x,f8.2)') l   ,Psi2max  
   psi2max = 0.
   
   !-----
   write(18,'(2x,I5,2x,f8.2)') l   ,Tempmax  
   Tempmax = 0.

   !-----
   write(19,'(2x,I5,2x,f8.4)') l   ,Phasemin  
   Phasemin = 0.

   !-----
   write(20,'(2x,I5,2x,I5)') l  ,ibest 

end do !l
!-------------------------------------------------------------------- End of run for each Pulse

!**********************************************************************************************************************
!                                        Arrays Deallocattion 
!**********************************************************************************************************************
!----------------------------------- Deallocate Arrays Thermal
!deallocate(temperature)     
!deallocate(KT)

!----------------------------------- Deallocate Arrays phase
!deallocate(deltaphase)
!deallocate(phasechange)                     

!----------------------------------- Deallocate Arrys Fields
!deallocate(Psi1)                               
!deallocate(Psi2)             
!deallocate(Psi3)             

!deallocate(Elec1)
!deallocate(Elec2) 
!deallocate(Elec3)  

!---------------
!end do !omegaf
!end do !Length
!end do !freq
!end do !E

!**********************************************************************************************************************
!                                        Printing Results     
!**********************************************************************************************************************

!------------------------------------------------For heat equation
do j=0,nr
   
   if (j<=int(nnrom*nr)) then
      r=j*deltar1
      deltar=deltar1
     else
      r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
      deltar=deltar2
   end if  

   write(2,'(2x,f25.10,5x,f25.10)')  r , temperature(1,j,nz/kn)

end do !j      						   

!------------------
do k=0,nz
   z=k*deltaz 
   
   if (mod(k,kn)==0) write(3,'(2x,f25.10,5x,f25.10)')  z , temperature(1,0,k/kn)

end do !k      						   

!------------------------------------------------ For Phase Equation
do j=0,nr

   if (j<=int(nnrom*nr)) then
      r=j*deltar1
      deltar=deltar1
     else
      r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
      deltar=deltar2
   end if  

   write(5,'(2x,f25.10,5x,2f25.10)') r , phasechange(nt/inn,j,nz)  

end do !j      						   

!------------------ 
do k=0,nz
   z=k*deltaz 
   
   write(6,'(2x,f25.10,5x,2f25.10)') z , phasechange(nt/inn,0,k)      

end do !k      						   

!------------------------------------------------ For feild equations    	
do i=31,nt,(4*inn) 
   
   !if (mod(i,4*inn)==0) then
      
      do j=0,nr
   
         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

         Elec12 = Elec1(i/inn,j,0) * conjg(Elec1(i/inn,j,0)) * 100      
         write( 8,'(2x,f25.20,5x,f25.20)') r , Elec12
   
      end do !j      						   

      write(8,*) !'ZONE I =',i   
   
   !end if !i  

end do !i

!------------------
do i=31,nt,(4*inn) 

   !if (mod(i,4*inn)==0) then     
      
      do k=0,nz
         z=k*deltaz 
         
         if (mod(k,kn)==0) then
         
	    Elec12 = Elec1(i/inn,0,k/kn) * conjg(Elec1(i/inn,0,k/kn)) * 100      
            write( 9,'(2x,f25.20,5x,f25.20)') z , Elec12
   
 	 end if !k  
		    
      end do !k      						   

      write(9,*) !'ZONE I =',i   
   
   !end if !i

end do !i

!==================
do i=31,nt,(4*inn)  
 
   !if (mod(i,4*inn)==0) then

      do j=0,nr

         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

         Elec22 = Elec2(i/inn,j,nz/kn) * conjg(Elec2(i/inn,j,nz/kn)) * 100      
         write(11,'(2x,f25.20,5x,f25.20)') r , Elec22
   
      end do !j      						   

      write(11,*) !'ZONE I =',i       
   
   !end if !i 

end do !i

!------------------
do i=31,nt,(4*inn) 
   
   !if (mod(i,4*inn)==0) then
      
      do k=0,nz
         z=k*deltaz 
         
         if (mod(k,kn)==0) then  
      
	    Elec22 = Elec2(i/inn,0,k/kn) * conjg(Elec2(i/inn,0,k/kn)) * 100      
            write(12,'(2x,f25.20,5x,f25.20)') z , Elec22
         
	 end if !k
      end do !k      						   

      write(12,*) !'ZONE I = ',i    
   
   !end if !i

end do !i

!==================
do i=31,nt,(4*inn) 
   
   !if (mod(i,4*inn)==0) then
   
      do j=0,nr

         if (j<=int(nnrom*nr)) then
            r=j*deltar1
            deltar=deltar1
           else
            r=int(nnrom*nr)*deltar1+(j-int(nnrom*nr))*deltar2
            deltar=deltar2
         end if  

         Elec32 = Elec3(i/inn,j,nz/kn) * conjg(Elec3(i/inn,j,nz/kn)) * 100      
         write(14,'(2x,f25.20,5x,f25.20)') r , Elec32
   
      end do !j      						   

      write(14,*) !'ZONE I =',i   
   
   !end if !i  

end do !i      						   

!================== Transverse profile
!do j=0,nr

   !if (j<=int(4.*nr/5.)) then
      !r=j*deltar1
      !deltar=deltar1
     !else
      !r=5*omegaf+(j-int(4.*nr/5.))*deltar2
      !deltar=deltar2
   !end if  

   !do f=0,360,6
      !fi=f*pi/180.
 
      !Psi32 = Psi3(nt,j,nz) * conjg(Psi3(nt,j,nz))      
      !write(100,'(2x,f25.20,5x,f25.20,5x,f25.20)') r , fi , Psi32
   
  !end do !f

!end do !j      						   

!================== Transverse profile
  
!do j=0,nr

!   if (j<=int(4.*nr/5.)) then
!      r=j*deltar1
!      deltar=deltar1
!     else
!      r=5*omegaf+(j-int(4.*nr/5.))*deltar2
!      deltar=deltar2
!   end if  

   
!   do f=-200,+200
!      x=f*(r/200.)

!      y=sqrt(r**2-x**2)
!      y1 = y  ; y2 = -y

!      Psi32 = Psi3(nt,j,nz) * conjg(Psi3(nt,j,nz))      
!      write(200,'(2x,f25.20,5x,f25.20,5x,f25.20)') x , y1 , Psi32
!      write(200,'(2x,f25.20,5x,f25.20,5x,f25.20)') x , y2 , Psi32
 
!   end do !i
!end do !j

!------------------
do i=31,nt,(4*inn)   
   
   !if (mod(i,(4*inn+31))==0) then
      
      do k=0,nz
         z=k*deltaz 
        
	 if (mod(k,kn)==0) then
      
	    Elec32 = Elec3(i/inn,0,k/kn) * conjg(Elec3(i/inn,0,k/kn)) * 100      
            write(15,'(2x,f25.20,5x,f25.20)') z , Elec32
		 
	 end if !k
		    
      end do !k    

      write(15,*) !'ZONE I = ',i    
   
   !end if !i

end do !i

!**********************************************************************************************************************
!                                      Closing Files and Ending the Program 
!**********************************************************************************************************************

!----------------- for temperature 

close(1)
close(2)
close(3)

!----------------- for phase equation

close(4)
close(5)
close(6)

!----------------- for field equations
close( 7)
close( 8)
close( 9)

close(10)
close(11)
close(12)

close(13)
close(14)
close(15)

!------------------
close(16)
close(17)
close(18)
close(19)
close(20)

!----------------------

write(*,*) 
write(*,*) '---- The results are stored in `.plt` format.                                  &
	         If a different format is required, users can set the desried extension in      &
			   "Determine Filenames & Open files" section of the code or rename the file      & 
			   manually and open it with their preferred software. ----!'	

			
write(*,*) 	
write(*,*) '---- Program Completed ----!'

end program Coupled_G_PW                     

!======================================================================================================================
        
 
