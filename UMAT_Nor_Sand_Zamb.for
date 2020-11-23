    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	!XXXXX This UMAT has been created by Luis Zambrano-Cruzatty based on 2 Nor-Sand   XXXX
	!XXXXX                   UMATs developed by Jorge Macedo 				          XXXX
	!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	
	subroutine UMAT(strss,statev,ddsdde,sse,spd,scd, &
       rpl,ddsddt,drplde,drpldt, &
       stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
       ndir,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
       celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc) 
		!implicit none
      include 'aba_param.inc'
	  
      !implicit real*8(a-h,o-z)
      character*80 cmname
       
      integer, intent (inout):: ntens, ndir, nshr, nstatv, nprops, noel, npt, &
      layer, kspt, kstep, kinc

      double precision, intent(inout):: strss(ntens), statev(nstatv), &
      ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens), &
      stran(ntens), dstran(ntens), time(2), predef(1), dpred(1), &
      props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision, intent(inout):: sse, spd, scd, rpl, drpldt, dtime, temp, &
      dtemp, pnewdt, celent
	  
	  !Parameters:
	  !ELASTIC PARAMETERS:
	  !G: Reference shear modulus
	  !p_ref: reference mean stress
	  !nG: Shear modulus exponent
	  !nu= Poisson ratio
	  !alpha_g: Strain rate factor for shear modulus (Xu&Zhang 2015 and Yamamuro et al 2011)
	  !alpha_K: Strain rate factor for bulk modulus (Xu&Zhang 2015)
	  !RefRate= Reference strain rate (2.77e-5 s-1: 1%/min)
	  
	  !CRITICAL STATE LINE: e_c=\Gamma-\lambda ln(p'_c)
	  !Gamma: Altitude of CSL, defined at 1 kPa
	  !lambda_c: Slope of CSL defined on natural logarithm scale
	  
	  !YIELD FUNCTION PARAMTERS (Jefferies et al 2015)
	  !M_tc: !critical friction ratio for triaxial compression
	  !p_i: Image mean effective stress (state variable)
	  
	  !HARDENING SOFTENING PARAMETERS
	  !CHI: Dilatancy coefficient (Jefferies&Shuttle 2002) [-]
	  !N: Volumetric coupling coefficient
	  !H: Hardening modulus [-]
	  !alpha_p=Strain rate factor for computing p_i (Isotache approach Mesri et al 1995)
	  !alpha_chi= Strain rate factor for dilatancy coefficient
	  
	  !NUMERICAL PARAMTERS
	  !switch_smooth: switch to activate the strain rate smoothening algorithm
	  !N_s: Moving average degree
	  !MAXSUBS: Max sub-stepping steps
	  
	  
	  !STATE VARIABLES AND NEEDED DATA
	  !M_i= Image stress ratio 
	  !p_i: Image mean effective stress [kPa]
	  !eps_rate: Strain rate [-/s]
	  !e_i: current void ratio [-]
	  !N_soft_i: Current number of strain rate values added
	  !SUM_rate: Cumulative strain rates
		
      double precision ::  G_0, p_ref, nG, nu, alpha_g !Elastic variables
      double precision ::  e_o, Gamma, lambda_c !CSL parameters
	  double precision :: M_tc, N !yield function parameters
	  double precision :: CHI_tc, CHIi, H, alpha_CHI, alpha_pi !Hardening parameters
	  double precision :: RefRate !reference strain rate
	  double precision :: bK, bk_0, G
	  double precision :: M_i, p_i, e, psi, SUM_rate, eps_rate, YSFt
	  integer :: MAXSUBS, N_S, N_soft_i, i
      double precision :: DE(6,6)
	  double precision :: p,q,eta, s1_double 
	  double precision, dimension(6) :: dsig, stress, sig, dEpsP, EpsP,Erate, &
										Erate0
	  logical:: switch_smooth, switch_yield

      ! get model parameters
      G_0=		   props(1) !Reference shear modulus [kPa]
      p_ref=	  -props(2) !reference mean effective stress [kPa]
      nG   =	   props(3) !Shear modulus exponent [-]
	  nu   =	   props(4) !Poisson ration [-]

	  
	  !Critical state line parameters
	  e_o=         props(5) !Initial void ratio (defines initial state parameter) [-]
	  Gamma=       props(6) !Altitude of CSL [-]
	  lambda_c=	   props(7) !CSL slope in ln space [-]
	  R=		   props(8) !OCR [-]
	  
	  !Yield surface and hardening parameters
	  M_tc=	       props(9) ! critical friction ratio [-]
	  N=	       props(10) ! volumetric coupling coefficient [-]
	  
	  !Hardening parameters
	  CHI_tc=	   props(11) ! Dilatancy coefficient [-]
	  H=		   props(12) ! Hardening modulus [-]

	 !Strain rate parameters
	  alpha_g=	   props(13) !Strain rate factor for shear modulus [-]
	  alpha_K=	   props(14) !Strain rate factor for bulk modulus [-]
	  alpha_CHI=   props(15) !strain rate factor for the Dilatancy coefficient [-]
	  alpha_pi=	   props(16) !Strain rate factor for the image mean effective stress [-]
	  RefRate=     props(17) !reference strain rate [1/s]

	  !numerical parameters
	  call dbltobool(props(18), switch_smooth)  ! switch for activating strain rate smoothening
	  N_S=           props(19)                  ! Degree of smoothening


	  ! State variables
	  G=             statev(1)               !current shear modulus
	  bK=            statev(2)               !current bulk modulus
	  e=             statev(3)               !current void ratio
	  psi=           statev(4)               !current state parameter
	  CHI_tce=		 statev(5)				 !current dilatancy coefficient
	  p_i =          statev(6)               !image pressure
	  pi_0=			 statev(23)				 !Initial image stress
	  M_i =          statev(7)               !image critical friction ratio
	  call dbltobool(statev(8),switch_yield) ! Point is plasticizing
	  do i=1,6
		  Erate0(i)= statev(8+i)             !previous strain rate
		  Epsp(I)=   statev(14+i)            !current plastic strain (I don't think I need this think about deleting this)
	  end do
	  SUM_rate=		 statev(21)              !Sum of strain rates for smoothing algorithm
	  N_soft_i=      statev(22)              !Current degree of smoothening	

	 X=1/DTIME
	 if (DTIME==0.0d0) then
		ERate= 0.0d0    ! Current strain rate
	 else
		ERate= X*DSTRAN ! Current strain rate
	 end if

	  
	  !Call sub-stepping procedure
	  Call Nor_Sand_Rate(noel, G, g_0, bk, nu, p_ref, nG, e_o, Gamma, lambda_c, M_tc, N, CHI_tc, H, &
		                 alpha_g, alpha_K, alpha_CHI, alpha_pi, RefRate, Erate0, Erate, R, &
		                 switch_smooth, N_S, e, psi, CHI_tce, p_i, pi_0, M_i, SUM_rate, N_soft_i, switch_yield, &
		                 dstran, dtime, strss, Sig, EpsP, DDSDDE)
	  !*
	  !*... stress state parameters update
	  !*
	  Do i=1,ntens
          strss(i) = Sig(i)
	  End Do
	  
	  !Return state variables
	  statev(1)=  G
	  statev(2)=  bK
	  statev(3)=  e	  
	  statev(4)=  psi
	  statev(5)=  CHI_tce
	  statev(6)=  p_i
	  statev(23)= pi_0
	  statev(7)=  M_i
	  statev(8)=  logic2dbl(switch_yield)
	  do i=1,6
		  statev(8+i)= Erate(I)
		  statev(14+i)= Epsp(I)
	  end do
	  statev(21)= Sum_rate
	  statev(22)= N_soft_i
	  !End UMAT routine
	  Return
	end subroutine UMAT
	

	!_________________________________________________________________________________________________________
	!*********************************************************************************************************
	!******************************************* Main subroutine *********************************************
	!*********************************************************************************************************
	!_________________________________________________________________________________________________________
	Subroutine Nor_Sand_Rate(noel, G, G_0, K, nu, p_ref, nG, e_o, Gamma, lambda_c, M_tc, N, CHI_tc, H, &
		                 alpha_g, alpha_K, alpha_CHI, alpha_pi, RefRate, Erate0, Erate, R, &
		                 switch_smooth, N_S, e, psi, CHI_tce, p_i, pi_0, M_i, SUM_rate, N_soft_i, switch_yield, &
		                 dEps, dtime, Sig0, Sig, EpsP, DDSDDE)
	
	!_________________________________________________________________________________
	!Sub-stepping algorithm procedure based in Sloan et al (2001). 
	!Elastic portion is recovered using the a Newton-Raphson procedure
	!_________________________________________________________________________________
	Implicit none
	!Variable declaration
	!input variables
	!Logical
	logical, intent(in):: switch_smooth
	!integers
	integer, intent(in):: noel, N_S
	!double precision
	double precision, intent(in):: G_0, nu, p_ref, nG, e_o, R, Gamma, lambda_c, M_tc, N, CHI_tc, H
	double precision, intent(in):: alpha_g, alpha_K, alpha_CHI, alpha_pi, RefRate, Erate0(6)
	double precision, intent(in):: dEps(6), dtime, Sig0(6)
	!Output variables
	!logical
	logical, intent(inout):: switch_yield
	!integers
	integer, intent(inout):: N_soft_i
	!double precision
	double precision, intent(inout):: G, K, e, psi, CHI_tce, p_i, pi_0, M_i
	double precision, intent(inout):: Erate(6), SUM_rate
	double precision, intent(inout):: Sig(6), EpsP(6), DDSDDE(6,6)
	!____________________________________________________________________________________
	!Local variables
	!Logical
	logical:: ApplyStrainRateUpdates, Locus, IsElasticUnloading, failed
	!integer
	integer:: icount, ITER
	!double precision
	double precision:: SSTOL, FTOL, SPTOL, LTOL, DTmin
	double precision:: IErateI, IErate0I, dErate(6), dErateS(6), DErateSS(6), ErateYield(6), dErate_eff
	doubleprecision:: p,q,eta, p_i0, dM, dp, CHIi, ptrial, qtrial, etatrial
	double precision:: DE(6,6), D1, D2, dSig(6), SigTrial(6)
	double precision:: dEpsvol, dEpsq, dEpsS(6), dEpsSS(6), dEpspS(6), dEpsp(6), Epspv, Epspq, dEpsTrial(6)
	double precision:: alpha, kappa, dummyArg(2)
	double precision:: T, DT, FT, F0
	double precision:: M_i1, M_i2, p_i1, p_i2, e_1, e_2, psi_1, psi_2, CHIi_1, CHIi_2, CHI_tce1, CHI_tce2,&
					   G_1, G_2, K_1, K_2, e_c, dSP1(2), dSP2(2)
	double precision:: dEpsP1(6), dEpsP2(6), Sig1(6), dSig1(6), dSig2(6), EpsP1(6)
	double precision:: RT_Dt, nSigma, qq
	double precision::  Theta, J3, J2, cos3Theta, Mtheta, km !(km is the smoothness parameter)

	  !______________________________________________________________________________________________________
	  !| Initialization	of variables                                                                         |
	  !|_____________________________________________________________________________________________________|
      SSTOL = 0.0001d0 !Tolerance Relative Error (10-3 to 10-5)
      FTOL = 0.000001d0 !Tolerance Error on the Yield surface (10-6 to 10-9)
      SPTOL = 0.0001d0 !Tolerance Softening Parameters (0.0001d0)
	  LTOL= 0.01d0 !Tolerance for elastic unloading
	  ITER=10 ! Number of iterations for stress correction
	  DTmin=0.000001d0
	  km=0.1
	  call getPandQ (Sig0,p,q,eta)
	  if (G==0.0d0) G=G_0*(p/p_ref)**nG
	  if (K==0.0d0) K=2*G*(1+nu)/(3*(1-2*nu))
	  if (e==0.0d0) e=e_o
	  if ((M_i==0.0d0).or.(p_i==0.0d0)) then
		call getP_M_Image(Sig0, M_tc, CHI_tc, R, e, lambda_c, Gamma, N, M_i, p_i, CHIi, psi)
		pi_0=p_i
	  endif
	  call TwoNormTensor(Erate,6,IErateI)
	  call TwoNormTensor(Erate0,6,IErate0I)
	  if (CHI_tce==0.0d0)  CHI_tce=CHI_tc
	  CHIi=CHI_tce/ (1.0d0- lambda_c * CHI_tce/M_tc)
	  !_______________________________________________________________________________________________________
	  
	  !______________________________________________________________________________________________________
	  !| Strain rate smoothening algorithm                                                                   |
	  !|_____________________________________________________________________________________________________|  

	  if (switch_smooth) then
		  N_soft_i=N_soft_i+1
		  if (N_soft_i<N_S) then !Not enough values
			  SUM_rate=SUM_rate+IErateI !Cumulate the strain rate
			  IErateI=SUM_rate/N_soft_i !Takes the average
		  else !Enough values			  
			  SUM_rate=SUM_rate*(1.0-1.0/N_S) !Approximate sum without one term
			  SUM_rate=SUM_rate+IErateI
			  IErateI=SUM_rate/N_S !Averaged strain rate
		  endif
		  call TwoNormTensor(Erate,6,IErate0I)! Norm of the uncorrected strain rate
		  Erate=(IErateI/IErate0I)*Erate !Corrected strain rate tensor
		  call TwoNormTensor(Erate0,6,IErate0I)
	  endif
	  !_______________________________________________________________________________________________________
	  
	  !______________________________________________________________________________________________________
	  !| Update elastic and state parameters due to strain rate                                              |
	  !|_____________________________________________________________________________________________________|	  
	  dErate=Erate-Erate0
	  dErate_eff=IErateI-IErate0I !Increment of effective strain rate
	  ! Degrade G
	  call check4crossing(IErate0I, IErateI, dErate_eff, RefRate, ApplyStrainRateUpdates)	  
	  if (abs(dErate_eff)==IErateI) then!Trim increment of strain rate
		 dErate_eff=dErate_eff-RefRate
	  endif	
	  if (ApplyStrainRateUpdates.and.(.not.switch_yield)) then
		  !Update G and K
		  call UpdateGandKdue2Erate(G_0, p, p_ref, nG, nu, alpha_G, alpha_K, IErateI, &
					     	IErate0I, dErate_eff, RefRate, G, K)
		  !Update M and pi
		  call UpdateMandpidue2Erate(Sig0, Chi_tc, M_tc, N, alpha_pi, alpha_CHI, RefRate, lambda_c, &
									Gamma, e, IErate0I, IErateI, dErate_eff, p_i, pi_0, M_i, &
									psi, CHI_tce, CHIi)

	  endif
	!_________________________________________________________________________________________________________  

	!_________________________________________________________________________________________________________
	!|Elastic procedure																						  |
	!|________________________________________________________________________________________________________|
	  !Assemble Elastic Matrix

        D1  = K+(4*G/3)
        D2  = K-(2*G/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = G 
        DE(5,5) = G
        DE(6,6) = G
	  ! Elastic stress increment
	  call MatVec(DE, 6, dEps, 6, dsig)
	  ! total Elastic stress
	  call AddVec(Sig0, dSig, 1d0, 1d0, 6, SigTrial)  
	  call getPandQ (SigTrial,pTrial,qTrial,etaTrial)
	  call getYieldFunctionNorSand (km,pTrial,qTrial,CHI_tce, CHIi, N, psi, M_tc, p_i,M_i,FT, Locus) !Evaluate yield function on trial stress	   
	  !Check Plasticity
	  if (FT <= FTOL) then !Elastic 
	      !DDSDDE is DE
		  DDSDDE = DE
		  call getDevVolStrain(dEps,dEpsvol,dEpsq)
		  Sig=SigTrial
		  !Update state parameters due to change in Lode's angle
	      e = e + dEpsvol * ( 1. + e) !Update void ratio
          call ElasticUpdating(Sig, dSig, G_0, nu, p_ref, nG, M_tc, p_i, N, &
								IErateI, IErate0I, RefRate, CHIi, e, Gamma, lambda_c, &
		                       psi, M_i,alpha_G, alpha_K, G, K)	
		  
		  if (switch_yield) pi_0=p_i	!In case of unloading	  
		  switch_yield=.false.
		  return
	!********************************************************************************************************
	!________________________________________________________________________________________________________
	!| Plastic strain. Have to integrate model using the elasto-plasticity framework						|
	!| The integration is based on  Sloan et al 2001                                                        |
	!|______________________________________________________________________________________________________|
    !*********************************************************************************************************
		  
	!________________________________________________________________________________________________________
	!| Retrieve alpha																						|
	!|______________________________________________________________________________________________________|
	  
	  else ! Plastic		  

		  call getYieldFunctionNorSand (km,p,q,CHI_tce, CHIi, N, psi, M_tc, p_i,M_i,F0, Locus) !Evaluate yield function on initial stress
		  !Check for elasto plastic transition

		  if (F0 < -FTOL) then !Elasto-plastic transition
			  call getElasticPartNewton(km, FTOL, Sig0, dEps, dErate, Erate0, IErate0I, IErateI, dErate_eff,  &
									refRate, G_0, nu, p_ref, nG, G, K, alpha_G, alpha_K, alpha_chi, alpha_pi, &
									Gamma, lambda_c, CHI_tce, Chi_tc, CHIi, e, psi, M_tc, N, p_i, pi_0, M_i,alpha,&
									SigTrial, dEpsS)
		  else !Check elastic unloading
			  if (F0 <= FTOL) then
				  !Check direction of stress path with respect to actual yield surface
				  call CheckElasticUnloading(LTOL, Sig0, dSig, p_i, M_i, M_tc, CHIi, CHI_tce, &
											N, psi, km, IsElasticUnloading)
				  if (IsElasticUnloading) then ! Must find intersections
					  pi_0=p_i
				  call getElasticPartNewton(km, FTOL, Sig0, dEps, dErate, Erate0, IErate0I, IErateI, dErate_eff,  &
									refRate, G_0, nu, p_ref, nG, G, K, alpha_G, alpha_K, alpha_chi, alpha_pi, &
									Gamma, lambda_c, CHI_tce, Chi_tc, CHIi, e, psi, M_tc, N, p_i, pi_0, M_i,alpha,&
									SigTrial, dEpsS)

				  else !Pure plasticity
					  alpha=0.0d0	
					  dEpsS=dEps
		              SigTrial=Sig0
				  end if
				  
			  else! Invalid stress state
				  ! Assume alpha=0 and solve for pure plasticity (for robustness)
				  alpha=0.0d0	
				  dEpsS=dEps
				  SigTrial=Sig0
			  end if			  
		  end if
	!________________________________________________________________________________________________________
	!_________________________________________________________________________________________________________	  
		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  !!!!!!                                   Sub-stepping Algorithm                             !!!!!!!!!
		  !!!!!!																				     !!!!!!!!!
		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! initiation of sub-stepping
		ErateYield=Erate0
        T = 0.
        dT = 1. !start with the whole step, reduce it if the relative error is larger than tolerance
        iCount = 0
        !sub-stepping
        do while (T < 1.0 )
          iCount = iCount + 1
		  if ( T+ dT > 1. ) dT = 1. - T
          dEpsSS=dT*DEpsS
		  DErateSS= dT* dErate
		  call TwoNormTensor(ErateYield, 6, IErate0I)
		  ErateYield=ErateYield+DErateSS
		  call TwoNormTensor(ErateYield, 6, IErateI)
		  dErate_eff=IErateI-IErate0I !Update strain rate increment 		  
	!________________________________________________________________________________________________________
	!| Compute first increment of stress																	|
	!|______________________________________________________________________________________________________|
		  !Store initial state parameters
		  M_i1=M_i
		  p_i1=p_i
		  e_1=e
		  psi_1=psi
		  CHIi_1=CHIi
		  CHI_tce1=CHI_tce
		  G_1=G
		  K_1=K
		  !First estimate of associated stress
		  call GetdSiganddSP (km, Locus, SigTrial, Epsp, dEpsSS, p_i1, pi_0, M_i1, psi_1, CHI_tce1, M_tc, CHI_tc, CHIi_1&
			         , N, e_o, e_1 ,Gamma, lambda_c, G_0, nG, p_ref, nu, alpha_G, alpha_K, alpha_pi, alpha_chi, &
							  IErate0I, IErateI, dErate_eff, RefRate, H, G_1, K_1, dEpsP1, dSig1, dSP1)

		  !Update for dSig1 and dSP1
		  Epsp1=Epsp+dEpsp1
		  Sig1=SigTrial+dSig1
	!_________________________________________________________________________________________________________
		  
	!________________________________________________________________________________________________________
	!| Compute second increment of stress																	|
	!|______________________________________________________________________________________________________|
		  !Store initial state parameters
		  M_i2=M_i1
		  p_i2=p_i1
		  e_2=e_1
		  psi_2=psi_1
		  CHIi_2=CHIi_1
		  CHI_tce2=CHI_tce1
		  G_2=G
		  K_2=K
		  
		  !Find second approximation evaluated at sig+dSig1 and SP+dSP1
		  call GetdSiganddSP (km, Locus, Sig1, Epsp1, dEpsSS, p_i2, pi_0, M_i2, psi_2, CHI_tce2, M_tc, CHI_tc, CHIi_2&
			         , N, e_o, e_2 ,Gamma, lambda_c, G_0, nG, p_ref, nu, alpha_G, alpha_K, alpha_pi, alpha_chi, &
							  IErate0I, IErateI, dErate_eff, RefRate, H, G_2, K_2, dEpsP2, dSig2, dSP2)
	!________________________________________________________________________________________________________

	!________________________________________________________________________________________________________
	!| Compute second average increment and Error          												    |
	!|______________________________________________________________________________________________________|
		  !Compute average stress and SP
		  Sig=SigTrial+0.5*(dSig1+dSig2)
		  
		  !Determine relative error
		  call TwoNormTensor2((dSig1-dSig2), 6, RT_dT)
		  call TwoNormTensor2(Sig, 6, nSigma)
		  RT_dT=RT_dT/(2.0d0*nSigma)
		  
		  !Check relative error
		  if ((RT_dT>SSTOL)) then !failed > sub-stepping
			  !iCount=iCount+1
			  qq=0.9d0*sqrt(SSTOL/RT_dT)
			  qq=max(qq,0.1)
			  dT=max(qq*dT,dTmin)
			  failed=.true.	
			  ErateYield=ErateYield-DErateSS
		  else !Succeeded
	!_______________________________________________________________________________________________________
			  
	!________________________________________________________________________________________________________
	!| Update parameters and correct stress drift           												|
	!|______________________________________________________________________________________________________|			  
			  !Update state parameters and stress
			  dEpspS=0.5d0*(dEpsp1+dEpsp2)
			  Epsp=Epsp+dEpspS
			  e=e_1
			  p_i=p_i+0.50d0*(dSP1(1)+dSP2(1))
			  e_c=Gamma-lambda_c*log(-p_i)
			  psi=e-e_c
			  CHI_tce=CHI_tce+0.5d0*(dSP1(2)+dSP2(2)) !update CHI_tce
			  CHIi=CHI_tce/(1-lambda_c*CHI_tce/M_tc)
			  call getMlode(Sig, M_tc, Theta, J3, J2, cos3Theta, Mtheta)
			  call GetMiwithPsi(Mtheta, M_tc,CHIi, N, psi,  M_i)
			  G=0.5*(G_1+G_2)
			  K=0.5*(K_1+K_2)
     		  !call UpdateGandKdue2confinement(Sig, 0.5d0*(dSig1+dSig2), G_0, nu, p_ref, nG, IErateI, &
									!			RefRate, alpha_G, alpha_K, G, K) !final update of elastic prop.
			  !test yield function in new stresses
			  call getPandQ(Sig, p,q,eta)
			  call getYieldFunctionNorSand (km,p,q,CHI_tce, CHIi, N, psi, M_tc, p_i,M_i,F0, Locus)
			  
			  if (abs(F0)>FTOL) then !stress back to yield surface
				  call stressCorrection(km, ITER, FTOL, F0, dEpsSS, Sig, G_0, nu, p_ref, nG, p_i, pi_0, M_i, M_tc, N, psi, CHIi, &
			      CHI_tc, CHI_tce, H, e_o, e, Gamma, Lambda_c, IErateI, RefRate, K, G, Locus, Epsp)					  
			  end if
			  qq=min((0.9d0*sqrt(SSTOL/RT_dT)),1.1d0)

			  if (failed) then
				  qq= min(qq , 1.0d0)
				  failed=.false.				  
			  end if
			  
			  T=T+dT
			  dT=qq*dT
			  dT=max(dT,dTmin)
			  dT=min(dT,(1.0d0-T))			  			  
		  end if		  
		end do
		!Assemble Elastic Matrix
        D1  = K+(4*G/3)
        D2  = K-(2*G/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = G 
        DE(5,5) = G
        DE(6,6) = G
		DDSDDE = DE
		switch_yield=.true.
		!Update stresses		
	  end if
	end subroutine Nor_Sand_Rate
                                                  
  
	
	
	
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	!_________________________________________________________________________________________________________
	!*********************************************************************************************************
	!******************************************* Internal subroutines*****************************************
	!*********************************************************************************************************
	!_________________________________________________________________________________________________________
subroutine GetdSiganddSP(km, Locus, Sig, Epsp, dEps, p_i, pi_0, M_i, psi, CHI_tce, M_tc, CHI_tc, CHIi, N, e_o, e &
						,Gamma, lambda_e, G_0, nG, p_ref, nu, alpha_G, alpha_K, alpha_pi, alpha_chi, &
						IErate0I, IErateI, dErate_eff, RefRate, H, G, K, dEpsP, dSig, dSP)
	!_____________________________________________________________________________________
	!Subroutine for correcting the stress to the yield surface
	!_____________________________________________________________________________________
	
	implicit none
	!input variables	
	double precision, dimension(6), intent(in):: dEps, Sig, Epsp
	double precision, intent(in):: G_0, nu, p_ref, nG, alpha_G, alpha_K, km
	double precision, intent(in):: M_tc, N, CHI_tc, H
	double precision, intent(in):: e_o, Gamma, lambda_e
	double precision, intent(in):: IErate0I, IErateI, dErate_eff, RefRate
	double precision, intent(in):: alpha_CHI, alpha_pi
	!output variables	
	double precision, dimension(6), intent(out)::dSig, dEpsP
	double precision, intent(inout):: p_i, M_i, psi, CHI_tce, CHIi, e, G, K, dSP(2), pi_0
	
	logical, intent(inout):: Locus
	!local variables
	double precision, dimension(6,6):: DE
	double precision, dimension(2):: DSPErate	
	double precision, dimension(6):: dFdSig, dSigel, aux, dPPdSig
	double precision, dimension(2):: dFdSP
	double precision, dimension(2,6):: dSPdEpsp, dSPdEpsrate
	double precision:: p, q, eta, psi_i, e_c, p_0, A, Hard, Gi, Ki, dCHI_tce
	double precision:: denom, numerator, lambda, aux2(2), Error, Dummyvar,&
						C1, C2, C3, Chier, Mitcer, Mitc, psir, F1, F2
	double precision:: theta, J3, J2, cos3theta, Mtheta, D1, D2, dEpsvol, dEpsq, dFdM
	logical:: ApplyStrainRateUpdates
	integer:: I, J
	
	!_Get the derivatives___________________________________________________________________________
	
	!step 1 get p and a
	call getPandQ(Sig, p, q, eta)

	!get dFdSig and dPPdSig for inner cap evaluated at Sig, M_i, p_i, psi_i
	call getdFdSig(km, Sig, p_i, M_i, M_tc, CHIi, CHI_tce, N, psi, dFdSig, dPPdSig)

	
	!get dFdSP evaluated at Sig, M_i, p_i, psi_i
	call getdFdSP(km, Sig, M_i, p_i, psi, CHIi, lambda_e, N, CHI_tce, M_tc, p, dFdSP)
	
	! get dSPdEpsp evaluated at Sig, M_i, p_i
	call getdSPdEpsp(Locus, Sig, Epsp, dEps, e_o, H, p_i, pi_0, p, M_i, M_tc,&
					CHIi, CHI_tce, psi, N, lambda_e, dSPdEpsp) 
	!_________________________________________________________________________________________________
	
	!Update parameters due to strain rate  ___________________________________________________________
	Gi=G
	Ki=K
	DSPErate=0.0d0
	dCHI_tce=0.0d0
	!call check4crossing(IErate0I, IErateI, dErate_eff, RefRate,ApplyStrainRateUpdates)
	!  if (ApplyStrainRateUpdates) then	!update for strain rate				
	!	dCHI_tce=Chi_tc*(1.0+alpha_chi*log10(IErateI/refrate))		
	!	DSPErate(2)=dCHI_tce-Chi_tce
	!	call UpdateGandKdue2Erate(G_0, p, p_ref, nG, nu, alpha_G, alpha_K, IErateI, &
	!								IErate0I, dErate_eff, RefRate, Gi, Ki)
	!  end if	  
	!_________________________________________________________________________________________________
	  
	!Ensemble Elastic constitutive matrix _____________________________________________________________
	  
    D1  = Ki+(4*Gi/3)
    D2  = Ki-(2*Gi/3)
    DE  = 0.0
    DE(1:3,1:3) = D2
    DE(1,1) = D1
    DE(2,2) = D1
    DE(3,3) = D1
    DE(4,4) = Gi 
    DE(5,5) = Gi
    DE(6,6) = Gi	  
	!___________________________________________________________________________________________________
	
	!Compute plastic multiplier ________________________________________________________________________
	call MatVec(DE, 6, dEps, 6, dSigel)
	!compute the plastic multiplier lambda
	!fist compute A =dFdSigTDedPdSig (dPdSig=dFdSig for associated flow rule)
	call MatVec(DE, 6, dPPdSig, 6, aux)
	A=0.0d0
	do I=1,6 !dot product
		A=A+dFdSig(I)*aux(I)
	end do
	
	!compute hardening parameter Hard=-dFdSp*dSPdEpsp*dPdSig_________________
	do I= 1,2
		aux2(I)=0.0d0	
	end do
	
	do J= 1,2 !dSPdEpsp*dPdSig_______________________________________________
		do I= 1, 6
			aux2(J)=aux2(J)+dSPdEpsp(J,I)*dPPdSig(I)
		end do		
	end do
	Hard=0.0d0
	do J=1,2
		Hard=Hard+dFdSP(J)*aux2(J)
	end do
	Hard=-Hard
	
	!Compute denominator_______________________________________________________
	denom=A+Hard
	
	!Compute dFdSigT*De*dEps= dFdSigT*dSigel___________________________________
	
	numerator=0.0d0
	do I=1,6
		numerator=numerator+dFdSig(I)*dSigel(I)
	end do
	!Compute strain rate contribution= dFdSP*dSP_Erate _________________________
	if (ApplyStrainRateUpdates) then
		numerator=numerator+dFdSP(1)*DSPErate(1)+dFdSP(2)*DSPErate(2)
	end if
	
	!compute lambda________________________________________________________
	lambda=numerator/denom

	!____________________________________________________________________________________________________
	
	!Compute plastic strain dEpsP=lambda*dPdSig _________________________________________________________
	
	dEpsP=lambda*dPPdSig
	!_____________________________________________________________________________________________________
	
	!compute dSig=dSigel-lambda De*dPPdSig _______________________________________________________________
	
	!De*dFdSig was calculated in aux
	do I=1,6
		aux(I)=lambda*aux(I)
		dSig(I)=dSigel(I)-aux(I) !dSig
	end do
	!______________________________________________________________________________________________________
	
	!Compute dSP=dSPdEpsp*dEpsp+dSPdEpsratedEpsrate*dEpsrate ______________________________________________
	dSP(1)=0.0d0
	dSP(2)=0.0d0
	do J=1,2
		do I=1,6!dSPdEpsp*dEpsp
			dSP(J)=dSP(J)+dSPdEpsp(J,I)*dEpsP(I)
		end do
	end do	
	dSP=dSP+DSPErate
	!______________________________________________________________________________________________________
	
	!Update state parameters ______________________________________________________________________________
	p_i=p_i+dSP(1)
	call getDevVolStrain(dEps, dEpsvol, dEpsq)
	e= e + dEpsvol * ( 1. + e)
	e_c=Gamma-lambda_e*log(-p_i)
	psi=e-e_c
	CHI_tce= CHI_tce+ DSPErate(2)
	dSP(2)=DSPErate(2)
	CHIi=CHI_tce/(1-lambda_e*CHI_tce/M_tc)
	call getMlode(Sig+dSig, M_tc, Theta, J3, J2, cos3Theta, Mtheta)
	call GetMiwithPsi(Mtheta, M_tc,CHIi, N, psi,  M_i)	
	G=Gi
	K=Ki
	!______________________________________________________________________________________________________
	end subroutine GetdSiganddSP
!	
!	
	subroutine stressCorrection(km, MAXIT, FTOL, F0, dEpsp, Sig, G_0, nu, p_ref, nG, p_i, pi_0, M_i, M_tc, N, psi, CHIi, &
			     CHI_tc, CHI_tce, H, e_o, e, Gamma, Lambda_e, IErateI, RefRate, K, G, Locus, Epsp, alpha_G, alpha_K)
	!_____________________________________________________________________________________
	!Subroutine for computing the change in stress (dSig) and state parameters (dSP)
	!_____________________________________________________________________________________
	
	implicit none
	!input variables	
	double precision, dimension(6), intent(in):: dEpsp
	double precision, intent(in):: FTOL
	double precision, intent(in):: G_0, nu, p_ref, nG, M_tc, N, CHIi, CHI_tc, CHI_tce, H, km
	double precision, intent(inout):: p_i, M_i, e, psi
	double precision, intent(in):: e_o, Gamma, lambda_e
	double precision, intent(in):: IErateI, RefRate, alpha_G, alpha_K
	integer, intent(in):: MAXIT
	!output variables	
	double precision, intent(inout):: F0, Sig(6), Epsp(6), G, K, pi_0

	logical, intent(in):: Locus
	!local variables
	double precision, dimension(6)::dEpspS
	double precision, dimension(2)::dSP	
	double precision, dimension(6,6):: DE
	double precision, dimension(6):: dFdSig, dPPdSig, dSigel, aux, zeta, Signew
	double precision, dimension(1,2):: dFdSP
	double precision, dimension(2,6):: dSPdEpsp, dSPdEpsrate
	double precision:: p, q, eta, e_c, A, Hard, dEpsvol, dEpsq
	double precision:: M_in, p_in, en, psin, D1, D2, Gn, Kn, dummy
	double precision:: denom, numerator, lambda, aux2(2), F2
	double precision:: J3, J2, cos3Theta, Mtheta, theta
	integer:: I, J, count
	count=0
	do while ((abs(F0)>FTOL) .and. (count<MAXIT))
		!step 1 get p and a
	!_Get the derivatives___________________________________________________________________________
	!step 1 get p and q
	call getPandQ(Sig, p, q, eta)

	!get dFdSig and dPPdSig for inner cap evaluated at Sig, M_i, p_i, psi_i
	call getdFdSig(km, Sig, p_i, M_i, M_tc, CHIi, CHI_tce, N, psi, dFdSig, dPPdSig)
	
	!get dFdSP evaluated at Sig, M_i, p_i, psi_i
	call getdFdSP(km, Sig, M_i, p_i, psi, CHIi, lambda_e, N, CHI_tce, M_tc, p, dFdSP)
	
	! get dSPdEpsp evaluated at Sig, M_i, p_i
	call getdSPdEpsp(Locus, Sig, Epsp, dEpsp, e_o, H, p_i, pi_0, p, M_i, M_tc, &
					CHIi, CHI_tce, psi, N, lambda_e, dSPdEpsp)
	
	!Ensemble Elastic constitutive matrix _____________________________________________________________
	  
    D1  = K+(4*G/3)
    D2  = K-(2*G/3)
    DE  = 0.0
    DE(1:3,1:3) = D2
    DE(1,1) = D1
    DE(2,2) = D1
    DE(3,3) = D1
    DE(4,4) = G 
    DE(5,5) = G
    DE(6,6) = G	  
	!___________________________________________________________________________________________________
	
	!Compute lambda ____________________________________________________________________________________
	!compute the plastic multiplier lambda
	!fist compute A =dFdSigTDedPdSig (dPdSig=dFdSig for associated flow rule)
	call MatVec(DE, 6, dPPdSig, 6, aux)
	A=0.0d0
	do I=1,6 !dot product
		A=A+dFdSig(I)*aux(I)
	end do
	
	!compute hardening parameter Hard=-dFdSp*dSPdEpsp*dPdSig_________________
	do I= 1,2
		aux2(I)=0.0d0	
	end do
	
	do J= 1,2 !dSPdEpsp*dPdSig_______________________________________________
		do I= 1, 6
			aux2(J)=aux2(J)+dSPdEpsp(J,I)*dPPdSig(I)
		end do		
	end do
	Hard=0.0d0
	do J=1,2
		Hard=Hard+dFdSP(1,J)*aux2(J)
	end do
	Hard=-Hard
	
	!Compute denominator_______________________________________________________
	denom=A+Hard
	
	!Compute dFdSigT*De*dEps= dFdSigT*dSigel___________________________________
	
	numerator=F0
	
	!compute lambda________________________________________________________
	lambda=numerator/denom

		
	!Update Sigma=Sig-lambda*Del*dFdSig
	
	!De*dFdSig was calculated in aux
	do I=1,6
		aux(I)=lambda*aux(I)
		Signew(I)=Sig(I)-aux(I)
	end do

	!Compute plastic strain: dEpspS=lambda* dFdSig
	do I=1, 6
		dEpspS(I)=lambda*dFdSig(I)
	end do
	
	Epsp= Epsp + dEpspS
	
	!Compute dSP=dSPdEpsp*dEpsp+dSPdEpsratedEpsrate*dEpsrate_______________________
	dSP(1)=0.0d0
	dSP(2)=0.0d0
	do J=1,2
		do I=1,6!dSPdEpsp*dEpsp
			dSP(J)=dSP(J)+dSPdEpsp(J,I)*dEpspS(I)
		end do
	end do
	
	!Update SP _____________________________________________________________________________________________
	p_in=p_i+dSP(1)
	call getDevVolStrain(dEpspS, dEpsvol, dEpsq)
	en= e + dEpsvol * ( 1. + e)
	e_c=Gamma-lambda_e*log(-p_i)
	psin=e-e_c
	call getMlode(SigNew,M_tc,theta,J3,J2,cos3Theta,Mtheta)
    call GetMiwithPsi(Mtheta, M_tc,CHIi, N, psi,  M_in)
	Gn=G
	Kn=K
	!call UpdateGandKdue2confinement(Signew, aux, G_0, nu, p_ref, nG, IErateI, RefRate, alpha_G, alpha_K, Gn, Kn)
	!________________________________________________________________________________________________________
	
	!Update F2 ______________________________________________________________________________________________
	call getPandQ(Signew, p,q,eta)
	call getYieldFunctionNorSand(km, p,q, CHI_tce,CHIi, N, psin,M_tc, p_in,M_in, F2, dummy)
	
	if (abs(F2)>abs(F0)) then
		lambda= 0.0d0
		do I=1,6
			lambda=lambda+dFdSig(I)*dFdSig(I)
		end do
		lambda=F0/lambda
		Signew=Sig-lambda*dFdSig
		p_in=p_i
		M_in=M_i
		en=e
		psin=psi
		Gn=G
		Kn=K		
	end if
		p_i=p_in
		M_i=M_in
		e=en
		psi=psin
		G=Gn
		K=Kn
		Sig=Signew
	!Update F0
	call getPandQ(Sig, p,q,eta)
	call getYieldFunctionNorSand(km, p,q,CHI_tce, CHIi, N,psi,M_tc, p_i,M_i,F0,dummy)
	count=count+1
	end do		
	end subroutine stressCorrection
!
!	
!
!
      subroutine getMlode (stress, Mtc, theta, J3, J2, cos3Theta, xM)
	  !_______________________________________________________________________________
	  !Function to obtain a path dependent critical stress ratio
	  !_______________________________________________________________________________
          implicit none
          double precision, intent(in) :: Mtc
          double precision, intent(inout):: stress(6)
          double precision, intent(out) :: theta, J3, J2, cos3Theta, xM
          double precision :: dJdSig(6), dJ3dSig(6)
          double precision ::  J3AJ3,  sin3Theta, pi
          PARAMETER (pi=3.14159265359D0)
          

          call ZERO1(dJdSig,6)
          call ZERO1(dJ3dSig,6)
          xM = 0.
          theta = 0.
          J3 = 0.
          J2 = 0.
          cos3Theta = 0.
          J3AJ3 = 0.
          sin3Theta = 0.
 
          call getInvariants(stress,J2,dJdSig,J3,dJ3dSig)
          if (J2 == 0.) then 
 	            J3AJ3 = 0.
          else
              J3AJ3 = J3/sqrt(J2**3)
          end if

          sin3Theta = 3.*sqrt(3.)/2. * J3AJ3 !+sin definition

          if (sin3Theta > 0.99) sin3Theta = 1.
          if (sin3Theta < -0.99) sin3Theta = -1.
 
          theta = 1./3. * asin(sin3Theta)
          if (theta > 0.523598) theta = 0.523598
          if (theta < -0.523598) theta = -0.523598
		  
		  cos3Theta = cos(3.*theta)
          if (-1.E-8 < cos3Theta < 1.E-8) cos3Theta = 0
          
          ! Jefferies \& Shuttle  2011
          xM = Mtc - Mtc**2./(3. + Mtc) * cos(-3.*theta /2. + pi / 4.) !The minus sign is for p<0		  
  
	end subroutine getMlode

      	  
	!XXXXXX STRESS INVARIANTS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      subroutine getPandQ (stress,p,q,eta)
        implicit none
        double precision, intent(in) :: stress(6)
        double precision, intent(out) :: p, q, eta
        double precision :: sigma(6)
        integer :: I
        sigma = stress
        q = 0.
!c        do while (q < 0.1)
          if (sigma(1) > 0.1)  sigma(1) = 0.1     ! tension cut-off in x-direction
          if (sigma(2) > 0.1)  sigma(2) = 0.1     ! tension cut-off in y-direction
          if (sigma(3) > 0.1)  sigma(3) = 0.1     ! tension cut-off in z-direction
          p  =(sigma(1) + sigma(2) + sigma(3))/3.    ! mean effective stress 
          q  =sqrt(0.5*  ( (sigma(1)-sigma(2))**2. &
            +(sigma(2)-sigma(3))**2. + (sigma(3)-sigma(1))**2. &
            + 6*( sigma(4)**2. + sigma(5)**2. + sigma(6)**2.)))  ! deviatoric stress
!c          if (q < 0.1) sigma(1) = sigma(1) + 0.03
!c        end do
        eta = q/p
!c        stress = -sigma
           
      end subroutine getPandQ
 	   
      subroutine getDevVolStrain (Eps, EpsV, EpsD)
        implicit none
        double precision, intent(in) :: Eps(6)
        double precision, intent(out) :: EpsV, EpsD
        EpsV = Eps(1) + Eps(2) + Eps(3)  ! volumetric strain
        EpsD = sqrt(2./3.)*sqrt( (Eps(1)-EpsV/3.)**2. &
               + (EpS(2)-EpsV/3.)**2. + (EpS(3)-EpsV/3.)**2. &
               + (Eps(4)**2.)/2. + (Eps(5)**2.)/2. + (Eps(6)**2.)/2.) ! deviatoric strain  
      end subroutine getDevVolStrain
 	  
    subroutine getInvariants (stress, J2, dJdSig,J3, dJ3dSig)
		implicit none
		double precision, intent(inout) :: stress(6)
		double precision, intent(out) ::J2,dJdSig(6),J3,dJ3dSig(6)
		double precision :: p, q, S(6), eta
		integer :: I
         
		call getPandQ (stress,p,q,eta)
 
		S(1) = stress(1) - p
		S(2) = stress(2) - p
		S(3) = stress(3) - p
		S(4) = stress(4)
		S(5) = stress(5)
		S(6) = stress(6)                                          
 	  
		J3 = S(1)*S(2)*S(3)-S(1)*S(6)**2-S(2)*S(5)**2 &
			-S(3)*S(4)**2+2*S(4)*S(6)*S(5)
 	 
		J2 = 1./6.*((stress(1)-stress(2))**2 &
			+(stress(2)-stress(3))**2 + (stress(3)-stress(1))**2 ) &
			+ stress(4)**2 + stress(5)**2 + stress(6)**2
			!if (abs(J2) < 0.0001) J2  =  0.0
 	 
		dJdSig(1) = S(1)
		dJdSig(2) = S(2)
		dJdSig(3) = S(3)
		dJdSig(4) = 2.*stress(4) !In the conventional tension as positive the sig here is +
		dJdSig(5) = 2.*stress(5)
		dJdSig(6) = 2.*stress(6)
 
		dJ3dSig(1) = -1./3.*S(1)*S(2) -1./3.*S(1)*S(3) +2./3.*S(2)*S(3) &
					-2./3.*S(6)**2 +1./3.*S(5)**2 +1./3.*S(4)**2
		dJ3dSig(2) = -1./3.*S(1)*S(2) +2./3.*S(1)*S(3) -1./3.*S(2)*S(3) &
					+1./3.*S(6)**2 -2./3.*S(5)**2 +1./3.*S(4)**2
		dJ3dSig(3) = 2./3.*S(1)*S(2) -1./3.*S(1)*S(3) -1./3.*S(2)*S(3) &
					+1./3.*S(6)**2 +1./3.*S(5)**2 -2./3.*S(4)**2
		dJ3dSig(4) = -2.*S(3)*S(4)+2.*S(6)*S(5)
		dJ3dSig(5) = -2.*S(2)*S(5)+2.*S(4)*S(6)
		dJ3dSig(6) = -2.*S(1)*S(6)+2.*S(4)*S(5)
 
    end subroutine getInvariants

	
	subroutine GetMiwithPsi(xM, M_tc, CHI, N, psi_i, M_i)
	!____________________________________________________________
	! Gets the static M_i with changes in the state parameter
	!Jeffereis and Shuttle 2011
	!____________________________________________________________
	implicit none
	double precision, intent(in):: xM, CHI, N, psi_i, M_tc
	double precision, intent(out):: M_i
	M_i=xM*(1.0d0-(N*CHI*abs(psi_i)/M_tc))
	end subroutine GetMiwithPsi
	
	
	subroutine getP_M_Image(Sig, M_tc, CHI_tc,R, e, lambda, Gamma, N, M_i, p_i, CHIi, psi_i)
	   !________________________________________________________________________
	   ! Function to obtain p_i and M_i, psi
	   !________________________________________________________________________
        implicit none
		!input variables   
        double precision, intent(in) :: Sig(6), M_tc, CHI_tc, lambda, Gamma, N, e, R
		!Output variables
		double precision, intent(inout):: M_i, p_i, CHIi, psi_i
		!local variables
		double precision:: p, q, eta, theta, J3, J2, cos3Theta, Mtheta,e_ci

        call getPandQ(Sig, p, q, eta)
		p_i=p/exp(1.0d0)
		p_i=p_i*R
		!Compute Mteta
		call getMlode(Sig, M_tc, theta, J3,J2,cos3Theta,Mtheta) !Correct M due to Lode's angle
		CHIi=CHI_tc/ (1.0d0- lambda* CHI_tc/M_tc) !Get CHI image
		e_ci=Gamma-lambda*log(-p_i) !Get critical voi ratio
		psi_i=e-e_ci	   
		call GetMiwithPsi(Mtheta, M_tc, CHIi, N, psi_i, M_i)
		   
	end subroutine getP_M_Image

  
	
	
subroutine getdFdSig(km, Sig, p_i, M_i, M_tc, CHIi, CHI_tce, N, psi, dFdSig, dPPdSig)
		!PartialFtoPartialSigma
		implicit none
        double precision, intent(in) ::Sig(6), p_i, M_tc, CHIi, N, psi, M_i, CHI_tce, km
		double precision, intent(out) :: dFdSig(6), dPPdSig(6)
        double precision :: dFdp, dpdSig(6), dqdSig(6), dFdM, dThetadSig(6), xM
		double precision :: dMdTheta, dFdq
        double precision :: p, q , eta,J2, J3, dJdSig(6), dJ3dSig(6), &
							theta, cos3Theta, F1(6), F2(6), F3(6), pi, Dmin
		double precision :: p_max, pr, pl, n_L, M_itc, D
		double precision :: dMidMtheta, sigma, iota, C_1, C_2, C_3, C_4, dMidtheta, dFdnL, dnLdTheta
		double precision :: dFdTheta, C_v, C_q, sin_theta
        integer :: I
		  
		PARAMETER (pi=3.14159265359D0)
		call getPandQ (Sig,p,q,eta)

		M_itc=M_tc*(1-CHIi*N*abs(psi)/M_tc)
		p_max=p_i/exp(-CHIi*psi/M_itc)
		pr=p_max*(1+km)
		pl=p_max*(1-km)
		n_L=M_i*(1-CHIi*psi/M_itc)
		
		!_____ Invariant derivatives______________________________________________________________
		call ZERO1(dPdSig,6)
		call ZERO1(dQdSig,6)      
		do I = 1, 3
				dPdSig(I) = 1./3.   
				if (q == 0 ) then
					dQdSig(I)=0
				else
					dQdSig(I) =  3. / 2. / q * (Sig(I) - p)
				end if
		end do
		do I = 4, 6
				dPdSig(I) = 0. 
				if (q == 0 ) then
					dQdSig(I)=0
				else
					dQdSig(I) =  3. / q * (Sig(I) )
				endif
		end do 
		call ZERO1(dJdSig,6)
		call ZERO1(dJ3dSig,6)
		call getInvariants (Sig, J2, dJdSig, J3, dJ3dSig)
		call getMlode (Sig,M_tc,theta,J3,J2,cos3Theta,xM)  
		  
		dMdTheta = -(3./2. * M_tc**2 /(3.+ M_tc))*sin((-3.0d0*theta/2)+pi/4.0d0)
		  
		call ZERO1(dThetadSig,6)
		  
		if (  sqrt(J2) < 0.01d0) then ! avoid dividing by 0
				do I = 1, 6
					dThetadSig(I) = 0.
				end do
		else
				if (-1.e-16 < cos3Theta < 1e-16) then  ! avoid dividing by 0
					if (cos3Theta == 0) then 
						do I = 1, 6
							dThetadSig(I) = 0.
						end do
					else
						do I = 1, 6
						! Need to confirm the correct sign + or - ? 
			dThetadSig(I) = sqrt(3.)/2./ cos3Theta/sqrt(J2**3.)*(dJ3dSig(I) &
						- 3./2.*J3/J2*dJdSig(I) ) 
						end do
					end if
				end if
		end if
		!_______________________________________________________________________________________________
		
		if (p<= pr) then !Use original surface
		!___________ Get dFdM*dMdTheta*dThetadSig=F1__________________	
			dFdM= p * (1.0d0+log(p_i/p))
			dMidMtheta=(1.0d0-(N*CHIi*abs(psi)/M_tc))
			call ZERO1(F1,6)
			do I= 1, 6
				F1(I)=dFdM* dMidMtheta * dMdTheta * dThetadSig(I)
			end do
		!_________ Get dFdp*dp*dSig=F2_________________________________
			dFdP=M_i*log(p_i/p)
			call ZERO1(F2,6)
			do I= 1, 6
				F2(I) = dFdP*dPdSig(I)
			enddo
		!________ Get dFdq*dqdSig= dqdSig= F3__________________________
			call ZERO1(F3,6)
			do I= 1, 6
				F3(I) = dQdSig(I)
			enddo
		!_______ Obtain dFdSig_________________________________________
			dFdSig=F1+F2+F3
!______________________________________________________________________________________________________
!______________________________________________________________________________________________________
		elseif (p<=pl) then ! Use the transition surface
			sigma=1+km
			iota=1-km
			C_1=(log(sigma)-km)/(4*km**3)
			C_2=((log(sigma)+1)/(2*(sigma-iota)))+3*C_1
			C_3=-3*(sigma**2)*C_1+2*sigma*C_2-log(sigma)-1
			C_4=-(sigma**3)*C_1+(sigma**2)*C_2-sigma*(C_3+log(sigma))
		!___________ Get (dFdM*dMdTheta+dFdetaL detaLdtheta) dThetadSig=F1__________________
			dFdM=(C_1*(p**3)/p_max**2)-(C_2*(p**2)/p_max)+(C_3*p)+p_max*C_4
			dMidtheta=(1.0d0-(N*CHIi*abs(psi)/M_tc))*dMdTheta
			dFdnL=p
			dnLdTheta=dMdTheta*(1-(CHIi*(N*abs(psi)+psi)/M_tc))
			dFdTheta=dFdM*dMidtheta+dFdnL*dnLdTheta
			call ZERO1(F1,6)
			do I= 1, 6
				F1(I)=dFdTheta* dThetadSig(I)
			end do
		!_________ Get dFdp*dpdSig=F2_________________________________
			dFdP=(3.0d0*M_i*C_1)*(p/p_max)**2-(2*M_i*C_2)*(p/p_max)+(n_L+M_i*C_3)
			call ZERO1(F2,6)
			do I= 1, 6
				F2(I) = dFdP*dPdSig(I)
			enddo
		!________ Get dFdq*dqdSig= dqdSig= F3__________________________
			call ZERO1(F3,6)
			do I= 1, 6
				F3(I) = dQdSig(I)
			enddo
		!_______ Obtain dFdSig_________________________________________
			dFdSig=F1+F2+F3	
!_______________________________________________________________________________________________________
!_______________________________________________________________________________________________________
		else !Use linear surface
		!___________ Get dFdnL*dnLdtheta*dThetadSig=F1__________________
			dFdnL=p
			dnLdTheta=dMdTheta*(1-(CHIi*(N*abs(psi)+psi)/M_tc))
			call ZERO1(F1,6)
			do I= 1, 6
				F1(I)=dFdnL*dnLdTheta* dThetadSig(I)
			end do
		!_________ Get dFdp*dpdSig=F2____________________________________
			dFdp=n_L
			call ZERO1(F2,6)
			do I= 1, 6
				F2(I) = dFdp*dPdSig(I)
			enddo
		!________ Get dFdq*dqdSig= dqdSig= F3__________________________
			call ZERO1(F3,6)
			do I= 1, 6
				F3(I) = dQdSig(I)
			enddo
		!_______ Obtain dFdSig_________________________________________
			dFdSig=F1+F2+F3				
		endif			
		Dmin=CHIi*psi
		!call ZERO1(dPPdSig,6)
		!dPPdSig=dFdSig
		!dPPdSig=dQdSig- Dmin*dPdSig
		call getDevVolStrain(dFdSig,C_v,C_q)
		sin_theta=(C_v-Dmin*C_q)/((C_v**2 + C_q**2)*(Dmin**2+ 1.0d0))**0.5d0
		dPPdSig=dFdSig
		!if (sin_theta >= 0.0d0) then ! Associated
		!	dPPdSig=dQdSig- Dmin*dPdSig					
		!endif!Non associated

		

 
end subroutine  getdFdSig


subroutine getdFdSP(km, Sig, M_i, p_i, psi, CHIi, lambda, N, CHI_tce, M_tc, p, dFdSP)
	!_____________________________________________________________________
	! obtain the derivative of the yield function with respect to 
	! the state parameters M_i and p_i
	!____________________________________________________________________
	
	implicit none
	!Input variables
	double precision, intent(in):: M_i, p_i, p, psi, CHI_tce, M_tc, CHIi, lambda, N, km
	double precision, dimension(6), intent(in):: Sig		
	!output variables
	double precision, dimension(2), intent(out):: dFdSP
	!local variables
	double precision:: dFdM, Mtheta, theta, J3, J2, cos3Theta,p_max, pr, pl
	double precision:: n_L, sigma, iota, C_1, C_2, C_3, C_4
	double precision:: dFdnL, dFdp_max, dMidpi, dnLdpi, dpmaxdpi, M_itc, dMidChi_tc, dn_LdChi_tc
	double precision:: dpmaxdChi_tc

	M_itc=M_tc*(1-CHIi*N*abs(psi)/M_tc)
	p_max=p_i/exp(-CHIi*psi/M_itc)
	pr=p_max*(1+km)
	pl=p_max*(1-km)
	n_L=M_i*(1-CHIi*psi/M_itc)
	call getMlode(Sig, M_tc, theta, J3, J2, cos3Theta, Mtheta)	
	if (p<= pr) then !Call dF1/dSP
		dFdM=p*(1.0d0+log(p_i/p))			
		dFdSP(1)=-dFdM*Mtheta*N*CHIi*abs(psi)*lambda/(M_tc*psi*p_i)
		dFdSP(1)=dFdSP(1)+(M_i*p/p_i)
		!Here Xhi_tc is set as the second state variable for handling the strain rates
		!dFdChitc=dFdMi*dMidChii*dChiidXhitc			
     	dFdSP(2)=-dfdM*Mtheta*N*abs(psi)/(M_tc*(1-Chi_tce*lambda/M_tc)**2.0d0)

	elseif (p<=pl) then ! dF2/dSp
		sigma=1+km
		iota=1-km
		C_1=(log(sigma)-km)/(4*km**3)
		C_2=((log(sigma)+1)/(2*(sigma-iota)))+3*C_1
		C_3=-3*(sigma**2)*C_1+2*sigma*C_2-log(sigma)-1
		C_4=-(sigma**3)*C_1+(sigma**2)*C_2-sigma*(C_3+log(sigma))
		!______ First term derivatives ______________________________________________
		dFdM=(C_1*(p**3)/p_max**2)-(C_2*(p**2)/p_max)+(C_3*p)+p_max*C_4
		dFdnL=p
		dFdp_max=-2.0d0*M_i*C_1*(p/p_max)**3+ M_i*C_2*(p/p_max)**2+ M_i*C_4
		!____________________________________________________________________________
		!_____Second term derivatives _______________________________________________
		dMidpi=-CHIi*N*Mtheta*lambda*abs(psi)/(M_tc*p_i*psi)
		dnLdpi=-(Mtheta*CHIi*lambda/(M_tc*p_i))*((N*abs(psi)/psi)+ 1.0d0)
		dpmaxdpi=exp(CHii*psi/M_itc)*(1.0d0+((CHIi*lambda/M_itc) *(1.0d0+(CHIi*N*abs(psi)/M_itc))))!FIX THIS. Approximating Mi_tc to M_tc
		!____________________________________________________________________________
		dFdSP(1)=dFdM*dMidpi+ dFdnL*dnLdpi+ dFdp_max*dpmaxdpi
		!____ Now find the dFdChi_tc_________________________________________________
		dMidChi_tc=-Mtheta*N*abs(psi)/(M_tc*(1-Chi_tce*lambda/M_tc)**2)
		dn_LdChi_tc= -Mtheta*(N*abs(psi)+psi)/(M_tc*(1.0d0-Chi_tce*lambda/M_tc)**2)
		dpmaxdChi_tc=p_i*psi*exp(CHii*psi/M_itc)*(1.0d0+CHIi*N*abs(psi)/M_itc)/(M_tc*(1-Chi_tce*lambda/M_tc)**2)
		!____________________________________________________________________________
		dFdSP(2)=dFdM*dMidChi_tc+dFdnL*dn_LdChi_tc+dFdp_max*dpmaxdChi_tc
	else !dF3/dSP
		dFdSP(1)=-p*(Mtheta*CHIi*lambda/(M_tc*p_i))*((N*abs(psi)/psi)+ 1.0d0)
		!____________________________________________________________________________
		dFdSP(2)=-p*Mtheta*(N*abs(psi)+psi)/(M_tc*(1.0d0-Chi_tce*lambda/M_tc)**2)				 
	endif		
end subroutine getdFdSP
	!
	subroutine getdSPdEpsp(Locus, Sig, Epsp, dEps, e_o, H, p_i, pi_0, p, M_i, M_tc, CHIi, &
							CHI_tce, psi, N, lambda, dSPdEpsp)
	!____________________________________________________________________
	! subroutine to get the derivative of the state parameters with respect
	! to the plastic strain
	!____________________________________________________________________
	
	implicit none
	!input variables
	double precision, dimension(6), intent(in):: Epsp, dEps, Sig
	double precision, intent(in):: M_i, p_i, psi, p, e_o, lambda, pi_0
	double precision, intent(in):: M_tc, N, CHIi, Chi_tce, H

	logical, intent(in):: Locus
	!output variables
	double precision, dimension(2,6), intent(out):: dSPdEpsp
	!local variables
	double precision:: theta, J3, J2, cos3Theta, Mtheta, dMitcdpsi
	double precision:: ps, pmax, dpidpsi, dpidMitc, dMidpsi, dpsidpi
	double precision:: dEpsqdEpsp(6), Epspq, Epspv, dpidEpspeq, dpsidEpsv, dpidEpspv
	double precision:: dEpspvdEpsp(6), dMidEpsv, dMidEpsq, M_itc, dEpspv, dEpspq
	integer:: I
	
	call getDevVolStrain(Epsp,Epspv,Epspq)
	!1depspvdepsp
	do I= 1,6
		depspvdepsp(I)=0.0d0
	end do
	depspvdepsp(1)=1.0d0
	depspvdepsp(2)=1.0d0
	depspvdepsp(3)=1.0d0
	
	dpsidEpsv=1.0d0+e_o
	
	!2 get the additional derivative depending on surface_______________________________________
	if (.not.locus) then !outer surface
		!2.1 get dSP/depsp(1)= dpidEpseq*dEpseqdEpsp ____________________________________________
		M_itc=M_tc*(1.0d0-(CHIi*N*abs(psi)/M_tc))
		dpidEpspeq=H*p_i*M_i*((p/p_i)*(p/p_i))*(exp(-CHIi*psi/M_itc)-(p_i/p))/M_itc
		call getdEpseqdEpsp(EpsP,EpsPq,dEpsqdEpsp)
		do I=1,6
			 dSPdEpsp(1,I)=dpidEpspeq*(1.0d0)*dEpsqdEpsp(I)
		end do
		!________________________________________________________________________________________
		!I assume p_i is the only state variable
		do I=1, 6
			dSPdEpsp(2,I)=0.0d0
		end do	
	!______________________________________________________________________________________________
	else !3 inner cap softening
		!3.1 get dSP/depsp(1)= dpidEpseq*dEpseqdEpsp__________________________________________________
		M_itc=M_tc*(1-CHii*N*abs(psi)/M_tc)
		call getDevVolStrain(dEps,dEpspv,dEpspq)
		if (dEpspq>0.0) then
		dpidEpspeq=-0.5d0*H*M_i/M_itc
		else
		dpidEpspeq=0.5d0*H*M_i/M_itc	
		endif
		call getdEpseqdEpsp(EpsP,EpsPq,dEpsqdEpsp)
		do I=1, 6
			dSPdEpsp(1,I)=dpidEpspeq*dEpsqdEpsp(I)
		end do
		!3.2 get dSPdEpsp(2)=dDmindpsi*dpsidEpsv*dEpscdEpsp
		do I=1, 6
			dSPdEpsp(2,I)=CHI_tce*dpsidEpsv*depspvdepsp(I)
		end do		
	endif
	end subroutine getdSPdEpsp
	!
	subroutine getdEpseqdEpsp(EpsP,EpsPEq,DEpsPEqDPS)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the equivalent plastic shear strain
      ! with respect the plastic strain
      !
      !**********************************************************************
 
      implicit none
 
      !Local Variables
      double precision :: k1, k2, k3
      double precision :: EpsPM
      double precision, dimension(3) :: EpsDev
      !In Variables
      double precision, intent(in), dimension(6) :: EpsP
      double precision, intent(in) :: EpsPEq
      !Out Variables
      double precision, intent(out), dimension(6):: DEpsPEqDPS
 
      k1 = 2.0d0/(3.0d0*EpsPEq)
      if (EpsPEq < 0.00000000001d0) then
        k1 = 0.0d0
      end if
      k2 = k1 * 1.0d0/3.0d0
      k3 = k1 * 2.0d0
 
      EpsPM = k2 * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM
 
      DEpsPEqDPS(1) = k2 * ( 2.0d0*EpsDev(1) - EpsDev(2) - EpsDev(3))
      DEpsPEqDPS(2) = k2 * (-EpsDev(1) + 2.0d0*EpsDev(2) - EpsDev(3))
      DEpsPEqDPS(3) = k2 * (-EpsDev(1) - EpsDev(2) + 2.0d0*EpsDev(3))
      DEpsPEqDPS(4) = k2 * EpsP(4)
      DEpsPEqDPS(5) = k2 * EpsP(5)
      DEpsPEqDPS(6) = k2 * EpsP(6)
 
	end subroutine getdEpseqdEpsp

	
	Subroutine UpdateMandpidue2Erate(Sig, Chi_tc, M_tc, N, alpha_pi, alpha_CHI, RefRate, lambda &
									,Gamma, e, IErate0I, IErateI, dErate_eff, p_i, pi_0, M_i, &
									 psi,CHI_tce, CHIi)
	implicit none
	!input variables	
	double precision, intent(in):: Sig(6), Chi_tc, M_tc, N
	double precision, intent(in):: alpha_pi, alpha_CHI, RefRate
	double precision, intent(in):: lambda, Gamma, e
	double precision, intent(in):: IErateI, IErate0I, dErate_eff	
	!out variables
	double precision, intent(inout):: p_i, M_i, psi, CHIi, CHI_tce, pi_0
	!local variables
	double precision:: e_ci, Aux, theta, J3, J2, cos3Theta, Mtheta, dp, Erate0, Epspq
	double precision:: dCHI_tce, CHIiold, dpidEpsqp0, dpidEpsqpf, dCHIi, Nsum, &
						C1, C2, C3, Mi_tc, ps, pmax, dpi_0, dpidpio, dp_i, dgdpi
	!Update Chi	

	CHI_tce=CHI_tc*(1.0d0+alpha_CHI*log10(IErateI/RefRate))
	p_i=pi_0*(1.0d0+alpha_pi*log10(IErateI/RefRate))
	CHIi=CHI_tce/(1- lambda*CHI_tce/M_tc)
	e_ci=Gamma-lambda*log(-p_i)
	psi=e-e_ci
	call getMlode(Sig, M_tc, theta, J3,J2,cos3Theta,Mtheta)
	call GetMiwithPsi(Mtheta, M_tc, CHIi, N, psi, M_i)
	end subroutine UpdateMandpidue2Erate
	
	
	
	subroutine UpdateGandKdue2Erate(G_0, p, p_ref, nG, nu, alpha_G, alpha_K, IErateI, &
									IErate0I, dErate_eff, RefRate, G, K)
	implicit none
	!input variables
	double precision, intent(in):: G_0, p, p_ref, nG, nu
	double precision, intent(in):: alpha_G, alpha_K, IErateI, IErate0I, dErate_eff, RefRate
	!output variables
	double precision, intent(inout):: G, K
	!local variables
	double precision:: dG, K_0
	
	!Update G
	G=G_0*((p/p_ref)**nG)*(1.0d0+alpha_K*log10(IErateI/RefRate))
	K=G_0*((p/p_ref)**nG)*2*(1+nu)/(3*(1-2*nu))*(1.0d0+alpha_G*log10(IErateI/RefRate))

	end subroutine UpdateGandKdue2Erate
	
subroutine UpdatePardue2StrainRate(alpha,IErate0I, IErateI, dNErate, Refrate, Par_0, Par)
	!***************************************************************************
	! Updates parameters using the logarithmic law
	!
	!***************************************************************************
	implicit none
	!input variables
	double precision, intent(in):: alpha, Refrate, dNErate, Par_0
	!inout variables
	double precision, intent(inout):: Par,  IErateI, IErate0I
	!local variables
	integer :: i, n, usecase
	double precision:: eratio, dPar, K
	logical:: usecase1
	
	if (dNErate> 0.0d0) then !Increasing strain rate
		eratio=dNerate/IErate0I
		if ((IErate0I<=RefRate).or.eratio>1.0d0) then !Singularity use case 2
			if (IErate0I<=RefRate) then			
				usecase=0
				IErate0I=Refrate
			else
				usecase=2
				eratio=dNerate/IErateI
			endif		
		else !use case 1
			usecase=1 !backward approximation
		endif		
	else !decreasing strain rate
		eratio=dNerate/IErateI
		if ((IErateI<=RefRate).or.eratio<-1.0d0)  then !Singularity use case 1
			if (IErateI<=RefRate) then
				usecase=0
				IErateI=Refrate
			else
				usecase=1
				eratio=dNerate/IErate0I
			endif			
		else
			usecase=2
		endif		
	endif

	if (abs(eratio)<0.2d0) then !Small ratio small error on approximations (See figure XXXX)
		n=3
	else if (abs(eratio)<0.8) then
		n=30
	else if (abs(eratio)<0.9) then
		n=100
	else !Big ratio creates a big error on the approximation of the parameter
		n=300
	end if

	dPar=0.0d0
	
	if (usecase==0) then !Use equation to approximate increment
		dPar=Par_0*alpha*log10(IErateI/Refrate)
		Par=Par+dPar
		return
	else if (usecase==1) then
		K=-1.0d0
	else
		K=1.0d0
	endif	
	do i= 1, n		
		dPar=dPar+(K**dble(i))*(eratio**dble(i))/dble(i)
	end do
	dPar=dPar*alpha*Par_0*0.434d0 ! See equation XXXX
	Par=Par+dPar
	end subroutine UpdatePardue2StrainRate	

    
    subroutine getYieldFunctionNorSand (km,p,q,CHI_tc, CHIi, N, psi, M_tc, p_image,M_image,yield, Locus)
	   !_______________________________________________________________
	   ! Yield function or plastic potential surface for Nor-Sand
	   !_______________________________________________________________
           implicit none
           double precision, intent(in) :: km, p, q, M_image, p_image, CHI_tc, psi, M_tc, CHIi, N
           double precision, intent(out) :: yield
		   logical, intent(out):: Locus
		   double precision:: p_max, pr, pl, n_L, sigma, iota, C_1, C_2, C_3, C_4, M_itc
		   Locus=.false.

		   M_itc=M_tc*(1-CHIi*N*abs(psi)/M_tc)
		   p_max=p_image/exp(-CHIi*psi/M_itc)
		   pr=p_max*(1+km)
		   pl=p_max*(1-km)
		   n_L=M_image*(1-CHIi*psi/M_itc)
		   
		   if (p<= pr) then !Call F1 (Normal Norsand yield surface)
			   yield  =  q + p*M_image *(1.0d0 + log(p_image/p))
		   elseif (p<=pl) then ! in between pr and pl evaluate cubic spline
			      sigma=1+km
				  iota=1-km
				  C_1=(log(sigma)-km)/(4*km**3)
				  C_2=((log(sigma)+1)/(2*(sigma-iota)))+3*C_1
				  C_3=-3*(sigma**2)*C_1+2*sigma*C_2-log(sigma)-1
				  C_4=-(sigma**3)*C_1+(sigma**2)*C_2-sigma*(C_3+log(sigma))
				  yield=q+(M_image*C_1*p**3/(p_max**2))-(M_image*C_2*p**2/p_max)+((n_L+M_image*C_3)*p)+(p_max*M_image*C_4) 			   
		   else !Linear surface
				  yield=q+n_l*p
		   endif
		   
		!! Evaluate main yield function
		!	yield1  =  q + p*M_image *(1.0d0 + log(p_image/p))
		!! Evaluate inner cap
		!	yield2  =  -p_image+p*exp(-CHI_tc*psi/M_tc)            
		!! Set yield to maximum value
		!	!yield= max(yield1, yield2)
		!	yield= yield1
		!! Activate inner cap?
		!	!Locus=((yield2>0.0d0).and.(yield2>yield1)) 
		!	!if ((abs(yield1)<=0.0001).and.(abs(yield2)<=0.0001)) Locus=.false.
	   end subroutine getYieldFunctionNorSand


 
      subroutine getElasticPartPegasus(km, F0, F1, alpha0, alpha1, FTOL, stress, dSig, &
								p_i, M_i, M_tc, CHI_tc, CHIi, N, psi, dstran, dEpsS, alpha, Sig)
	  !_________________________________________________________________
	  !Pegasus algorithm to determine elastic part in elastic unloading
	  !_________________________________________________________________
        implicit none
		!Input variable          
        double precision,intent(in):: p_i, M_i, CHI_tc, psi, M_tc, FTOL, km, CHIi, N
		!output variables
        double precision,intent(inout)::dSig(6), alpha0, alpha1, stress(6)
        double precision,intent(inout)::F0, F1, alpha, dEpsS(6), dstran(6), Sig(6)
		!local variables
        double precision :: FNew,SigNew(6),pNew,qNew,etaNew     
        integer :: check, iteration, I, MAXITS
		logical:: Locus

        check = 0
        !get elastic contribution (Pegasus method)
        FNew= 1000
		MAXITS=10
        iteration = 0
        do while ( abs(FNew) > FTOL ) 
          alpha = alpha1 - F1*(alpha1-alpha0)/(F1-F0)
          do I = 1,6
            SigNew(I) = stress(I) + alpha*dSig(I)
          end do
          call getPandQ (SigNew, pNew, qNew, etaNew)
          call getYieldFunctionNorSand (km,pNew,qNew, CHI_tc,CHIi, N, psi, M_tc, p_i, M_i, FNew, Locus)
          if ( (FNew*F0) < 0. ) then
            alpha1 = alpha
            F1 = Fnew
          else
            F1 = F1*F0/(F0+FNew)
			alpha0=alpha
			F0=FNew
		  end if
		  
          iteration = iteration + 1
          if (iteration > MAXITS) then !control maximum iteration
            FNew = 0.
            alpha = 0.
            check = 999 
          else
            check = 1
          end if
        end do 
              
        if (alpha /= alpha) then
            alpha = 0.
            check = -1
        elseif (alpha > 1.) then
              alpha = 0.
              check = -1
        elseif (alpha < 0.) then 
              alpha = 0.
              check = -1
        end if
        if (check == 1.) then ! update stress for elastic part
            dEpsS = (1-alpha) * dstran           
            Sig= stress + alpha*dSig
        end if
      end subroutine getElasticPartPegasus
!
      subroutine getElasticPartUnload(km, Locus, F0, stress, dSig, Sig, dstran, &
						CHI_tc, CHIi, N, psi, M_tc,M_i, p_i, FTOL, dEpsS, alpha)
	  !__________________________________________________________________
	  !Elastic unloading Pegasus scheme
	  !__________________________________________________________________
        implicit none
		!input variables		 		
		double precision, intent(in) :: M_i, p_i, CHI_tc, CHIi, N, psi, M_tc, km
        double precision,intent(in):: FTOL
		logical, intent(in):: Locus
        double precision,intent(inout)::stress(6), dEpsS(6), dstran(6)
        double precision,intent(inout):: F0, alpha, Sig(6), dSig(6)
        double precision :: FNew,Sig1(6),pNew,qNew,eta,  dEps_e(6)      
        double precision :: alpha0, alpha1, F1, Fsave
        double precision :: dalpha
        integer :: i,j, Is, NSUB, MAXITS
 
        NSUB = 10
        MAXITS = 3
        alpha0 = 0.
        alpha1 = 1.
        Fsave = F0
        do i =1, MAXITS
          dalpha = (alpha1 - alpha0)/NSUB
          do j = 1, NSUB    
            alpha = alpha0 + dalpha
			do Is = 1,6
              Sig1(Is) = stress(Is) + alpha*dSig(Is)
			  dEps_e(Is)= alpha*dstran(Is)
			end do
			!get invariant stresses of stress trial
			call getPandQ (Sig1,pNew,qNew,eta)           
            !test yield function
			call getYieldFunctionNorSand (km,pNew,qNew, CHI_tc,CHIi, N, psi, M_tc,p_i,M_i,FNew, Locus) 
			
            if ( FNew > FTOL ) then
              alpha1 = alpha
              if ( F0 < -FTOL ) then
                F1 = FNew
                go to 200
              else
                alpha0 = 0.0d0
                F0 = Fsave
                EXIT
              end if
            else
              alpha0 =  alpha
              F0 = FNew
            end if
          end do
        end do
        !did not found the intersection with yield surface
        alpha=0.0d0
		dEpsS=dstran
		Sig=stress
		
		return
        !found the intersection, call Pegasus scheme to determine the exact value
200		call getElasticPartPegasus(km,F0, F1, alpha0, alpha1, FTOL, stress, dSig, &
								p_i, M_i, M_tc, CHI_tc,CHIi, N, psi, dstran, dEpsS, alpha, Sig)

      end subroutine getElasticPartUnload
!
!	
	subroutine getElasticPartNewton(km, FTOL, Sig_0, dEps, dErate, Erate0, IErate0I, IErateI, dErate_eff,  &
									refRate, G_0, nu, p_ref, nG, G, K, alpha_G, alpha_K, alpha_chi, alpha_pi, &
									Gamma, lambda, CHI_tce, Chi_tc, CHIi, e, psi, M_tc, N, p_i, pi_0, M_i,alphaNewton,&
									Sig, dEpsS)
	!__________________________________________________________________
	!Routine to obtain elastic part of loading path using the
	!The Newton Raphson method
	!__________________________________________________________________
	implicit none
	!input variables
	double precision, intent(in):: FTOL, Sig_0(6), dEps(6), RefRate, km
	double precision, intent(in):: Chi_tc, M_tc, N, pi_0, alpha_chi, alpha_pi
	double precision, intent(in):: G_0, nu, p_ref, nG, alpha_G, alpha_K
	double precision, intent(in):: Gamma, Lambda, IErateI
	!out variables
	double precision, intent(inout):: alphaNewton, e, p_i, chi_tce, psi, Chii, G, K, M_i
	double precision, intent(out):: Sig(6), dEpsS(6)
	double precision, intent(inout)::dErate(6), Erate0(6), IErate0I, dErate_eff
	!local variables
	double precision:: FNewton, e_new, p_inew, chi_tcenew, psi_new, e_c, Chii_new, M_inew, Fp
	double precision:: dEpsT(6), dErateT(6), ErateT(6), IErateTI, dSigEl(6), SigT(6)
	double precision:: p, q, eta, Gi, Ki, dEpsVol, dEpsq, Mtheta, theta, J3, J2, cos3Theta, DeltaK, deltaG
	double precision:: dSPRate(2), D1, D2, DE(6,6), dDdK(6,6), dDdG(6,6), AlphaD(6,6), dSigdalpha(6), &
						dFdSig(6), dFdSP(2)
	integer:: I
	logical:: ApplyStrainRateUpdating, Locus
	
	Locus=.false.
	
	!Initialize the data	
	alphaNewton=1.0d0
	FNewton=999.0
	do while ((abs(FNewton)>FTOL).or.(FNewton<0.0))
		e_new=e
		dEpsT=alphaNewton*dEps !Guess of alpha
		derateT=alphaNewton*dErate
		ErateT=Erate0+derateT
		call TwoNormTensor(ErateT,6, IErateTI) !Effective strain rate		
		call getDevVolStrain(dEpsT, dEpsVol, dEpsq)
		dErate_eff=IErateTI-IErate0I
		call check4crossing(IErate0I, IErateTI, dErate_eff, RefRate, ApplyStrainRateUpdating)
		
		!__________________________________________________________________________________________
		call getPandQ(Sig_0, p, q, eta)
		if (ApplyStrainRateUpdating) then !update state parameters due to strain rate
		p_inew=pi_0*(1.0+alpha_pi*log10(IErateTI/refRate))
		chi_tcenew=CHI_tc*(1.0+alpha_chi*log10(IErateTI/refRate))
		Gi=(G_0*(p/p_ref)**nG)*(1.0+alpha_G*log10(IErateTI/refRate))
		Ki=(G_0*(p/p_ref)**nG)*(2.*(1.+nu)/(3.*(1.-2.*nu)))*(1.0+alpha_K*log10(IErateTI/refRate))		
		else !Check This
		Gi=G
		Ki=K		
		p_inew=p_i
		Chi_tcenew=CHi_tce
		end if
		
		e_new=e_new+ dEpsVol*(1.+e_new) !Updated void ratio
		e_c=Gamma-lambda*log(-p_inew)
		psi_new=e_new-e_c
		Chii_new=chi_tcenew/(1.0-lambda*chi_tcenew/M_tc)
		call getMlode(Sig_0, M_tc,theta, J3, J2, cos3Theta, Mtheta)
		call GetMiwithPsi(Mtheta, M_tc, Chii_new, N, psi_new, M_inew)
		! Get the increment of state parameters due to strain rate________________________________
		dSPRate(1)=p_inew-p_i
		dSPRate(2)=Chi_tcenew-CHi_tce		
		!__________________________________________________________________________________________
		!Get the increment of stress
		D1  = Ki+(4*Gi/3)
        D2  = Ki-(2*Gi/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = Gi 
        DE(5,5) = Gi
        DE(6,6) = Gi
		call MatVec(DE, 6, dEpsT, 6, dSigEl)
		SigT=Sig_0+dSigEl
		!___________________________________________________________________________________________
		!Get the derivative of F with respect to the state parameters
		call getdFdSP(km, SigT, M_inew, p_inew, psi_new, CHIi_new, lambda, N,&
						CHI_tcenew, M_tc, p, dFdSP)
		Fp=(dFdSP(1)*dSPRate(1)+dFdSP(2)*dSPRate(2))/alphaNewton
		!___________________________________________________________________________________________
		!Get dFdalpha*dSigdAlpha
		dDdK=0.0
		dDdK(1:3,1:3) = 1.0
		dDdG=0.0
		dDdG(1:3,1:3)=-2./3.
		dDdG(1,1)=4./3.
		dDdG(2,2)=4./3.
		dDdG(3,3)=4./3.
		dDdG(4,4)=1.0
		dDdG(5,5)=1.0
		dDdG(5,5)=1.0
		DeltaK=Ki-K
		DeltaG=Gi-G
		AlphaD=alphaNewton*(DeltaK*dDdK+DeltaG*dDdG)
		AlphaD=DE+AlphaD
		call MatVec(AlphaD, 6, dEps, 6, dSigdAlpha)
		!____________________________________________________________________________________________
		!Get derivative with respect to the stress tensor
		call getdFdSig(km, SigT, p_inew, M_inew, M_tc, CHIi_new, CHI_tcenew, N, psi_new, &
						dFdSig, dFdSig)
		do I=1,6 !dot product
			Fp=Fp+dFdSig(I)*dSigdAlpha(I)
		enddo
		!_____________________________________________________________________________________________
		!Get new Alpha
		call getPandQ(SigT, p, q, eta)
		call getYieldFunctionNorSand(km, p, q, CHI_tcenew, Chii_new, N, psi_new, M_tc, p_inew,M_inew,FNewton, Locus)
		AlphaNewton=AlphaNewton-FNewton/Fp		
	end do
	!output
	dEpsS=(1.0-alphaNewton)*dEps
	dErate=(1.0-alphaNewton)*dErate
	Erate0=ErateT
	IErate0I=IErateTI
	dErate_eff=IErateI-IErate0I
	e=e_new
	p_i=p_inew
	psi=psi_new
	Chi_tce=CHI_tcenew
	Chii=Chii_new
	M_i=M_inew
	G=Gi
	K=Ki
	Sig=SigT
	end subroutine getElasticPartNewton
	
	 Subroutine CheckElasticUnloading(LTOL, Sig, dSig, p_i, M_i, M_tc, CHIi, CHI_tce,&
									  N, psi, km, IsElasticUnloading)
	!_____________________________________________________________________	
	! It returs true if the stress path correspond to elastic unloading 	
	!_____________________________________________________________________

	implicit none

	!Local Variables
	double precision ::  NormDSig, NormDFDSig, DSigdpDFDSig, CTeta, dummy(6)
	double precision, dimension(6) :: dFdSig
	!In Variables
	double precision, intent(in), dimension(6) :: Sig, dSig
	double precision, intent(in) :: M_tc, N, p_i, psi, CHIi, CHI_tce, LTOL, M_i, km
	
	!Out Variables
	Logical, intent(out) :: IsElasticUnloading
	
	
	!Calulate yield function derivative with respect to stress
	call getdFdSig(km, Sig, p_i, M_i, M_tc, CHIi, CHI_tce, N, psi, dFdSig, dummy)

	!calculate the 2-norms of the tensors
	call TwoNormTensor2(dSig, 6, NormDSig)
	call TwoNormTensor2(dFdSig, 6, NormDFDSig)
	!Calculate Dsig:DFDSig
	call MatAdpMatB(DSig, DFDSig, 6, DSigdpDFDSig)
	  
	!Cos of angle between stress path and normal vector to yield surface
	CTeta=DSigdpDFDSig/(NormDSig*NormDFDSig)
	  
	  If (CTeta < -LTOL) then !Elastic Unloading
		  IsElasticUnloading=.true.
	  else 
		  IsElasticUnloading=.false.
	  end if
	end subroutine CheckElasticUnloading
	
	subroutine UpdateGandKdue2confinement(Sig, dSig, G_0, nu, p_ref, nG, IErateI, RefRate, alpha_G, alpha_K, G, K)
	implicit none
	!In variables
	double precision, intent(in):: Sig(6), dSig(6), G_0, nu, p_ref, nG
	double precision, intent(in):: IErateI, RefRate, alpha_g, alpha_K
	!out variables
	double precision, intent(inout):: G, K
	!local variables
	double precision:: p, q, eta, dp, DeltaG, FacG, FacK
	call getPandQ (Sig, p, q, eta)
	DeltaG=G_0*nG*((p/p_ref)**(nG-1))/p_ref
	if (IErateI>RefRate) then
		FacG=1+alpha_G*log10(IErateI/RefRate)
		FacK=1+alpha_K*log10(IErateI/RefRate)
	else
		FacG=1.
		FacK=1.
	endif
	call getPandQ (dSig, dp, q, eta)
	G=G+DeltaG*FacG*dp
	FacK=FacK*2*(1+nu)/(3*(1-2*nu))
	K=K+DeltaG*FacK*dp
	end subroutine UpdateGandKdue2confinement
	
	subroutine ElasticUpdating(Sig, dSig, G_0, nu, p_ref, nG, M_tc, p_i, N, &
								IErateI, IErate0I, RefRate, CHIi, e, Gamma, lambda &
		                        ,psi, M_i, alpha_G, alpha_K,G, K)
	!______________________________________________________________________________
	!Updates the state parameter for an elastic load path
	!______________________________________________________________________________
	implicit none
	!Input variables
	double precision, intent(in):: Sig(6), dSig(6), M_tc, p_i, N, CHIi
	double precision, intent(in):: G_0, nu, p_ref, nG, IErateI,IErate0I,  RefRate
	double precision, intent(in):: e, Gamma, lambda, alpha_G, alpha_K
	!out variables
	double precision, intent(inout):: psi, M_i, G, K
	!local variables
	double precision:: theta, J3, J2, cos3Theta, Mtheta, e_c, p, q, eta
	logical :: ApplyStrainRateUpdating
	
	call getMlode (Sig, M_tc,theta,J3,J2,cos3Theta,Mtheta)
	call getPandQ(Sig, p,q,eta)
	e_c=Gamma-lambda*log(-p_i)
	psi=e-e_c
	call GetMiwithPsi(Mtheta, M_tc, CHIi, N, psi, M_i)
	!Update Only for strain rate reference 
	call check4crossing(IErate0I, IErateI, (IErateI-IErate0I), RefRate, ApplyStrainRateUpdating)
	if (ApplyStrainRateUpdating) then
	G=G_0*((p/p_ref)**nG)*(1.0d0+alpha_G*log10(IErateI/RefRate))
	K=G_0*((p/p_ref)**nG)*2.0*(1.0+nu)/(3.0*(1.0-2.0*nu))*(1.0d0+alpha_K*log10(IErateI/RefRate))
	else
	G=G_0*((p/p_ref)**nG)
	K=G_0*((p/p_ref)**nG)*2.0*(1.0+nu)/(3.0*(1.0-2.0*nu))
	end if
	end subroutine ElasticUpdating

!/////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	!_________________________________________________________________________________________________________
	!*********************************************************************************************************
	!******************************************* Math subroutines*********************************************
	!*********************************************************************************************************
	!_________________________________________________________________________________________________________
      subroutine ZERO2(A,N)
         implicit none
         integer, intent(in) :: N
         double precision, intent(inout) :: A (N,N)
         integer :: I,J
 	     do I = 1,N
 	         do J=1,N
 	             A(I,J) = 0.
 	         end do
 	     end do
      end subroutine ZERO2
 
      subroutine ZERO1(A,N)
          implicit none
          integer, intent(in) :: N
          double precision, intent(inout) :: A(N)
          integer :: I
 	      do I = 1, N
 	        A(I) = 0.
 	      end do
      end subroutine ZERO1
 	          
      subroutine addMatComp (A, B, C)
           implicit none
           double precision, intent(in) :: A(6,6), B(6,6)
           double precision, intent(out) :: C (6,6)
           integer :: I , J
 	      do I = 1, 6
              do J = 1, 6
                      C(I,J) = A(I,J) + B(I,J)
              end do
           end do
      end subroutine addMatComp

      subroutine AddVec (A, B, coef1, coef2, N, C)
          implicit none
          double precision, intent(in) :: A(N), B(N)
          double precision, intent(in) :: coef1, coef2
          double precision, intent(out) :: C(N)
          integer :: I, J, N

          call ZERO1(C,N)
          do I = 1, N
             C(I) = coef1*A(I) + coef2*B(I)
          end do
      end subroutine AddVec
      
      subroutine MatVec (A, M, B, N, C)
          implicit none
          double precision, intent(in) :: A(M,N), B(N)
          double precision, intent(out) :: C (N)
          integer :: I, J, M, N

          call ZERO1(C,N)
          do I = 1, M
              do J = 1, N
              	C(I) = C(I) + A(I,J) * B(J)
              end do
          end do
	end subroutine MatVec
	
Subroutine TwoNormTensor(Tensor, N, TwoNorm)
!***********************************************************************
!
!     Calculate 2NormTensor = sqrt(A:A)
!
! I   Tensor  : (Square or vecetor of dimension N)
! I   N     :   Number of elments
! O   2Norm : Resulting norm
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0	  
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+(0.5d0*Tensor(I)*Tensor(I))
	  end do
	  TwoNorm=sqrt(TwoNorm)
	  TwoNorm=abs(TwoNorm)
	
	end subroutine TwoNormTensor

Subroutine TwoNormTensor2(Tensor, N, TwoNorm)
!***********************************************************************
!
!     Calculate 2NormTensor = sqrt(A:A)
!
! I   Tensor  : (Square or vecetor of dimension N)
! I   N     :   Number of elments
! O   2Norm : Resulting norm
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0	  
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+(2.0d0*Tensor(I)*Tensor(I))
	  end do
	  TwoNorm=sqrt(TwoNorm)
	  TwoNorm=abs(TwoNorm)
	
	end subroutine TwoNormTensor2
	Subroutine MatAdpMatB(xMatA, xMatB,N, DProduct)
!***********************************************************************
!
!     Calculate the scalar produc of (A:B)
!
! I   xMatA xMATB  : (Tensors N)
! I   N     :   Number of elments
! O   DProduct : Resulting scalar produc
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension xMatA(N), xMatB(N)
!***********************************************************************
	  X=N/2
      DProduct=0.0d0	  
	  Do I=1,X
		  DProduct=DProduct+xMatA(I)*xMatB(I)
	  end Do
	  Do I=X+1,N
		  DProduct=DProduct+2.0d0*xMatA(I)*xMatB(I)
	  end do  
	
	end subroutine MatAdpMatB
	
	subroutine dbltobool(A,B)
	!******************************************************************
	! Takes a double which values are either 1.0 or 0.0 and returns a *
	! Boolean
	!******************************************************************
	implicit none
	double precision, intent(in):: A
	logical, intent(out):: B
	if (A<1.0) then
		B=.false.
	else
		B=.true.
	endif
	end subroutine dbltobool

	function logic2dbl(a)
	  logical, intent(in) :: a

	  if (a) then
		logic2dbl = 1.d0
	  else
		logic2dbl = 0.d0
	  end if
	end function logic2dbl
	
	Subroutine check4crossing(IErate0I, IErateI, dErate_eff,RateRef, Apply)
		implicit none
		double precision, intent(inout):: IErate0I, IErateI, dErate_eff, RateRef
		logical:: cond1, cond2, cond3
		logical, intent(out)::Apply
		  Apply=.false.
		  cond1=(IErate0I<RateRef).and.(IErateI>RateRef)
		  if (cond1) IErate0I=RateRef	
		  cond2=(IErate0I>RateRef).and.(IErateI<RateRef)
		  if (cond2) IErateI=RateRef
		  dErate_eff=IErateI-IErate0I
		  cond3=(IErate0I>RateRef).and.(abs(dErate_eff)>0.0d0)
		  if (cond1.or.cond2.or.cond3) Apply=.true.
	end Subroutine check4crossing
