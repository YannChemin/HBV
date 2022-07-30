*******************************************************************************
*     HBV-7 Karkheh								*
*    										*
*******************************************************************************
*	-
*1=Doab					*
*		2=Pole_Char							*
*		3=Dobe_Merek							*
*		4=Ghor_B							*
*		5=Darttot_TS							*
*		6=Holilan							*
*		7=Kaka_R_Pole_d_Cham_A						*
*		8=Joligir							*
*		9=Pole_z_Paye_P							*
*										*
*	- Version 1.2 (3-2-2005): incorporation of sub-basin dependent initial	*
*		conditions,correction coefficients,elevation differences (between*
*		subbasin and station) and surface area fractions forests and fields *
*	- Version 1.3.1 (14-2-2005): incorporation of additional model 		*
*		performance criteria						*
*	- Version 1.3.2 (5-9-2005):	incorporation of output infordmation for all*
*		sub-basins and additional performance criterion			*
*	- Version 1.4 (30-9-2005): incorporation of uncertainty analysis for	*
*		climate change purposes and corresponding output information (see*
*		file 10)							*
*										*
*******************************************************************************
      program hbv15
	integer ntot,ttot,dtot,ytot,btot,rtot,TH,TL,TY,TS,tstart,q				! no. runs,no. time steps,no. time steps data,no. years,no. subbasins,no. discharge 'stations',temporal range moving average high flow,temporal range moving average low flow,days per year,seconds per day,warming-up period
      parameter (ntot=45000,ttot=916,dtot=1910,ytot=5,btot=8,rtot=8
     $	,TH=3,TL=15,TY=365,TS=86400,tstart=274)
	integer m(rtot),n,t,b,mtot(rtot),u,qotot(rtot),yotot(rtot)			! selected run number,run number,time step,subbasin number,maximum m,time integrator moving average,total number of observed discharges,total number of observed discharge years
	integer etotot(rtot)
	real MTDJF,MTMAM,MTJJA,MTSON,STDJF,STMAM,STJJA,STSON,					! parameters temperature for climate change (M=mean,S=standard deviation,DFJ...SON=seasons)
     $	MPDJF,MPMAM,MPJJA,MPSON,SPDJF,SPMAM,SPJJA,SPSON,					! parameters precipitation for climate change (M=mean,S=standard deviation,DFJ...SON=seasons)
     $	EDJF,EMAM,EJJA,ESON,												! parameters evaporation for climate change (DFJ...SON=seasons)
     $	TCON,CEVPFO,ECALT,PCALT,TCALT,SFCF,FOSFCF,RFCF,TT,TTI,DTTM,			! no. sec/day,10 HBV parameters
     $	CFMAX,FOCFMAX,CFR,WHC,RQH,RQL,CNS,CRVE,QDB,ATOT,KL,KU,				! 4 HBV parameters,relative lower limit high flow,relative upper limit low flow,lower limit Nash-Sutcliffe,upper limit volume error,critical low discharge,surface area Meuse basin,frequency factor for lower return value,frequency factor for upper return value
     $	LAG0103,LAG0204,LAG0304,LAG0504,LAG0604,LAG0407,LAG0708,			! fraction total discharge current time step (7 links)
     $	LAG0716,LAG0916,LAG1016,LAG1116,LAG1216,LAG1316,LAG1514,			! fraction total discharge current time step (7 links)
     $	ssm1,ssw1,sgw1,ssp1,smw1,sgwmax											! initial values for 5 storages
      parameter (MTDJF=0.0,MTMAM=0.0,MTJJA=0.0,MTSON=0.0,
     $	STDJF=0.,STMAM=0.,STJJA=0.,STSON=0.,
     $	MPDJF=0.0,MPMAM=0.0,MPJJA=0.0,MPSON=0.0,
     $	SPDJF=0.,SPMAM=0.,SPJJA=0.,SPSON=0.,
     $	EDJF=0.0564,EMAM=0.0350,EJJA=0.0384,ESON=0.0762,
     $	TCON=86.4,CEVPFO=1.15,ECALT=0.1,PCALT=0.2,TCALT=0.6,
     $	SFCF=1.0,FOSFCF=1.0,RFCF=1.0,TT=-0.5,TTI=2.0,
     $	DTTM=0.54391,CFMAX=3.5,FOCFMAX=0.6,CFR=0.05,WHC=0.1,
     $	RQH=5.,RQL=0.5,CNS=-100.00,CRVE=1000.0,QDB=100.,ATOT=21000.,
     $	KL=1.30456,KU=3.13668,												! lower is 10 years and upper is 100 years
     $	LAG0102=0.87,LAG0205=0.30,LAG0304=0.31,LAG0507=0.74,
     $	LAG0607=0.57,LAG0708=0.7,LAG0810=0.52,LAG0910=0.37,
     $	LAG1011=0.82,LAG1113=0.70,LAG1213=0.39,
     $	ssm1=15.,ssw1=15.,sgw1=15.,ssp1=0.,smw1=0.,sgwmax=5.0)
      real po(btot,dtot),etpo(btot,dtot),qo(rtot,dtot),tm(btot,dtot),		! observed precipitation,potential evapotranspiration,discharge,temperature
     $	pcc(ttot),tcc(ttot),ecc(ttot),tseas(ttot),eto(btot,dtot),					! transformation factors for climate change,season (1=DJF,2=MAM,3=JJA,4=SON)
     $	ssm(btot,ttot+1),ssw(btot,ttot+1),sgw(btot,ttot+1),				! storage soil moisture,surface water,ground water
     $	ssp(btot,ttot+1),smw(btot,ttot+1),sgwx(btot,ttot+1),						! storage snow pack,melting water
     $	p(btot,ttot),ta(btot,ttot),s(btot,ttot),r(btot,ttot),				! precipitation corrected,temperature corrected,snowfall,rainfall,
     $	sm(btot,ttot),sr(btot,ttot),in(btot,ttot),etp(btot,ttot),			! snow melt,refreezing rate,infiltration,potential evapotranspiration corrected
     $	qd(btot,ttot),qin(btot,ttot),qs(btot,ttot),qf(btot,ttot),			! direct discharge,indirect discharge,slow flow,fast flow
     $	qt(btot,ttot),qc(btot,ttot),qr(rtot,ttot),eta(btot,ttot),			! total discharge,capillary rise,river discharge,actual evapotranspiration
     $	A(btot),FFO(btot),FFI(btot),DEP(btot),DET(btot),DEE(btot),			! surface area,fraction forest,fraction field,elevation differences (100 m) subbasin and station for precipitation,temperature and evapotranspiration
     $	LFC(btot),HFC(btot),LBETA(btot),HBETA(btot),LLP(btot),				! lower (L) and upper (H) boundaries of HBV parameter ranges for uncertainty analysis
     $	HLP(btot),LALFA(btot),HALFA(btot),LKF(btot),HKF(btot),				! lower (L) and upper (H) boundaries of HBV parameter ranges for uncertainty analysis
     $	LKS(btot),HKS(btot),LPERC(btot),HPERC(btot),LCFLUX(btot),			! lower (L) and upper (H) boundaries of HBV parameter ranges for uncertainty analysis
     $	HCFLUX(btot),														! lower (L) and upper (H) boundaries of HBV parameter ranges for uncertainty analysis
     $	RTDJF(ntot),RTMAM(ntot),RTJJA(ntot),RTSON(ntot),					! random value U(0,1) for uncertainty analysis
     $	RPDJF(ntot),RPMAM(ntot),RPJJA(ntot),RPSON(ntot),					! random value U(0,1) for uncertainty analysis
     $	RFC(btot,ntot),RBETA(btot,ntot),RLP(btot,ntot),						! random value U(0,1) for uncertainty analysis
     $	RALFA(btot,ntot),RKF(btot,ntot),RKS(btot,ntot),						! random value U(0,1) for uncertainty analysis
     $	RPERC(btot,ntot),RCFLUX(btot,ntot),									! random value U(0,1) for uncertainty analysis
     $	TDJF(ntot),TMAM(ntot),TJJA(ntot),TSON(ntot),						! random value for climate parameters for uncertainty analysis
     $	PDJF(ntot),PMAM(ntot),PJJA(ntot),PSON(ntot),						! random value for climate parameters for uncertainty analysis
     $	FC(btot,ntot),BETA(btot,ntot),LP(btot,ntot),ALFA(btot,ntot),		! random value for HBV parameter for uncertainty analysis
     $	KF(btot,ntot),KS(btot,ntot),PERC(btot,ntot),CFLUX(btot,ntot),		! random value for HBV parameter for uncertainty analysis
     $	difq(rtot),dq(rtot),dqh(rtot),dql(rtot),mqo(rtot),meto(rtot),					! cumulative difference calculated and observed discharge,cumulative squared difference observed and calculated discharge (average,high flow,low flow),mean observed discharge
     $	difet(rtot),deta(rtot),veto(rtot),
     $	mqoh(rtot),mqol(rtot),vqo(rtot),vqoh(rtot),vqol(rtot),				! mean observed discharge (high flow,low flow),variance observed discharge (average,high flow,low flow)
     $	qrh(rtot,ttot),qoh(rtot,ttot),qrl(rtot,ttot),qol(rtot,ttot),		! moving average of discharge (calculated high flow,observed high flow,calculated low flow,observed low flow)
     $	dqha(rtot),dqla(rtot),qdo(ytot),qdr(ytot),							! cumulative relative absolute difference calculated and observed moving average discharge (high flow,low flow),observed annual discharge deficit (m3),simulated annual discharge deficit (m3)
     $	NS(rtot,ntot),NSH(rtot,ntot),NSL(rtot,ntot),RVE(rtot,ntot),			! Nash-Sutcliffe coefficient (average,high flow,low flow),relative volume error
     $  NSET(rtot,ntot),RVEET(rtot,ntot),				!****************
     $	RMAEH(rtot,ntot),RMAEL(rtot,ntot),DQD(ntot),						! relative mean absolute error (high flow,low flow),average error in annual simulated discharge deficit (m3)
     $	qodef(ntot),qrdef(ntot),											! average annual discharge deficit in m3 (observed and simulated)
     $	MFC(rtot,ntot),MBETA(rtot,ntot),MLP(rtot,ntot),						! selected values for 3 HBV parameters
     $	MALFA(rtot,ntot),MKF(rtot,ntot),MKS(rtot,ntot),						! selected values for 3 HBV parameters
     $	MPERC(rtot,ntot),MCFLUX(rtot,ntot),									! selected values for 2 HBV parameters
     $	MNS(rtot,ntot),MNSH(rtot,ntot),MNSL(rtot,ntot),MNSET(rtot,ntot),		! selected values for Nash-Sutcliffe coefficient (average,high flow,low flow)
     $	MRVE(rtot,ntot),MRMAEH(rtot,ntot),MRMAEL(rtot,ntot),		! selected values for relative volume error,relative mean absolute error (high flows,low flow)
     $	MRVEET(rtot,ntot),
     $	MDQD(ntot),MREVE(rtot,ntot),										! selected values for average error in annual simulated discharge deficit,average relative difference between observed and calculated return values for KL and KU
     $	ata(ttot),ap(ttot),mqr(ntot),mta(ntot),mp(ntot),				! areally averaged temperature and precipitation (surface area weighted and corrected for elevation), mean of simulated discharge, temperature and precipitation
     $	sqr(ntot),sqo(ntot),sta(ntot),sp(ntot),								! standard deviation of simulated discharge, observed discharge, temperature and precipitation
     $	tyear(ttot),qom(rtot,ytot),qrm(rtot,ytot),					! year number for time step,annual maximum observed discharge,annual maximum calculated discharge
     $	mqom(rtot),mqrm(rtot),sqom(rtot),sqmr(rtot),REVE(rtot,ntot)			! mean of annual maximum observed discharges,mean of annual maximum calculated discharges,standard deviation of annual maximum observed discharges,standard deviation of annual maximum calculated discharges,average relative difference between observed and calculated return values for KL and KU

*	open input and output files
	open (1,file='precip.dat')							! precipitation
	open (2,file='evap.dat')							! potential evapotranspiration
	open (3,file='dischargeobs.dat')						! observed discharge
	open (4,file='temp.dat')							! temperature
	open (7,file='param_sto.dat')							! lower and upper boundaries of HBV parameter ranges,surface area,fraction forest,fraction field,elevation differences (100 m) subbasin and station for precipitation,temperature and evapotranspiration
	open (8,file='etobs.dat')	!********************
	open (11,file='Output\01-Doab.csv')		! output stochastic Doab
	open (12,file='Output\02-Pole_Chehr.csv')	! output stochastic Pole Chehr
	open (13,file='Output\03-Doabe_M.csv')		! output stochastic Doabe M
	open (14,file='Output\04-Ghor_B.csv')		! output stochastic Ghor B
	open (15,file='Output\05-Holilan.csv')		! output stochastic Holilan
	open (16,file='Output\06-Pole_D.csv')		! output stochastic Pole D
	open (17,file='Output\07-Jelogir.csv')		! output stochastic Jelogir
	open (18,file='Output\08-Paye_P.csv')		! output stochastic Paye P
	open (21,file='Basinout\01det-Doab.csv')	! output deterministic Doab
	open (22,file='Basinout\02det-Pole_Chehr.csv')	! output deterministic Pole Chehr
	open (23,file='Basinout\03det-Doabe_M.csv')	! output deterministic Doabe M
	open (24,file='Basinout\04det-Ghor_B.csv')	! output deterministic Ghor B
	open (25,file='Basinout\05det-Holilan.csv')	! output deterministic Holilan
	open (26,file='Basinout\06det-Pole_D.csv')	! output deterministic Pole D
	open (27,file='Basinout\07det-Jelogir.csv')	! output deterministic Jelogir
	open (28,file='Basinout\08det-Paye_P.csv')	! output deterministic Paye P
	open (31,file='ETout\01et-Doab.csv')		! output ET Doab
	open (32,file='ETout\02et-Pole_Chehr.csv')	! output ET Pole Chehr
	open (33,file='ETout\03et-Doabe_M.csv')		! output ET Doabe M
	open (34,file='ETout\04et-Ghor_B.csv')		! output ET Ghor B
	open (35,file='ETout\05et-Holilan.csv')		! output ET Holilan
	open (36,file='ETout\06et-Pole_D.csv')		! output ET Pole D
	open (37,file='ETout\07et-Jelogir.csv')		! output ET Jelogir
	open (38,file='ETout\08et-Paye_P.csv')		! output ET Paye P
	open (41,file='EToutput\01etsto-Doab.csv')	! output ETstochastic Doab
	open (42,file='EToutput\02etsto-Pole_Chehr.csv')! output ETstochastic Pole Chehr
	open (43,file='EToutput\03etsto-Doabe_M.csv')	! output ETstochastic Doabe M
	open (44,file='EToutput\04etsto-Ghor_B.csv')	! output ETstochastic Ghor B
	open (45,file='EToutput\05etsto-Holilan.csv')	! output ETstochastic Holilan
	open (46,file='EToutput\06etsto-Pole_D.csv')	! output ETstochastic Pole D
	open (47,file='EToutput\07etsto-Jelogir.csv')	! output ETstochastic Jelogir
	open (48,file='EToutput\08etsto-Paye_P.csv')	! output ETstochastic Paye P

      
*     read input files (mm)
   30 read (1,*,end=35) po
      read (2,*,end=35) etpo
      read (3,*,end=35) qo
      read (4,*,end=35) tm
      read (8,*,end=35) eto
      goto 30 
   35 continue
	do 40 b=1,btot
		read (7,*,end=40) LFC,HFC,LBETA,HBETA,LLP,HLP,LALFA,HALFA,LKF,
     $	   HKF,LKS,HKS,LPERC,HPERC,LCFLUX,HCFLUX,A,FFO,FFI,DEP,DET,DEE
   40	continue
      close(7)
	do 45 b=1,rtot
		m(b)=1
   45 continue

	do 1100 n=1,ntot
	

*	| independent sampling of random values U(0,1) for HBV parameters
	do 50 b=1,btot
		call random_seed ()
		call random_number(harvest=RFC)
		call random_number(harvest=RBETA)
		call random_number(harvest=RLP)
		call random_number(harvest=RALFA)
		call random_number(harvest=RKF)
		call random_number(harvest=RKS)
		call random_number(harvest=RPERC)
		call random_number(harvest=RCFLUX)
*	| transformation of U(0,1) values to HBV parameter values
		FC(b,n)=LFC(b)+(HFC(b)-LFC(b))*RFC(b,n)
		BETA(b,n)=LBETA(b)+(HBETA(b)-LBETA(b))*RBETA(b,n)
		LP(b,n)=LLP(b)+(HLP(b)-LLP(b))*RLP(b,n)
		ALFA(b,n)=LALFA(b)+(HALFA(b)-LALFA(b))*RALFA(b,n)
		KF(b,n)=LKF(b)+(HKF(b)-LKF(b))*RKF(b,n)
		KS(b,n)=LKS(b)+(HKS(b)-LKS(b))*RKS(b,n)
		PERC(b,n)=LPERC(b)+(HPERC(b)-LPERC(b))*RPERC(b,n)
		CFLUX(b,n)=LCFLUX(b)+(HCFLUX(b)-LCFLUX(b))*RCFLUX(b,n)

*	water contents (mm) and fluxes (mm/d)/ (m3/s)
		ssp(b,1)=ssp1
		smw(b,1)=smw1
		ssm(b,1)=min(ssm1,FC(b,n))
		ssw(b,1)=ssw1
		sgw(b,1)=sgw1
		ssp(b,2)=ssp(b,1)
		smw(b,2)=smw(b,1)
		ssm(b,2)=ssm(b,1)
		ssw(b,2)=ssw(b,1)
		sgw(b,2)=sgw(b,1)
		qf(b,1)=KF(b,n)*(ssw(b,1)**(1+ALFA(b,n)))
		qs(b,1)=KS(b,n)*sgw(b,1)
		qt(b,1)=(qf(b,1)+qs(b,1))*A(b)/tcon
   50	continue

*	criteria
	do 100 b=1,rtot
		difq(b)=0.
		dq(b)=0.
		dqh(b)=0.
		dql(b)=0.
		dqha(b)=0.
		dqla(b)=0.
		deta(b)=0.	!****************
		mqo(b)=0.
		meto(b)=0.	!****************
		do 60 t=1,ttot
			qrh(b,t)=0.
			qoh(b,t)=0.
			qrl(b,t)=0.
			qol(b,t)=0.
   60		continue
		mqoh(b)=0.
		mqol(b)=0.
		vqo(b)=0.
		vqoh(b)=0.
		vqol(b)=0.
		veto(b)=0.	!****************
		qotot(b)=0.
		etotot(b)=0.	!****************
            difet(b)=0.
		do 70 u=1,ytot
			qom(b,u)=0.
			qrm(b,u)=0.
   70		continue
		yotot(b)=0
		mqom(b)=0.
		mqrm(b)=0.
		sqom(b)=0.
		sqmr(b)=0.
  100 continue
	do 150 t=1,ytot
		qdo(t)=0
		qdr(t)=0
  150	continue
	DQD(n)=0

*	water fluxes (m3/s)
*	qr(2,1)=qt(1,1)
* 	qr(4,1)=qt(3,1)
*	qr(5,1)=qr(2,1)
*	qr(7,1)=qr(5,1)
*	qr(8,1)=qt(8,1)
     	
*	DYNAMIC PART
	do 500 t=2,ttot
*		tseas(t)=ceiling(mod(real(t)+761.,365.25)/91.3125)
*		if(tseas(t).eq.1) then
*			tcc(t)=TDJF(n)
*			pcc(t)=PDJF(n)
*			ecc(t)=EDJF
*		elseif(tseas(t).eq.2) then
*			tcc(t)=TMAM(n)
*			pcc(t)=PMAM(n)
*			ecc(t)=EMAM
*		elseif(tseas(t).eq.3) then
*			tcc(t)=TJJA(n)
*			pcc(t)=PJJA(n)
*			ecc(t)=EJJA
*		elseif(tseas(t).eq.4) then
*			tcc(t)=TSON(n)
*			pcc(t)=PSON(n)
*			ecc(t)=ESON
*		endif
		tcc(t)=0.
		pcc(t)=0.
		ecc(t)=0.
	do 250 b=1,btot

*	|snow pack balance components (mm/d)
		p(b,t)=(1+0)*po(b,t)*(1+PCALT*DEP(b))
		ta(b,t)=tm(b,t)+0-TCALT*DET(b)
		if(ta(b,t).lt.(TT-TTI/2)) then
			s(b,t)=p(b,t)*(FFO(b)*FOSFCF+FFI(b))*SFCF
		elseif (ta(b,t).ge.(TT-TTI/2).and.ta(b,t).lt.(TT+TTI/2)) then
			s(b,t)=p(b,t)*((TT+TTI/2)-ta(b,t))/TTI*(FFO(b)*FOSFCF+
     $			FFI(b))*SFCF
			r(b,t)=p(b,t)*(ta(b,t)-(TT-TTI/2))/TTI*RFCF
		else
			r(b,t)=p(b,t)*RFCF
		endif
		sm(b,t)=min(max(CFMAX*(ta(b,t)-(TT+DTTM)),0.),ssp(b,t))
		sr(b,t)=min(max(CFR*(FFO(b)*FOCFMAX+FFI(b))*CFMAX*((TT+DTTM)-
     $		ta(b,t)),0.),smw(b,t))
		in(b,t)=max(smw(b,t)+sm(b,t)+r(b,t)-sr(b,t)-WHC*ssp(b,t),0.)

*	|soil moisture balance components (mm/d)
		qd(b,t)=max((in(b,t)+ssm(b,t)-FC(b,n)),0.)
		qin(b,t)=((ssm(b,t)/FC(b,n))**BETA(b,n))*(in(b,t)-qd(b,t))
		etp(b,t)=(1+0*ecc(t))*etpo(b,t)*(FFO(b)*CEVPFO+FFI(b))
     $		*(1-ECALT*DEE(b))
		eta(b,t)=min(etp(b,t),(etp(b,t)*ssm(b,t)/(LP(b,n)*FC(b,n))))
		qc(b,t)=CFLUX(b,n)*(FC(b,n)-ssm(b,t))/FC(b,n)

*	|surface water balance components (mm/d)
		qf(b,t)=KF(b,n)*(ssw(b,t)**(1+ALFA(b,n)))

*	|ground water balance components (mm/d)
		qs(b,t)=KS(b,n)*sgw(b,t)

*     |total discharge (m3/s)
		qt(b,t)=(qs(b,t)+qf(b,t))*A(b)/tcon

*	|snow pack balance (mm)
		ssp(b,t+1)=ssp(b,t)+s(b,t)+sr(b,t)-sm(b,t)
		smw(b,t+1)=smw(b,t)+r(b,t)-sr(b,t)+sm(b,t)-in(b,t)

*	|surface water balance (mm)
		if(sgw(b,t)>=sgwmax) then

		ssw(b,t+1)=max(ssw(b,t)+max((qd(b,t)+qin(b,t)
     $          -0.0),0.)
     $	        -KF(b,n)*(ssw(b,t)**(1+ALFA(b,n)))
     $          -min(ssw(b,t),qc(b,t)),0.)
               else
                ssw(b,t+1)=max(ssw(b,t)+max((qd(b,t)+qin(b,t)
     $          -PERC(b,n)),0.)
     $	        -KF(b,n)*(ssw(b,t)**(1+ALFA(b,n)))
     $          -min(ssw(b,t),qc(b,t)),0.)

                 
               endif

*	|soil moisture balance (mm)
		if(ssw(b,t+1).eq.0.) then
			qc(b,t)=ssw(b,t)+max((qd(b,t)+qin(b,t)-PERC(b,n)),0.)
     $			-KF(b,n)*(ssw(b,t)**(1+ALFA(b,n)))
		else
			qc(b,t)=min(ssw(b,t),qc(b,t))
		endif
		ssm(b,t+1)=ssm(b,t)+in(b,t)-qd(b,t)-qin(b,t)+qc(b,t)-eta(b,t)

*	|ground water balance (mm)

		if(sgw(b,t)>=sgwmax) then

                sgw(b,t+1)=(1-KS(b,n))*sgw(b,t)
                        
                else
                sgwx(b,t+1)=(1-KS(b,n))*sgw(b,t)+min((qd(b,t)+qin(b,t)),
     $		PERC(b,n))
                
                sgw(b,t+1)=min(sgwx(b,t+1),sgwmax)

                endif
  

  250 continue

*     |river discharge (m3/s)
		qr(2,t)=qt(2,t)
		qr(4,t)=qt(4,t)
    	      qr(5,t)=qt(5,t)
     		qr(7,t)=qt(7,t)
    		qr(8,t)=qt(8,t)
     		qr(1,t)=qt(1,t)
		qr(3,t)=qt(3,t)
		qr(6,t)=qt(6,t)
		

		
  500 continue

*     criteria
	do 1000 b=1,rtot
*	| calculation of Nash-Sutcliffe coefficient (NS),weighted high flow NS,weighted low flow NS,relative volume error (RVE),high flow relative mean absolute error (RMAEH),low flow relative mean absolute error (RMAEL),average relative difference between observed and calculated return values for KL and KU (REVE)
		do 940 t=tstart,ttot
			if (qo(b,t).ge.0.) then
				mqo(b)=mqo(b)+qo(b,t)
				qotot(b)=qotot(b)+1
			endif
			if (eto(b,t)>0.) then
				meto(b)=meto(b)+eto(b,t)	!*******cumulative observed ETA
				etotot(b)=etotot(b)+1		!*******total observed ETA
			endif

  940		continue
		mqo(b)=mqo(b)/qotot(b)
		meto(b)=meto(b)/etotot(b)			!********mean observed ETA

		do 950 t=tstart,ttot
			if (qo(b,t).ge.0.) then
				difq(b)=difq(b)+(qr(b,t)-qo(b,t))
				dq(b)=dq(b)+(qo(b,t)-qr(b,t))**2
				dqh(b)=dqh(b)+qo(b,t)*((qo(b,t)-qr(b,t))**2)
				dql(b)=dql(b)+(1/(1+qo(b,t)))*((qo(b,t)-qr(b,t))**2)
				vqo(b)=vqo(b)+(qo(b,t)-mqo(b))**2
				vqoh(b)=vqoh(b)+qo(b,t)*((qo(b,t)-mqo(b))**2)
				vqol(b)=vqol(b)+(1/(1+qo(b,t)))*((qo(b,t)-mqo(b))**2)
			endif
			if (eto(b,t)>0.) then
				difet(b)=difet(b)+(eta(b,t)-eto(b,t))		!*******cumulative difference calculated and observed ETA
				deta(b)=deta(b)+(eto(b,t)-eta(b,t))**2		!*******cumulative squared difference observed and calculated ETA
				veto(b)=veto(b)+(eto(b,t)-meto(b))**2		!*******variance observed ETA
			endif

  950		continue
		do 960 t=(tstart+TH),(ttot-TH+1)
			do 955 u=(t-TH+1),(t+TH-1)
				qrh(b,t)=qrh(b,t)+qr(b,u)
				if (qo(b,u).ge.0.) then
					qoh(b,t)=qoh(b,t)+qo(b,u)
				endif
  955			continue
			qrh(b,t)=qrh(b,t)/(2*TH-1)
			qoh(b,t)=qoh(b,t)/(2*TH-1)
			if (qoh(b,t).gt.(RQH*mqo(b))) then
				dqha(b)=dqha(b)+abs(qrh(b,t)-qoh(b,t))
				mqoh(b)=mqoh(b)+qoh(b,t)
			endif
  960		continue
		do 970 t=(tstart+TL),(ttot-TL+1)
			do 965 u=(t-TL+1),(t+TL-1)
				qrl(b,t)=qrl(b,t)+qr(b,u)
				if (qo(b,u).ge.0.) then
					qol(b,t)=qol(b,t)+qo(b,u)
				endif
  965			continue
			qrl(b,t)=qrl(b,t)/(2*TL-1)
			qol(b,t)=qol(b,t)/(2*TL-1)
			if (qol(b,t).le.(RQL*mqo(b))) then
				dqla(b)=dqla(b)+abs(qrl(b,t)-qol(b,t))
				mqol(b)=mqol(b)+qol(b,t)
			endif
  970		continue
		do 975 t=1,ttot
			tyear(t)=ceiling((real(t)-0.75)/365.25)
  975		continue
		do 985 u=1,ytot
			do 980 t=tstart,ttot
				if (tyear(t).eq.u) then
					qom(b,u)=max(qom(b,u),qo(b,t))
					qrm(b,u)=max(qrm(b,u),qr(b,t))
				endif					 
  980			continue
			if(qom(b,u).gt.0.) then
				yotot(b)=yotot(b)+1
			endif
			mqom(b)=mqom(b)+qom(b,u)
			mqrm(b)=mqrm(b)+qrm(b,u)
  985		continue
		mqom(b)=mqom(b)/yotot(b)
		mqrm(b)=mqrm(b)/ytot
		do 990 u=1,ytot
			if (qom(b,u).gt.0.) then
				sqom(b)=sqom(b)+(qom(b,u)-mqom(b))**2
			endif
			sqmr(b)=sqmr(b)+(qrm(b,u)-mqrm(b))**2
  990		continue
		sqom(b)=sqrt(sqom(b)/(max(yotot(b)-1,1)))
		sqmr(b)=sqrt(max(1.,sqmr(b))/(ytot-1))
				
		NS(b,n)=1-dq(b)/vqo(b)
		NSET(b,n)=1-deta(b)/veto(b)	!*******Nash-Sutcliffe coefficient-ETA
		NSH(b,n)=1-dqh(b)/vqoh(b)
		NSL(b,n)=1-dql(b)/vqol(b)
		RVE(b,n)=100*difq(b)/(mqo(b)*qotot(b))
		RVEET(b,n)=100*difet(b)/(meto(b)*etotot(b))	!*******relative volume error-ETA
		RMAEH(b,n)=dqha(b)/mqoh(b)
		RMAEL(b,n)=dqla(b)/mqol(b)
		REVE(b,n)=100*(((mqrm(b)-mqom(b)+(sqmr(b)-sqom(b))*KL)/
     $		(mqom(b)+sqom(b)*KL))+((mqrm(b)-mqom(b)+(sqmr(b)-sqom(b))
     $		*KU)/(mqom(b)+sqom(b)*KU)))/2
 1000 continue
*	NS(1,n)=NS(2,n)
*	NS(3,n)=NS(2,n)
*	NS(4,n)=NS(2,n)
*	NS(6,n)=NS(2,n)
*	NS(9,n)=NS(2,n)
*	NS(11,n)=NS(2,n)
*	NS(12,n)=NS(2,n)
*	NS(13,n)=NS(2,n)
*	NS(14,n)=NS(2,n)
*	NS(15,n)=NS(2,n)
*	NS(16,n)=NS(2,n)
*	RVE(1,n)=RVE(2,n)
*	RVE(3,n)=RVE(2,n)
*	RVE(4,n)=RVE(2,n)
*	RVE(6,n)=RVE(2,n)
*	RVE(9,n)=RVE(2,n)
*	RVE(11,n)=RVE(2,n)
*	RVE(12,n)=RVE(2,n)
*	RVE(13,n)=RVE(2,n)
*	RVE(14,n)=RVE(2,n)
*	RVE(15,n)=RVE(2,n)
*	RVE(16,n)=RVE(2,n)
*	| calculation of discharge deficit for complete basin
	qodef(n)=0.
	qrdef(n)=0.
	do 1020 t=1,ytot
		do 1010 u=(1+(t-1)*TY),(TY*t)
			qdo(t)=qdo(t)+max(QDB-qo(8,u),0.)*TS
			qdr(t)=qdr(t)+max(QDB-qr(8,u),0.)*TS
 1010		continue
		DQD(n)=DQD(n)+(qdo(t)-qdr(t))/qdo(t)
		qodef(n)=qodef(n)+qdo(t)
		qrdef(n)=qrdef(n)+qdr(t)
 1020	continue
	DQD(n)=100*DQD(n)/ytot
	qodef(n)=qodef(n)/ytot
	qrdef(n)=qrdef(n)/ytot

*	| selection of parameter sets based on boundaries for criteria (CNS and CRVE)
	do 1040 b=1,btot
*		if(NS(b,n).ge.CNS.and.abs(RVE(b,n)).le.CRVE) then
			MFC(b,m(b))=FC(b,n)
			MBETA(b,m(b))=BETA(b,n)
			MLP(b,m(b))=LP(b,n)
			MALFA(b,m(b))=ALFA(b,n)
			MKF(b,m(b))=KF(b,n)
			MKS(b,m(b))=KS(b,n)
			MPERC(b,m(b))=PERC(b,n)
			MCFLUX(b,m(b))=CFLUX(b,n)
*		endif
 1040 continue
	do 1050 b=1,rtot
*		if(NS(b,n).ge.CNS.and.abs(RVE(b,n)).le.CRVE) then
			MNS(b,m(b))=NS(b,n)
			MNSH(b,m(b))=NSH(b,n)
			MNSL(b,m(b))=NSL(b,n)
			MNSET(b,m(b))=NSET(b,n)		!**************
			MRVE(b,m(b))=RVE(b,n)
			MRVEET(b,m(b))=RVEET(b,n)	!**************
			MRMAEH(b,m(b))=RMAEH(b,n)
			MRMAEL(b,m(b))=RMAEL(b,n)
			MREVE(b,m(b))=REVE(b,n)
			MDQD(m(b))=DQD(n)
			mtot(b)=m(b)
			m(b)=m(b)+1
*		endif
 1050 continue

**	| calculation of mean and standard deviation for discharge, temperature and precipitation for complete basin
*	do 1060 t=1,ttot
*		ata(t)=0.
*		ap(t)=0.
*		do 1055 b=1,btot
*			ata(t)=ata(t)+(A(b)/ATOT)*ta(b,t)
*			ap(t)=ap(t)+(A(b)/ATOT)*p(b,t)
* 1055		continue
* 1060	continue
*	mqr(n)=0.
*	mp(n)=0.
*	mta(n)=0.
*	do 1065 t=1,ttot
*		mqr(n)=mqr(n)+qr(2,t)
*		mta(n)=mta(n)+ata(t)
*		mp(n)=mp(n)+ap(t)
* 1065	continue
*	mqr(n)=mqr(n)/ttot
*	mp(n)=mp(n)/ttot
*	mta(n)=mta(n)/ttot
*	sqr(n)=0.
*	sqo(n)=0.
*	sta(n)=0.
*	sp(n)=0.
*	do 1070 t=1,ttot
*		sqr(n)=sqr(n)+(qr(2,t)-mqr(n))**2
*		sqo(n)=sqo(n)+(qo(2,t)-mqo(2))**2
*		sta(n)=sta(n)+(ata(t)-mta(n))**2
*		sp(n)=sp(n)+(ap(t)-mp(n))**2
* 1070 continue
*	sqr(n)=sqrt(sqr(n)/(ttot-1))
*	sqo(n)=sqrt(sqo(n)/(ttot-1))
*	sta(n)=sqrt(sta(n)/(ttot-1))
*	sp(n)=sqrt(sp(n)/(ttot-1))
	print '(I10)',n
 1100	continue
  
**     OUTPUT FILES

*      write (10,'(9(A10,A1),A10)') 'mp',',','sp',',','mta',',','sta',','
*     $	,'mqo',',','mqr',',','sqo',',','sqr',',','qodef',',','qrdef'
*	do 8900 n=1,ntot
*		write (10,'(9(F16.4,A1),F16.4)') mp(n),',',sp(n),',',mta(n),
*     $		',',sta(n),',',mqo(2),',',mqr(n),',',sqo(n),',',sqr(n),
*     $		',',qodef(n),',',qrdef(n)
* 8900	continue
*      write (11,'(14(A10,A1),A10)') ' ',',','FC(mm)',',','BETA(-)',',',
*     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
*     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
*     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
*	do 9000 n=1,ntot
*		write (11,'(A10,A1,13(F10.4,A1),F10.4)') 'Meuse',',',FC(2,n),
*     $		',',BETA(2,n),',',LP(2,n),',',ALFA(2,n),',',KF(2,n),',',
*     $		KS(2,n),',',PERC(2,n),',',CFLUX(2,n),',',NS(2,n),',',
*     $		NSH(2,n),',',NSL(2,n),',',RVE(2,n),',',RMAEH(2,n),',',
*     $		RMAEL(2,n)
*		write (11,'(A10,A1,13(F10.4,A1),F10.4)') 'Vesdre',',',FC(5,n),
*     $		',',BETA(5,n),',',LP(5,n),',',ALFA(5,n),',',KF(5,n),',',
*     $		KS(5,n),',',PERC(5,n),',',CFLUX(5,n),',',NS(5,n),',',
*     $		NSH(5,n),',',NSL(5,n),',',RVE(5,n),',',RMAEH(5,n),',',
*     $		RMAEL(5,n)
*		write (11,'(A10,A1,13(F10.4,A1),F10.4)') 'Ourthe',',',FC(7,n),
*     $		',',BETA(7,n),',',LP(7,n),',',ALFA(7,n),',',KF(7,n),',',
*     $		KS(7,n),',',PERC(7,n),',',CFLUX(7,n),',',NS(7,n),',',
*     $		NSH(7,n),',',NSL(7,n),',',RVE(7,n),',',RMAEH(7,n),',',
*     $		RMAEL(7,n)
*		write (11,'(A10,A1,13(F10.4,A1),F10.4)') 'Ambleve',',',FC(8,n),
*     $		',',BETA(8,n),',',LP(8,n),',',ALFA(8,n),',',KF(8,n),',',
*     $		KS(8,n),',',PERC(8,n),',',CFLUX(8,n),',',NS(8,n),',',
*     $		NSH(8,n),',',NSL(8,n),',',RVE(8,n),',',RMAEH(8,n),',',
*     $		RMAEL(8,n)
*		write (11,'(A10,A1,13(F10.4,A1),F10.4)') 'Lesse',',',FC(10,n),
*     $		',',BETA(10,n),',',LP(10,n),',',ALFA(10,n),',',KF(10,n),
*     $		',',KS(10,n),',',PERC(10,n),',',CFLUX(10,n),',',NS(10,n),
*     $		',',NSH(10,n),',',NSL(10,n),',',RVE(10,n),',',RMAEH(10,n),
*     $		',',RMAEL(10,n)
* 9000	continue
       do 9100 t=1,ttot

		write (21,'(3(F8.2,A1),F8.2)') qo(1,t),',',qr(1,t),',',qf(1,t)
		write (22,'(3(F8.2,A1),F8.2)') qo(2,t),',',qr(2,t),',',po(2,t)
		write (23,'(3(F8.2,A1),F8.2)') qo(3,t),',',qr(3,t),',',po(3,t)
		write (24,'(3(F8.2,A1),F8.2)') qo(4,t),',',qr(4,t),',',po(4,t)
		write (25,'(3(F8.2,A1),F8.2)') qo(5,t),',',qr(5,t),',',po(5,t)
		write (26,'(3(F8.2,A1),F8.2)') qo(6,t),',',qr(6,t),',',po(6,t)
		write (27,'(3(F8.2,A1),F8.2)') qo(7,t),',',qr(7,t),',',po(7,t)
		write (28,'(3(F8.2,A1),F8.2)') qo(8,t),',',qr(8,t),',',po(8,t)

		write (31,'(2(F8.2,A1),F8.2)') eta(1,t),',',eto(1,t)
		write (32,'(2(F8.2,A1),F8.2)') eta(2,t),',',eto(2,t)
		write (33,'(2(F8.2,A1),F8.2)') eta(3,t),',',eto(3,t)
		write (34,'(2(F8.2,A1),F8.2)') eta(4,t),',',eto(4,t)
		write (35,'(2(F8.2,A1),F8.2)') eta(5,t),',',eto(5,t)
		write (36,'(2(F8.2,A1),F8.2)') eta(6,t),',',eto(6,t)
		write (37,'(2(F8.2,A1),F8.2)') eta(7,t),',',eto(7,t)
		write (38,'(2(F8.2,A1),F8.2)') eta(8,t),',',eto(8,t)




 9100 continue
*	| Doab
	write (11,'(A30)') 'Doab'
	write (11,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (41,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9500 q=1,mtot(1)
		write (11,'(14(F12.6,A1),F12.6)') MFC(1,q),',',MBETA(1,q),
     $		',',MLP(1,q),',',MALFA(1,q),',',MKF(1,q),',',MKS(1,q),
     $		',',MPERC(1,q),',',MCFLUX(1,q),',',MNS(1,q),',',
     $		MNSH(1,q),',',MNSL(1,q),',',MRVE(1,q),',',MRMAEH(1,q),
     $		',',MRMAEL(1,q),',',MREVE(1,q)
     		write (41,'(10(F12.6,A1),F12.6)') MFC(1,q),',',MBETA(1,q),
     $		',',MLP(1,q),',',MALFA(1,q),',',MKF(1,q),',',MKS(1,q),
     $		',',MPERC(1,q),',',MCFLUX(1,q),',',MNSET(1,q),',',MRVEET(1,q)
 9500	continue
*	| Pole_c
	write (12,'(A30)') 'Doab'
	write (12,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (42,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9510 q=1,mtot(1)
		write (12,'(14(F12.6,A1),F12.6)') MFC(2,q),',',MBETA(2,q),
     $		',',MLP(2,q),',',MALFA(2,q),',',MKF(2,q),',',MKS(2,q),
     $		',',MPERC(2,q),',',MCFLUX(2,q),',',MNS(2,q),',',
     $		MNSH(2,q),',',MNSL(2,q),',',MRVE(2,q),',',MRMAEH(2,q),
     $		',',MRMAEL(2,q),',',MREVE(2,q)
     		write (42,'(10(F12.6,A1),F12.6)') MFC(2,q),',',MBETA(2,q),
     $		',',MLP(2,q),',',MALFA(2,q),',',MKF(2,q),',',MKS(2,q),
     $		',',MPERC(2,q),',',MCFLUX(2,q),',',MNSET(2,q),',',MRVEET(2,q)
9510	continue
*	| Doab_m
	write (13,'(A30)') 'Doab'
	write (13,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (43,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9520 q=1,mtot(1)
		write (13,'(14(F12.6,A1),F12.6)') MFC(3,q),',',MBETA(3,q),
     $		',',MLP(3,q),',',MALFA(3,q),',',MKF(3,q),',',MKS(3,q),
     $		',',MPERC(3,q),',',MCFLUX(3,q),',',MNS(3,q),',',
     $		MNSH(3,q),',',MNSL(3,q),',',MRVE(3,q),',',MRMAEH(3,q),
     $		',',MRMAEL(3,q),',',MREVE(3,q)
     		write (43,'(10(F12.6,A1),F12.6)') MFC(3,q),',',MBETA(3,q),
     $		',',MLP(3,q),',',MALFA(3,q),',',MKF(3,q),',',MKS(3,q),
     $		',',MPERC(3,q),',',MCFLUX(3,q),',',MNSET(3,q),',',MRVEET(3,q)
9520	continue
*	| Ghor_b
	write (14,'(A30)') 'Doab'
	write (14,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (44,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9530 q=1,mtot(1)
		write (14,'(14(F12.6,A1),F12.6)') MFC(4,q),',',MBETA(4,q),
     $		',',MLP(4,q),',',MALFA(4,q),',',MKF(4,q),',',MKS(4,q),
     $		',',MPERC(4,q),',',MCFLUX(4,q),',',MNS(4,q),',',
     $		MNSH(4,q),',',MNSL(4,q),',',MRVE(4,q),',',MRMAEH(4,q),
     $		',',MRMAEL(4,q),',',MREVE(4,q)
     		write (44,'(10(F12.6,A1),F12.6)') MFC(4,q),',',MBETA(4,q),
     $		',',MLP(4,q),',',MALFA(4,q),',',MKF(4,q),',',MKS(4,q),
     $		',',MPERC(4,q),',',MCFLUX(4,q),',',MNSET(4,q),',',MRVEET(4,q)
9530	continue
*	| Holilan
	write (15,'(A30)') 'Doab'
	write (15,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (45,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9540 q=1,mtot(1)
		write (15,'(14(F12.6,A1),F12.6)') MFC(5,q),',',MBETA(5,q),
     $		',',MLP(5,q),',',MALFA(5,q),',',MKF(5,q),',',MKS(5,q),
     $		',',MPERC(5,q),',',MCFLUX(5,q),',',MNS(5,q),',',
     $		MNSH(5,q),',',MNSL(5,q),',',MRVE(5,q),',',MRMAEH(5,q),
     $		',',MRMAEL(5,q),',',MREVE(5,q)
     		write (45,'(10(F12.6,A1),F12.6)') MFC(5,q),',',MBETA(5,q),
     $		',',MLP(5,q),',',MALFA(5,q),',',MKF(5,q),',',MKS(5,q),
     $		',',MPERC(5,q),',',MCFLUX(5,q),',',MNSET(5,q),',',MRVEET(5,q)
9540	continue

*	| Pole_D
	write (16,'(A30)') 'Doab'
	write (16,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (46,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9550 q=1,mtot(1)
		write (16,'(14(F12.6,A1),F12.6)') MFC(6,q),',',MBETA(6,q),
     $		',',MLP(6,q),',',MALFA(6,q),',',MKF(6,q),',',MKS(6,q),
     $		',',MPERC(6,q),',',MCFLUX(6,q),',',MNS(6,q),',',
     $		MNSH(6,q),',',MNSL(6,q),',',MRVE(6,q),',',MRMAEH(6,q),
     $		',',MRMAEL(6,q),',',MREVE(6,q)
     		write (46,'(10(F12.6,A1),F12.6)') MFC(6,q),',',MBETA(6,q),
     $		',',MLP(6,q),',',MALFA(6,q),',',MKF(6,q),',',MKS(6,q),
     $		',',MPERC(6,q),',',MCFLUX(6,q),',',MNSET(6,q),',',MRVEET(6,q)
9550	continue

*	| Jelogir
	write (17,'(A30)') 'Doab'
	write (17,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (47,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9560 q=1,mtot(1)
		write (17,'(14(F12.6,A1),F12.6)') MFC(7,q),',',MBETA(7,q),
     $		',',MLP(7,q),',',MALFA(7,q),',',MKF(7,q),',',MKS(7,q),
     $		',',MPERC(7,q),',',MCFLUX(7,q),',',MNS(7,q),',',
     $		MNSH(7,q),',',MNSL(7,q),',',MRVE(7,q),',',MRMAEH(7,q),
     $		',',MRMAEL(7,q),',',MREVE(7,q)
     		write (47,'(10(F12.6,A1),F12.6)') MFC(7,q),',',MBETA(7,q),
     $		',',MLP(7,q),',',MALFA(7,q),',',MKF(7,q),',',MKS(7,q),
     $		',',MPERC(7,q),',',MCFLUX(7,q),',',MNSET(7,q),',',MRVEET(7,q)
9560	continue

*	| Paye_pole
	write (18,'(A30)') 'Doab'
	write (18,'(14(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NS(-)',',','NSH(-)',
     $		',','NSL(-)',',','RVE(%)',',','RMAEH(-)',',','RMAEL(-)'
     $		,',','REVE(%)'
	write (48,'(10(A10,A1),A10)') 'FC(mm)',',','BETA(-)',',',
     $		'LP(-)',',','ALFA(-)',',','KF(1/d)',',','KS(1/d)',',',
     $		'PERC(mm/d)',',','CFLUX(mm/d)',',','NSET(-)',',','RVEET(%)'

	do 9570 q=1,mtot(1)
		write (18,'(14(F12.6,A1),F12.6)') MFC(8,q),',',MBETA(8,q),
     $		',',MLP(8,q),',',MALFA(8,q),',',MKF(8,q),',',MKS(8,q),
     $		',',MPERC(8,q),',',MCFLUX(8,q),',',MNS(8,q),',',
     $		MNSH(8,q),',',MNSL(8,q),',',MRVE(8,q),',',MRMAEH(8,q),
     $		',',MRMAEL(8,q),',',MREVE(8,q)
     		write (48,'(10(F12.6,A1),F12.6)') MFC(8,q),',',MBETA(8,q),
     $		',',MLP(8,q),',',MALFA(8,q),',',MKF(8,q),',',MKS(8,q),
     $		',',MPERC(8,q),',',MCFLUX(8,q),',',MNSET(8,q),',',MRVEET(8,q)
9570	continue


	close (1)
	close (2)
	close (3)
	close (4)
	close (8)
	close (11)
	close (12)
	close (13)
	close (14)
	close (15)
	close (16)
	close (17)
	close (18)
	close (21)
	close (22)
	close (23)
	close (24)
	close (25)
	close (26)
	close (27)
	close (28)
	close (31)
	close (32)
	close (33)
	close (34)
	close (35)
	close (36)
	close (37)
	close (38)
	close (41)
	close (42)
	close (43)
	close (44)
	close (45)
	close (46)
	close (47)
	close (48)
      end
