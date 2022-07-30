#include<math.h>
#include<sys/param.h> /*MIN() && MAX()*/

int hbv_model(int n,int t,int b,float **p,float **po,float pcalt,float *dep,float **ta,float **tm,float tcalt,float *det,float tt,float tti,float **s,float *ffo,float foscf,float *ffi,float sfcf,float **r,float rfcf,float **sm,float cfmax,float dttm,float **ssp,float cfr,float focfmax,float **smw,float **inp,float **sr,float whc,float **qd,float **ssm,float **fc,float **qin,float **beta,float **etp,float **etpo,float ecevpfo,float ecalt,float *dee,float **lp,float **qc,float **cflux,float **qf,float **kf,float **ssw,float **alpha,float **qs,float **ks,float **sgw,float **qt,float *area,float tcon,float **sgwx,float sgwmax,float **perc,float **eta){
	// snow pack balance components (mm/d)
	p[b][t]		= (1+0)*po[b][t]*(1+pcalt*dep[b]);		// Error
	ta[b][t]	= tm[b][t]+0-tcalt*det[b] ;			// Error
	if( ta[b][t] < (tt-tti/2)){
		s[b][t]	= p[b][t]*(ffo[b]*foscf+ffi[b])*sfcf;
	} else if (ta[b][t] >= (tt-tti/2) && ta[b][t] < (tt+tti/2)){
		s[b][t]	= p[b][t] * ((tt+tti/2)-ta[b][t])/tti*(ffo[b]*foscf+ffi[b])*sfcf;
		r[b][t]	= p[b][t] * ( ta[b][t] - (tt-tti/2 ) ) / tti * rfcf;
	}else{
		r[b][t]	= p[b][t] * rfcf;
	}
	sm[b][t]	= MIN(MAX(cfmax*(ta[b][t]-(tt+dttm)),0.),ssp[b][t]);
	sr[b][t]	= MIN(MAX(cfr*(ffo[b]*focfmax+ffi[b])*cfmax*((tt+dttm)-ta[b][t]),0.),smw[b][t]);
	inp[b][t]	= MAX(smw[b][t]+sm[b][t]+r[b][t]-sr[b][t]-whc*ssp[b][t],0.);
	// soil moisture balance components (mm/d)
	qd[b][t]	= MAX((inp[b][t]+ssm[b][t]-fc[b][n]),0.);
	qin[b][t]	= pow(ssm[b][t]/fc[b][n],beta[b][n]) * (inp[b][t]-qd[b][t]);
	etp[b][t]	= etpo[b][t] * (ffo[b] * ecevpfo+ffi[b]) * (1-ecalt*dee[b]);
	eta[b][t]	= MIN( etp[b][t], ( etp[b][t]*ssm[b][t] / ( lp[b][n] * fc[b][n] ) ) );
	qc[b][t]	= cflux[b][n] * ( fc[b][n]-ssm[b][t] ) / fc[b][n];
	// surface water balance components (mm/d)
	qf[b][t]	= kf[b][n] * pow( ssw[b][t], 1+alpha[b][n] );
	// ground water balance components (mm/d)
	qs[b][t]	= ks[b][n] * sgw[b][t];
	// total discharge (m3/s)
	qt[b][t]	= ( qs[b][t] + qf[b][t] ) * area[b] / tcon;
	// snow pack balance (mm)
	ssp[b][t+1]	= ssp[b][t] + s[b][t] + sr[b][t] - sm[b][t];
	smw[b][t+1]	= smw[b][t] + r[b][t] - sr[b][t] + sm[b][t] - inp[b][t];
	// surface water balance (mm)
	if(sgw[b][t] >= sgwmax){
		ssw[b][t+1] = MAX(ssw[b][t]+MAX((qd[b][t]+qin[b][t]),0.)-qf[b][t]-MIN(ssw[b][t],qc[b][t]),0.);
	}else{
		ssw[b][t+1] = MAX(ssw[b][t]+MAX((qd[b][t]+qin[b][t]-perc[b][n]),0.)-qf[b][t]-MIN(ssw[b][t],qc[b][t]),0.);
	}
	//soil moisture balance (mm)
	if(ssw[b][t+1]==0.){
		qc[b][t] = ssw[b][t]+MAX( (qd[b][t]+qin[b][t]-perc[b][n]),0.)-qf[b][t];
	}else{
		qc[b][t] = MIN( ssw[b][t],qc[b][t] );
	}
	ssm[b][t+1] = ssm[b][t]+inp[b][t]-qd[b][t]-qin[b][t]+qc[b][t]-eta[b][t];
	//ground water balance (mm)
	if(sgw[b][t] >= sgwmax){
		sgw[b][t+1] = (1-ks[b][n])*sgw[b][t];
	}else{
		sgwx[b][t+1] = (1-ks[b][n])*sgw[b][t]+MIN( qd[b][t]+qin[b][t] , perc[b][n] );
	}
	sgw[b][t+1] = MIN(sgwx[b][t+1] , sgwmax);
	return (1);
}
