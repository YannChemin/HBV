#include "arrays.h"
// open output files
int hbv_report(int btot,float *mtot,int ttot,float **qo,float **qr,float **po,float **eta,float **eto,float **mfc,float **mbeta,float **mlp,float **malpha,float **mkf,float **mks,float **mperc,float **mcflux,float **MNSET,float **MRVEET,float **MNS,float **MNSH,float **MNSL,float **MRVE,float **MRMAEH,float **MRMAEL,float **MREVE){
	int i,t,q;
	FILE *f1[btot];
	f1[0] = fopen ("Output01-Doab.csv", "w");			// output stochastic Doab
	f1[1] = fopen ("Output02-Pole_Chehr.csv", "w");			// output stochastic Pole Chehr
	f1[2] = fopen ("Output03-Doabe_M.csv", "w");			// output stochastic Doabe M
	f1[3] = fopen ("Output04-Ghor_B.csv", "w");			// output stochastic Ghor B
	f1[4] = fopen ("Output05-Holilan.csv", "w");			// output stochastic Holilan
	f1[5] = fopen ("Output06-Pole_D.csv", "w");			// output stochastic Pole D
	f1[6] = fopen ("Output07-Jelogir.csv", "w");			// output stochastic Jelogir
	f1[7] = fopen ("Output08-Paye_P.csv", "w");			// output stochastic Paye P
	FILE *f2[btot];
	f2[0] = fopen ("Basinout01det-Doab.csv", "w");			// output deterministic Doab
	f2[1] = fopen ("Basinout02det-Pole_Chehr.csv", "w");		// output deterministic Pole Chehr
	f2[2] = fopen ("Basinout03det-Doabe_M.csv", "w");		// output deterministic Doabe M
	f2[3] = fopen ("Basinout04det-Ghor_B.csv", "w");		// output deterministic Ghor B
	f2[4] = fopen ("Basinout05det-Holilan.csv", "w");		// output deterministic Holilan
	f2[5] = fopen ("Basinout06det-Pole_D.csv", "w");		// output deterministic Pole D
	f2[6] = fopen ("Basinout07det-Jelogir.csv", "w");		// output deterministic Jelogir
	f2[7] = fopen ("Basinout08det-Paye_P.csv", "w");		// output deterministic Paye P
	FILE *f3[btot];
	f3[0] = fopen ("ETout01et-Doab.csv", "w");			// output ET Doab
	f3[1] = fopen ("ETout02et-Pole_Chehr.csv", "w");		// output ET Pole Chehr
	f3[2] = fopen ("ETout03et-Doabe_M.csv", "w");			// output ET Doabe M
	f3[3] = fopen ("ETout04et-Ghor_B.csv", "w");			// output ET Ghor B
	f3[4] = fopen ("ETout05et-Holilan.csv", "w");			// output ET Holilan
	f3[5] = fopen ("ETout06et-Pole_D.csv", "w");			// output ET Pole D
	f3[6] = fopen ("ETout07et-Jelogir.csv", "w");			// output ET Jelogir
	f3[7] = fopen ("ETout08et-Paye_P.csv", "w");			// output ET Paye P
	FILE *f4[btot];
	f4[0] = fopen ("EToutput01etsto-Doab.csv", "w");		// output ETstochastic Doab
	f4[1] = fopen ("EToutput02etsto-Pole_Chehr.csv", "w");		// output ETstochastic Pole Chehr
	f4[2] = fopen ("EToutput03etsto-Doabe_M.csv", "w");		// output ETstochastic Doabe M
	f4[3] = fopen ("EToutput04etsto-Ghor_B.csv", "w");		// output ETstochastic Ghor B
	f4[4] = fopen ("EToutput05etsto-Holilan.csv", "w");		// output ETstochastic Holilan
	f4[5] = fopen ("EToutput06etsto-Pole_D.csv", "w");		// output ETstochastic Pole D
	f4[6] = fopen ("EToutput07etsto-Jelogir.csv", "w");		// output ETstochastic Jelogir
	f4[7] = fopen ("EToutput08etsto-Paye_P.csv", "w");		// output ETstochastic Paye P

	//PRINT TO FILE
	for (i=0;i<btot;i++){//Sub Basin loop
		for (t=0;t<ttot;t++){
			fprintf(f2[i],"%f %f %f \n",qo[i][t],qr[i][t],po[i][t]);
			fprintf(f3[i],"%f %f \n",eta[i][t],eto[i][t]);
		}
		//Write sub basins output datasets
		fprintf(f1[i],"fc(mm),beta(-),lp(-),alpha(-),kf(1/d),ks(1/d),perc(mm/d),cflux(mm/d),NS(-),NSH(-),NSL(-),RVE(pc),RMAEH(-),RMAEL(-),REVE(pc)\n");
		fprintf(f4[i],"fc(mm),beta(-),lp(-),alpha(-),kf(1/d),ks(1/d),perc(mm/d),cflux(mm/d),NSET(-),RVEET(pc)\n");
		for (q=0;q<(int)mtot[1];q++){
			fprintf(f1[i],"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",mfc[i][q],mbeta[i][q],mlp[i][q],malpha[i][q],mkf[i][q],mks[i][q],mperc[i][q],mcflux[i][q],MNS[i][q],MNSH[i][q],MNSL[i][q],MRVE[i][q],MRMAEH[i][q],MRMAEL[i][q],MREVE[i][q]);
			fprintf(f4[i],"%f %f %f %f %f %f %f %f %f %f\n",mfc[i][q],mbeta[i][q],mlp[i][q],malpha[i][q],mkf[i][q],mks[i][q],mperc[i][q],mcflux[i][q],MNSET[i][q],MRVEET[i][q]);
		}
	}
	//CLOSE FILES
	for (i=0;i<btot;i++){//Sub Basin loop
		fclose(f1[i]);
		fclose(f2[i]);
		fclose(f3[i]);
		fclose(f4[i]);
	}
	return (1);
}
