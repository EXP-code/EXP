extern int ncom;	/* defined in LINMIN */
extern double *pcom,*xicom,(*nrfunc)();

double f1dim(double x)
{
	int j;
	double f,*xt,*dvector();
	void free_dvector();

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_dvector(xt,1,ncom);
	return f;
}
