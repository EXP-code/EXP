double gammp(double a, double x)
{
  double gamser,gammcf,gln;
  void gser(double *gamser, double a, double x, double *gln),gcf(double *gammcf, double a, double x, double *gln),nrerror();

  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMP");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }
}

