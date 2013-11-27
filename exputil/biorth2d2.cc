				// By recurrance relation

void CBDisk::potl_recur(int np, int m, double r, Vector& a)
{
  int n = np-1;
  double pfac = pow(r, (double)m+1.0e-20);
  
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

  a[1] = pfac*cur;

  if (n==0) return;
  
  double curl1 = cur;
  double curl2;
  
  fac *= r2 - 1.0;
  cur *= fac*(2*m+1);

  a[2] = pfac*cur;

  if (n==1) return;

  for (int nn=2; nn<=n; nn++) {
    curl2 = curl1;
    curl1 = cur;
    cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
      (1.0 + (double)(2*m-1)/nn)*curl2;
    a[nn+1] = pfac*cur;
  }
  
  return;
}



