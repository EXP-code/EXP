/*
  Recursive factorial computation
*/

/* Test routine
main()
{
	double		factrl();
	int		n;

	printf("n: ");
	scanf("%d",&n);
	printf("\n%d!=%f\n",n,factrl(n));
}
*/

double factrl(int n)
{
	if (n<=1)
		return 1;
	else
		return n*factrl(n-1);
}

