#include <oxstd.oxh>


Case(const casenum, const L, const N) //takes partial correlations, variance, and number of regressors

{
   	decl i;	  //variance covariance matrix

	decl beta = zeros(L,1);
	
	if(casenum==1 || casenum==2 || casenum== 7|| casenum==8)
	{
	
	for (i=0; i<N; ++i)
	{
	 beta[i][0] = 0.6;
	}

	 }


	 
	if(casenum==3 || casenum==4 || casenum== 9|| casenum==10)
	{
	
	for (i=0; i<N; ++i)
	{
	 beta[i][0] = 0.2;
	}

	 }


	 if(casenum==5 || casenum==6 || casenum== 11|| casenum==12 || casenum==25 || casenum==26)
	{
	
	for (i=0; i<N; ++i)
	{
	 beta[i][0] = (i+2)*0.1;
	}

	 }

	 if(casenum==13 || casenum==14 || casenum==19 || casenum==20)
	 {
	 for (i=0; i<N; ++i)
	{
	 beta[i][0] = 0.6*(-1)^(i+1);
	}
	}


	if(casenum==15 || casenum==16 || casenum==21 || casenum==22)
	 {
	 for (i=0; i<N; ++i)
	{
	 beta[i][0] = 0.2*(-1)^(i+1);
	}
	}


	if(casenum==17 || casenum==18 || casenum==23 || casenum==24 || casenum==27|| casenum==28||casenum==60)
	 {
	 for (i=0; i<N; ++i)
	{
	 beta[i][0] = (i+2)*0.1*(-1)^(i+1);
	}
	}

	if(casenum==50 ||casenum==51)
	{

	
	  for (i=0; i<N; ++i)
	{
	 beta[i][0] = (i+1)*0.5;
	}

	beta[0][0] = 0.575;
	beta[1][0] = 0.84;
	beta[2][0] = 1.25;
	beta[3][0] = 1.55;
	beta[4][0] = 1.75;
	//beta[5][0] = 2.1;
	}


		if(casenum==52||casenum==53)
	{

	
	  for (i=0; i<N; ++i)
	{
	 beta[i][0] = (i+1)*0.5;
	}

	beta[0][0] = -0.575;
	beta[1][0] = 0.84;
	beta[2][0] = -1.25;
	beta[3][0] = 1.55;
	beta[4][0] = -1.75;
	//beta[5][0] = 2.1;
	}


	if(casenum==54||casenum==55)
	{

	
	  for (i=0; i<N; ++i)
	{
	 beta[i][0] = (i+1)*0.5;
	}

	beta[0][0] = 0.4;
	beta[1][0] = 0.6;
	beta[2][0] = 0.9;
	beta[3][0] = 1.1;
	beta[4][0] = 1.25;
	//beta[5][0] = 2.1;
	}

		if(casenum==56||casenum==57)
	{

	
	  for (i=0; i<N; ++i)
	{
	 beta[i][0] = (i+1)*0.5;
	}

	beta[0][0] = -0.4;
	beta[1][0] = 0.6;
	beta[2][0] = -0.9;
	beta[3][0] = 1.1;
	beta[4][0] = -1.25;
	//beta[5][0] = 2.1;
	}


	  
	


	


	return beta;
}



/* main()
{
	  decl L = 80;
	  
	  decl beta = Case(1, 80, 10);
	  println(beta);
}  */


