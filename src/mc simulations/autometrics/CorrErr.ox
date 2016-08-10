#include <oxstd.oxh>


VarCoVar(const p, const v, const N, const n, const allcorr) //takes partial correlations, variance, and number of regressors, number of relevant variables, and whether all variables are correlated or not

{
   	decl s;	  //variance covariance matrix
	decl i, j, k;


	if (allcorr)
	{
	
	 s = constant(p, N, N);
	}

	else
	{
	s = zeros(N,N);

	for(k=0; k<n; k++)
	{

	for(j=0; j<n; j++)
	 {
	  s[k][j] = p;

	 }

	}
	}
	//set variance = 1
	for(k=0; k<N; k++)
	{
	  s[k][k] = 1;

	}


	decl corrmatrix = zeros(N,N);

	for(k=0; k<N; k++)
	{
	for(j=0; j<N; j++)
	 {
	 corrmatrix[k][j] = s[k][j]/(sqrt(s[k][k])*sqrt(s[j][j]));
	   

 
	}
	}
//	println("variance covariance matrix is", s);
//	println("corr matrix is", corrmatrix);
	return corrmatrix;
}

Choleski(const s)	 //takes variance covariance matrix
{
 decl c = choleski(s);
 //println("choleski decomp is", c);
 return c;

 
 
}




CorrErr(const N, const t, const c) //takes number of regressors and number of observations - right now for mu =0, and choleski decomposed matrix
{
	decl mu; //matrix of means
	decl m_z = zeros(t, N);

	mu = zeros(1, N);

	decl z, x;

	decl m;
	for(m=0; m<t; m++)

	{
	 z = rann(1,N);
	 x= mu + z*c';
	 m_z[m][] = x;

	}

	//m_z = rann(t,N);
	//println("m_z is", m_z);
	//println("M_z covariane structure is", variance(m_z));
	//	println("M_z correlation structure is", correlation(m_z));
	return m_z;
			   //this function creates the matrix of errors

}


CreateCorrErr(const p, const v, const numrel, const NumregTot, const T, const allcorr)	//takes error covariance, error variance, num relevant variables, total number of regressors, number of observations
 {
	 decl VcV = VarCoVar(p, v, NumregTot, numrel, allcorr) ;
	 //println(VcV);
	 decl C = Choleski(VcV);
	 //println(C);
	 decl Errors = CorrErr(NumregTot,T,C);
	 //println(Errors);
	 return Errors;
	 //println("Variance covariance matrix is", VcV);
}




/*main()
{
	//  decl errors = VarCoVar(0.99, 1, 80, 5);
	  //println("errors are", errors);
	//  decl chol = Choleski(errors);
	//  decl V = chol*chol';
	  //println("choleski", chol);
	// decl correlatederrors = CorrErr(80, 100, chol);
	// decl corrcalc = variance(correlatederrors);

	 decl errors = CreateCorrErr(0.9, 1, 5, 10, 20);
	decl correrrors = correlation(errors);
	println("correlation calculation is", correrrors);
;
	 
	 // print(errors);
	  //print(chol);
}

*/



