#include <oxstd.oxh>
#import <modelbase>
#import <simulator>

#include <oxprob.h>
#include <oxdraw.h>
#import <ranmc>
#import <packages/PcGive/pcgive_ects>
#include "CorrErr.ox"
#include "CaseDescriptions.ox"

main()
{
	decl orthogonalregressors = TRUE;
	decl test = FALSE;
	decl errcov, errvar;
	decl T, cRep, L, p_a;
	decl i, j, k, m, f, h, g, Xeps, auto, N, q, vt, psi;
	decl time = timer();

	decl casenum =  1;
	
	decl w;
	T = 100;									// Sample size (currently fixed for fixed regs)
	h = 20;										// Number of regressors
	N = 5;										// Number of relevant regressors
	cRep = 1000;								// Number of replications
	p_a = 0.01;					  				// Reduction significance level
	errcov = 0.9;
	decl allcorr = FALSE;						// Covariance between regressors
	errvar = 1;

	if (casenum == 7 ||	casenum == 8  || casenum == 9  || casenum == 10 || casenum == 11 || casenum == 12 || casenum == 19 || casenum == 20
	|| casenum == 21 || casenum == 22 || casenum == 23 || casenum == 24 || casenum == 25 || casenum == 26 || casenum == 27 || casenum == 28
	|| casenum == 50 || casenum == 51 || casenum == 52 || casenum == 53 || casenum == 54 || casenum == 55 || casenum == 56 || casenum == 57)
	{
	 orthogonalregressors = FALSE;
	}

	if (casenum == 2 ||	casenum == 4  || casenum == 6  || casenum == 8  || casenum == 10 || casenum == 12 || casenum == 14 || casenum == 16
	|| casenum == 18 || casenum == 20 || casenum == 22 || casenum == 24 || casenum == 51 || casenum == 53 || casenum == 55 || casenum==57 )
	{
	 L = 120;
	}
	else
	{
	 L= 80;
	}

	if (casenum == 25 || casenum == 26 || casenum == 27 || casenum == 28 ||	casenum == 29 || casenum == 30 || casenum == 54
	||  casenum == 55 || casenum == 56 || casenum == 57)
	{
	 allcorr = TRUE;
	}

//Variance of regressors
	println("This is case number ", casenum);
	println("In this simulatation...");
	println("Observations: ", T);
	println("Number of regressors: ", L);
	println("Number of relevant regressors: ", N);
	println("Number of MC repetitions: ", cRep);
	println("The significance level is: ", p_a);
	println("The variance of each of the regressors is ", errvar);
	if (orthogonalregressors)
	{
	println("The regressors are all orthogonal");
	}

	
	if (!orthogonalregressors)
	{
	println("The regressors are not orthogonal");
	
	if (!allcorr)
	{
	 println("There is only correlation between the relevant variables");	
	 println("The covariance between each of the relevant variables is ", errcov);
	}
	 else
	 {
	  println("There is correlation between all regressors");
	  println("The covariance between each of the regressors is ", errcov);
	 }
	}
	
	decl cx_gum, asx_gum = { "Constant" }, mx_gum, y, eps, mx_gumh;
	asx_gum ~= { "Ly" };
	for (k = 1; k <= L; ++k)
	{
		asx_gum ~= { sprint("x", k) };
	}

	 cx_gum = sizeof(asx_gum);
	 mx_gum = constant(.NaN, T, L+2);

	decl rho_c = 0.5;
	decl rho = rho_c*ones(L,1), gamma = 0.5;
	decl x = zeros(T+h,L), yh = zeros(T+h,1);
	decl v, p;
	decl coeff;
	decl c_alpha = quant(1-(p_a/2), T);
	decl beta0 = <1>;
	decl beta = Case(casenum, L, N);

//create indicator variables to indicate which variables are included in model
	decl mTruePar =  beta0|gamma|beta;
	decl mTrueVar = zeros(rows(mTruePar), 1);
	decl mmse_mc = zeros(cx_gum, cRep);
	decl mestimatedcoeffs_mc = zeros(cx_gum, cRep);
	decl mtstat_mc = zeros(cx_gum, cRep);
	decl mestcoefssols_mc = zeros(N+1, cRep);
	decl mstderrols_mc = zeros(N+1, cRep);
	decl mtstatsols_mc = zeros(N+1, cRep);
	decl minclols_mc = zeros(N+1, cRep);
 	decl msquarederr_mc	= zeros(N+1, cRep);
	decl tbar_mc = zeros(cx_gum, cRep);
	decl noncentrals_mc = zeros(cx_gum, cRep);

//a separate matrix not including constant (used in calculations)
	decl truepar = mTruePar[1:6];

	for(m=0; m<L; m++)
	{
	if (mTruePar[m][0]!=0)	
	 mTrueVar[m][0] = 1;	 
	}

//Calculate noncentralities
	decl mpsi = zeros(N+1, 1);
	decl varx = 1;
	for (p=0; p<N+1; p++)
	{
	 mpsi[p] = truepar[p]/sqrt(1/(T*varx));
	}

//Calculate theoretical retention probabilities
	decl retprob = zeros(N+1, 1);
	for (p=0; p<N+1; p++)
	{
	 retprob[p] = 1-probt((c_alpha-fabs(mpsi[p])), T);
	}
	
	println("\nThe true parameters of the model are: ");
	for (i=0; i<cx_gum; i++)
	{
	 if (mTruePar[i] != 0)
	 {
	 println(asx_gum[i], "=", mTruePar[i]);
	 }
	}

	decl TrueParValues = new Database();
	TrueParValues.Append(mTruePar, "TrueParameterValues");
	TrueParValues.SaveCsv(sprint("/Users/f/projects/masters-thesis/data/simulated data/Case", casenum, "/TrueParValues"));

//Create database for modelling
	decl model = new PcGive();
	model.Create(1/*freq*/, 1, 1, 1, T);		
	model.Append(constant(.NaN, T, 2), {"y","eps"});
	model.Append(mx_gum, asx_gum[ : L+1]);
	model.Deterministic(-1);
	model.SetPrint(FALSE);

//Set up matrix design for evaluation
	decl mindresults = zeros(cx_gum, cRep);	//indicator results
	decl aspreadsheeti; //name of spreadsheet exported
	decl mindvar_mc = zeros(cx_gum, 1);
	decl m_retentionrate = zeros(cx_gum, 1);

	if (orthogonalregressors)
	{	
	v = rann(T+h,L);
	vt = v*sqrt((1-rho_c^2));
	for (p = 1; p < T+h; p++)	  			// the possible regressors stay the same for each of the replications
		{									//the only thing that changes are eps!	which then changes the y variables
			x[p][] = x[p-1][].*0.5 + vt[p][];
		}
	}
	else
	{
	v = CreateCorrErr(errcov, errvar,N, L, T+h, allcorr);
	for (p = 1; p < T+h; p++)	  			// the possible regressors stay the same for each of the replications
		{									//the only thing that changes are eps!	which then changes the y variables
			x = v;
		}
	}

	if (test)
	{
	decl tempx = x[][:4];
	for (p=0; p<N; p++)
	{
	 x[][p] = tempx[][4-p];
	}
	}

	for (i = 0; i < cRep; i++)	  			// cRep simulations within each loop
		{
		mindvar_mc = zeros(cx_gum, 1);

//SIMULATION MODEL DESIGN
		eps = rann(T+h,1);	// N[0,1]
		
		for (p = 1; p < T+h; p++)	  			// cRep simulations within each loop
		{
			yh[0][] = 0;
			yh[p][] = beta0 + gamma*yh[p-1][] + x[p][]*beta + eps[p][];	 	// generate the dependent variable
		}
		
		y = yh[h:][];
		mx_gumh = ones(T+h,1)~lag0(yh,1)~x;
		mx_gum = mx_gumh[h:][];
		aspreadsheeti = sprint("repnum", i+1, ".csv");
		decl repi = new Database();
		repi.Append(y, "y");
		repi.Append(mx_gum, asx_gum);
		repi.SaveCsv(sprint("/Users/f/projects/masters-thesis/data/simulated data/Case", casenum, "\RepNum",i+1));
		
			model.Renew(y, "y");  						// Renews the dependent variable for each draw
			model.Renew(mx_gum, asx_gum[0 : L+1]);		// Regressors may be unchanged			
			model.DeSelect();
			model.Select(Y_VAR, {"y", 0, 0});
			model.Select(U_VAR, {"Constant", 0, 0});
			model.Select(X_VAR, {"Ly", 0, 0});

			for (k = 1; k <= L; ++k)					// Replace L with j+1 and remove autometrics selection for DGP
			{
				model.Select(X_VAR, {sprint("x", k), 0, 0});
			}

			model.SetSelSample(1, 1, T, 1);
			model.SetMethod(M_OLS);
			model.Autometrics(p_a, "none", 1);
			model.AutometricsSet("pvalue_tests", 0);
			model.Estimate();

			decl vcoefs_num = model.GetParCount();
			decl vcoefs_est = model.GetPar();
			decl vstderr_est = model.GetStdErr();
			decl ascoefs_est = model.GetParNames();
			decl vgum_sel = strfind(asx_gum, ascoefs_est);  // creturns a matrix of the index where the ascoefs_est are found in asx_gum
			mindvar_mc[vgum_sel][0] = 1;					// assigns a value of 1 to the index's found in vgum_sel - i.e. indicator matrix 
			mestimatedcoeffs_mc[vgum_sel][i] = vcoefs_est;  // this "aligns" the coefficients with their names in asx_gum				
			mtstat_mc[vgum_sel][i] = vcoefs_est ./ vstderr_est;
			mtstat_mc[vgum_sel][i] = vcoefs_est ./ vstderr_est;
			mindresults[][i] = mindvar_mc;
			tbar_mc[][i]= mtstat_mc[][i] -(densn(c_alpha - mtstat_mc[][i]) - densn(-c_alpha - mtstat_mc[][i]))/(1-probn(c_alpha-mtstat_mc[][i])+probn(-c_alpha-mtstat_mc[][i]));
			noncentrals_mc[][i] = mtstat_mc[][i] - (densn(c_alpha - tbar_mc[][i])-densn(-c_alpha - tbar_mc[][i]))/(1-probn(c_alpha - tbar_mc[][i])+probn(-c_alpha-tbar_mc[][i]));
			mmse_mc[][i] = (mTruePar - mestimatedcoeffs_mc[][i]).^2;



//BEGINS MCMC SIMULATIONS WITHOUT VARIABLE SELECTION
			model.DeSelect();
			model.Select(Y_VAR, {"y", 0, 0});
			model.Select(U_VAR, {"Constant", 0, 0});
			model.Select(X_VAR, {"Ly", 0, 0});

			for (k = 1; k <= N; ++k)			//Replace L with j+1 and remove autometrics selection for DGP
			{
				model.Select(X_VAR, {sprint("x", k), 0, 0});
			}

			model.SetSelSample(1, 1, T, 1);
			model.SetMethod(M_OLS);
			model.Estimate();
			decl parvals = model.GetPar();
			parvals = dropr(parvals, 6);
			decl parnames = model.GetParNames();
			parnames = dropr(parnames, 6);
			decl parstderr = model.GetStdErr();
			parstderr = dropr(parstderr, 6);
			decl tstat = parvals ./ parstderr;
			decl isparsign = zeros(N+1, 1);
			decl selectedpars = zeros(N+1, 1);


//Parameter estimates are given estimates if parameter was chosen, otherwise = 0;
			for (p=0; p<N+1; p++)
			{
			if (fabs(tstat[p]) > c_alpha)
			{
			 selectedpars[p] = parvals[p]; 
			}
			}
			decl squarederr_UMSE = (selectedpars-truepar).^2;
			
			for (p=0; p<N+1; p++)
			{
			 if (fabs(tstat[p])>c_alpha)
			 {
			 isparsign[p] = 1;
			 }
			}

			mestcoefssols_mc[][i] = parvals;
			minclols_mc[][i] = isparsign;
			msquarederr_mc[][i] = squarederr_UMSE;
			
	}

	
	
//SUMMARY RESULTS FOR AUTOMETRICS	

	//Unconditional MSE
		decl mUMSE = zeros(cx_gum, 1);
		mUMSE = meanr(mmse_mc);
		decl RUSMSE = sqrt(mUMSE);

	//Conditional MSE
		decl mCMSE = zeros(cx_gum, 1);		
		mCMSE = sumr(mmse_mc .* (mestimatedcoeffs_mc .!= 0)) ./ sumr(mestimatedcoeffs_mc .!= 0);
		decl RCSMSE = sqrt(mCMSE);
	
		m_retentionrate = meanr(mindresults);
		decl m_potcalc = constant(2,1,1);
		decl m_gaugecalc = constant(2,1,1);

	//Statistic for MSE for retained irrelvant variables
		decl m_retratmse = zeros(cx_gum, 1);
		for (p=1;p<cx_gum; p++)
		{
		 if(mTruePar[p] == 0)
		 {
		  if(sumr(mindresults[p][])>0)
		  {
		  m_retratmse[p] = sumr(mmse_mc[p][])/sumr(mindresults[p][]);
		  }
		  else
		  {
		   m_retratmse[p] = 0;
		  }
		 }		 
		}
		decl sum_retratemse = sumc(m_retratmse);
		decl TMSE;
		
		for(j=2; j<cx_gum; j++)
		{
		 if (mTrueVar[j][0] == 1)

			 m_potcalc~=m_retentionrate[j][0];
		 else	
			m_gaugecalc~=m_retentionrate[j][0];		 
		}
		
		 decl m_potency, m_gauge;
		 m_potency = m_potcalc[1:];
		 m_gauge = m_gaugecalc[1:];

		decl format = m_retentionrate ~ RUSMSE ~ RCSMSE	  ;
		print( "%r", asx_gum, "%c", {"Ret Rates", "     RUSMSE", "     RCSMSE"}, format) ;		
		println("\nPotency: ", meanr(m_potency));
		println("\nGauge: ", meanr(m_gauge));

//SUMMARY RESULTS WITH SELECTION FROM LDGP
  
	// Unconditional MSE
		decl mUMSE_OLS = zeros(N+1, 1);
		mUMSE_OLS = meanr(msquarederr_mc);
		decl RUIMSE = sqrt(mUMSE_OLS);

	// Conditional MSE
		decl mCMSE_OLS = zeros(N+1, 1);
		mCMSE_OLS = sumr(msquarederr_mc .* (minclols_mc .!= 0)) ./ sumr(minclols_mc .!= 0);
		decl RCIMSE = sqrt(mCMSE_OLS);
		decl m_retentionrate_OLS = zeros(N+1,1);
		m_retentionrate_OLS = meanr(minclols_mc);
		decl m_potcalc_OLS = constant(2,1,1);
		m_potcalc_OLS = meanc(m_retentionrate_OLS[1:N]);
		decl format2 = 	m_retentionrate_OLS~RUIMSE~RCIMSE;

		println("\n\n\n\nResults from selection from LDGP:");		
		print( "%r", asx_gum[1:], "%c", {"Ret Rates", "     RUIMSE", "     RCIMSE"}, format2 ) ;
		println("\nPotency:", m_potcalc_OLS);

		savemat("ExportTest.xlsx", format2~m_potcalc_OLS,{"Ret Rates", "     RUIMSE", "     RCIMSE"} );
	
		decl meannoncentrals = meanr(noncentrals_mc);	
		println("\n\n\n");
		println("\nfinished in ", timespan(time));	
}
