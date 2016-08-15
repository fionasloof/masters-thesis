#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <modelbase>
#import <simulator>
#import <packages/PcGive/pcgive_ects>

// In sample uses all 313 observations
// But because of lags can only estimate the model actually using observations 53 onways
// There are 261 observations in sub sample - the 262th is the first to be nowcast. But indexing starts at 0 so T=261 is the first to be nowcast

main()
{
	decl i, j, k, l;
	decl T = 313;	
	decl fludata =	new PcGive();
	decl p_a = 0.01;
	decl c_alpha = quant(1-(p_a/2), T);
	decl numforecasts = 52 ;					// number of forecasts minus the lags
	decl m_nowcasts = zeros(numforecasts, 2);
	fludata.Load("/Users/f/projects/masters-thesis/src/flu trends/InSample/correlate-ILI-FullSample_lags.csv");
	decl names = fludata.GetAllNames();

	
	//Variables to be included in GUM
	fludata.Select(Y_VAR, {"ILI", 0, 0});

	for (i=2; i<columns(names); i++) 			// start at i=2 beause i=0 is date and i=1 is FluData (i.e. y variable)
	{
	   fludata.Select(X_VAR, {names[i], 0, 0});	// third value is number of lags

	}

	fludata.SetSelSample(1, 54, T, 1); 			// estimate using up to observation T
	fludata.Autometrics(p_a, "none", 1);
	fludata.AutometricsSet("pvalue_tests", 0);
	fludata.Estimate();

	decl vcoefs_est = fludata.GetPar();
	decl vcoefs_names = fludata.GetParNames();
	decl vcoefs_stderr = fludata.GetStdErr();
	decl vcoefs_tstat = zeros(rows(names),1);

	decl coefmatrix = zeros(rows(names), 5); 	// this matrix will contain the original estimates, stderrs, t-stats, and bias corrected values	
	decl vgum_sel = strfind(names, vcoefs_names); // returns a matrix of the indexes where the ascoefs_est are found in asx_gum
	coefmatrix[vgum_sel][0] = vcoefs_est;
	coefmatrix[vgum_sel][1] = vcoefs_stderr;
	decl beta_1step = zeros(rows(names), 1)	;
	decl beta_2step	= zeros(rows(names),1) ;

	for (j=0; j<rows(names); j++)
	{
	 	if	 (coefmatrix[j][0] == 0)
		{
		coefmatrix[j][2] = 0;
		coefmatrix[j][3] = 0;		 
		}
		else
		{
	 	coefmatrix[j][2] = coefmatrix[j][0] / coefmatrix[j][1]; 
		decl db = coefmatrix[j][0];
		decl dt = coefmatrix[j][2];
		decl dr = (densn(c_alpha-dt)-densn(-c_alpha-dt)) / (1-probn(c_alpha-dt)+probn(-c_alpha-dt));
		decl dtbar = dt - dr;
		decl drbar = (densn(c_alpha-dtbar)-densn(-c_alpha-dtbar)) / (1-probn(c_alpha-dtbar)+probn(-c_alpha-dtbar));
		coefmatrix[j][3] = fabs(dt).> c_alpha .? db .* (1 - (dr ./ dt)) .: 0;
		coefmatrix[j][4] = fabs(dt).> c_alpha .? db .* (1 - (drbar ./ dt)) .: 0;
   	  	}
	}

	  decl observations = fludata.GetAll();

	for (k=0; k<52; k++)
   	{
	decl nowcastobs = observations[261+k][];
	decl nowcast_calc = nowcastobs.*coefmatrix[][4]';
	decl nowcast = sumr(nowcast_calc);
	m_nowcasts[k][0] =	nowcastobs[1];
	m_nowcasts[k][1] = nowcast;
	}
	
	println(coefmatrix);
	println(m_nowcasts);
	decl se = (m_nowcasts[][0] - m_nowcasts[][1]).^2;
	decl mse = meanc(se);
	println("mse is", mse);
	decl results = new Database();
	results.Append(m_nowcasts, {"True", "Fitted"});
	 
}
