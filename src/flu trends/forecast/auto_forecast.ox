#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <modelbase>
#import <simulator>
#import <packages/PcGive/pcgive_ects>

// FIRST TRANSFORM THE DATA USING Datatransform.ox
// saved spreadsheet is _lags

main()
{
	decl i, k;
	decl T = 261; 				// there are 261 observations in sub sample - the 262th is the first to be nowcast. But indexing starts at 0 so T=261 is the first to be nowcast
	decl p_a = 0.2;
	decl numforecasts = 52;		// a years worth of forecasts
	decl m_nowcasts = zeros(numforecasts, 2);

	decl fludata =	new PcGive();
	fludata.Load("/Users/f/projects/masters-thesis/src/flu trends/forecast/correlate-ILI_1_forecast_Ox_lags.csv");	 //load data set that has been transformed to include all lags

	decl names = fludata.GetAllNames();

	fludata.Select(Y_VAR, {"ILI", 0, 0});
	for (i=2; i<columns(names); i++) // start at i=2 beause i=0 is date and i=1 is FluData (i.e. y variable)
	{
	 fludata.Select(X_VAR, {names[i], 0, 0});	 // third value is number of lags
	}

 //Select initial model

	fludata.SetSelSample(1, 1, 261, 1); 
	fludata.Autometrics(p_a, "none", 1);
	fludata.AutometricsSet("pvalue_tests", 0);
	fludata.Estimate();

	decl vcoefs_est = fludata.GetPar();
	decl vcoefs_names = fludata.GetParNames();


	//want a matrix that has all the search terms coefficients i.e. coeffcients = 0 for those not chosen
	
	decl coefmatrix = zeros(rows(names), 1);
	decl vgum_sel = strfind(names, vcoefs_names); //returns a matrix of the indexes where the ascoefs_est are found in asx_gum

	coefmatrix[vgum_sel][] = vcoefs_est;
	decl observations = fludata.GetAll();

	decl nowcastobs = observations[T][];
	decl nowcast_calc = nowcastobs.*coefmatrix';
	decl nowcast = sumr(nowcast_calc);

	println(nowcastobs);
	println("now cast is", nowcast);

	m_nowcasts[0][0] =	nowcastobs[1];
	m_nowcasts[0][1] = nowcast;	// First column has true values, second column has forecasts

	//Next step: only use the selected variables and re-estimate using OLS for each new observations, finding new fitted value

	decl selectedvar = fludata.GetParNames();
	println(selectedvar);


	fludata.DeSelect();
	fludata.Select(Y_VAR, {"ILI", 0, 0});
	fludata.Select(U_VAR, {"Constant", 0, 0});

		for (i=0; i<rows(selectedvar); i++) //start at i=2 beause i=0 is date and i=1 is FluData (i.e. y variable)
		{
	    fludata.Select(X_VAR, {selectedvar[i], 0, 0});	 //third value is number of lags
		}

	decl resids ;
		for (k=1;k<52;k++)
		{

		fludata.SetMethod(M_OLS);
		fludata.SetSelSample(1,1,T+k,1);
		fludata.Estimate();
		resids = fludata.GetResiduals();
		nowcastobs = observations[T+k][];
		decl nowcast = nowcastobs[1]+resids[rows(resids)-1];
		m_nowcasts[k][0] =	nowcastobs[1];
		m_nowcasts[k][1] = nowcast;	
    	}

	
	decl se = (m_nowcasts[][0] - m_nowcasts[][1]).^2;
	decl mse = meanc(se);
	println("mse is", mse);



	decl results = new Database();

	results.Append(m_nowcasts, {"True", "Fitted"});

}
