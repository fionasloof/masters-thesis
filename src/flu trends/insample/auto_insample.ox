#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <modelbase>
#import <simulator>
#import <packages/PcGive/pcgive_ects>

main()
{
	decl i, j, k, l;
	decl T = 313; // There are 261 observations in sub sample - the 262th is the first to be nowcast. But indexing starts at 0 so T=261 is the first to be nowcast
	decl c_alpha = 0.01;
	decl numforecasts = 52;//a years worth of forecasts
	decl m_nowcasts = zeros(numforecasts, 2);
	
	decl fludata =	new PcGive();
	fludata.Load("/Users/f/projects/masters-thesis/src/flu trends/InSample/correlate-ILI-FullSample_lags.csv");			 //load data set that has been transformed to include all lags
	decl names = fludata.GetAllNames();

	fludata.Select(Y_VAR, {"ILI", 0, 0});

	for (i=2; i<columns(names); i++) //start at i=2 beause i=0 is date and i=1 is FluData (i.e. y variable)
	{
	   fludata.Select(X_VAR, {names[i], 0, 0});	 //third value is number of lags
	}

	fludata.SetSelSample(1, 1, T, 1); //estimate using up to observation T
	fludata.Autometrics(c_alpha, "none", 1);
	fludata.AutometricsSet("pvalue_tests", 0);
	fludata.Estimate();

	decl vcoefs_est = fludata.GetPar();
	decl vcoefs_names = fludata.GetParNames();
	decl vcoefs_stderr = fludata.GetStdErr();
	decl vcoefs_tstat = zeros(rows(names),1);

	println(vcoefs_est);
	println(vcoefs_names);
	println(vcoefs_tstat);
	//want a matrix that has all the search terms coefficients i.e. coeffcients = 0 for those not chosen
	
	decl coefmatrix = zeros(rows(names), 4);
	decl vgum_sel = strfind(names, vcoefs_names); 	//returns a matrix of the indexes where the ascoefs_est are found in asx_gum
	coefmatrix[vgum_sel][1] = vcoefs_est;
	coefmatrix[vgum_sel][2] = vcoefs_tstat;

	decl observations = fludata.GetAll();

	for (k=0; k<numforecasts; k++)
   	{
	decl nowcastobs = observations[261+k][];
	decl nowcast_calc = nowcastobs.*coefmatrix[][1]';
	decl nowcast = sumr(nowcast_calc);
	
	println(nowcastobs);
	println("now cast is", nowcast);

	m_nowcasts[k][0] =	nowcastobs[1];
	m_nowcasts[k][1] = nowcast; 	//First column has true values, second column has forecasts 
	 }

	 println(m_nowcasts);

	 
	decl se = (m_nowcasts[][0] - m_nowcasts[][1]).^2;
	decl mse = meanc(se);
	println("mse is", mse);
	
}
