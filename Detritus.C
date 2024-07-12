// Detritus Test Code - A crushed up mix of different peak fitting and background subtraction codes
// Mix includes: PeakFitDemo.C | PeakFitDemoRR.C | Bridge.C | CalibFitter.C | PeakFind_one.C | PeakFind.C | GateSigBG_bckgnd.C | GateSigBG_sum2.C | GateSigBG.C | New.C


/*
 *	OPEN ROOT WITH THE .root FILE YOU WANT TO USE - IT THEN IS CALLED _file0 AUTOMATICALLY (as seen later on)
 *
 *	(remember to have made histogram "h" before using PeakFit)
 *	> .L PeakFitDemo.C;
 *	> PeakFit(h);
 *
 *	
 *	remember that ROOT numbers bins from 1 (not 0)
 */


/*	What this code needs to do:
 |		- Have ability to exclude an area (the peak itself)
 |		- Fit the background (with the peak excluded)
 |		- "Fix"/set the fit parameters for fitting the peak to those of the background fit (leaving parameter 0 to vary)
 |		- Fit the desired peak
 |
 |	Inputs needed:
 |		- Upper and Lower bound of the background range
 |		- Centroid OR Upper and Lower bound of centroid location may be better
 |		- Sigma? - maybe not needed
 |
*/

// Setting a Boolian here for use later.
bool reject_peak;

//-------------------------------------------------------------------------------------------------------------
// The Doubles:
//-------------------------------------------------------------------------------------------------------------
double	Step(double *dim, double *par){

	double	x	= dim[0];
	double	area	= par[0];
	double	cent	= par[1];	
	double	sigm	= par[2];
	double	step	= par[3];

	return	abs(step) * area/(sigm*TMath::Sqrt(2*TMath::Pi())) * TMath::Erfc((x - cent) / (sigm * TMath::Sqrt(2)));

}
//-------------------------------------------------------------------------------------------------------------
double 	Peak_step(double *dim, double *par){	// peak with step

	double	x	= dim[0];
	double	area	= par[0];
	double	cent	= par[1];
	double	sigm	= par[2];

	double	step	= 0;
	if(par[3])
		step	= Step(dim,par);	// sending stuff to the Step function

	return	area * (1./(sigm*TMath::Sqrt(2*TMath::Pi()))) * TMath::Gaus(x,cent,sigm)	+ step; // Step commented out as it is not needed at the moment

}
//-------------------------------------------------------------------------------------------------------------
double 	Peak(double *dim, double *par){		// peak without step

	double	x	= dim[0];
	double	area	= par[0];
	double	cent	= par[1];
	double	sigm	= par[2];

	double	step	= 0; // no step for this one

	return	area * (1./(sigm*TMath::Sqrt(2*TMath::Pi()))) * TMath::Gaus(x,cent,sigm)	;//+ step; // Step commented out as it is not needed at the moment

}
//-------------------------------------------------------------------------------------------------------------
double	BG(double *dim, double *par){

	double	x	= dim[0];
	double	p0	= par[0];
	double	p1	= par[1];
	double	p2	= par[2];

	return	x * x * p2 + x * p1 + p0;	// quadratic equation for fit the background


}
//-------------------------------------------------------------------------------------------------------------
double	BG2(double *dim, double *par){

	double	x	= dim[0];
	double	p0	= par[0];
	double	p1	= par[1];
	double	p2	= par[2];
	double	cent	= par[3];

	// The below commented out bit works BUT it doesn't look too nice
	// IF I HAVE DONE THIS RIGHT, THE REGION ABOUT THE PEAK SHOULD BE AVOIDED FOR THE FITTING? - yes, but the fit has a chunk going to 0 for this range
	if(reject_peak && x > (cent-10) && x < (cent+10)){	// +/- 10 bins or units of energy, not sure, from the peak.
		TF1::RejectPoint();	// 
		return 0;
	}

	return	x * x * p2 + x * p1 + p0;	// quadratic equation for fit the background

}
//-------------------------------------------------------------------------------------------------------------
double	SinglePeak(double *dim, double *par){

	double	x	= dim[0];	// variable x (iterates over, not a parameter per se)?
	double	area	= par[0];	// area of peak
	double	cent	= par[1];	// centroid of peak
	double	sigm	= par[2];	// sigma (crosssection or +/- where 60 ish% of counts lie) of peak
	double	step	= par[3];	// step in bavkground (allowing for a different height in background either side of the peak)
	double	p0	= par[4];	//	}
	double	p1	= par[5];	//	 }-- Parameters from the background fit
	double	p2	= par[6];	//	}

	double	*par_peak	= new double[4];
	double	*par_bg		= new double[3];

	par_peak[0]	= area;		//	Parameters going to the Peak function 
	par_peak[1]	= cent;		//
	par_peak[2]	= sigm;		//
	par_peak[3]	= step;		//

	par_bg[0]	= p0;		//	Parameters going to the BG function
	par_bg[1]	= p1;		//
	par_bg[2]	= p2;		//

	double	result	= 0;	
	result 	+= Peak(dim,par_peak);
	result	+= BG(dim,par_bg);

	return	result;

	// // If not cutting up into different functions it looks like this: (also used in GaussPeakFit.C & CalibFitter.C & PeakFind_one.C)
	// // BUT THE ONE BELOW HAS A Linear BACKGROUND (first line) I THINK
	// return	m * x	+ c	
	// 		+	abs(step) * area/(sigm*TMath::Sqrt(2*TMath::Pi())) * TMath::Erfc((x - cent) / (sigm * TMath::Sqrt(2)))
	// 		+	area * (1./(sigm*TMath::Sqrt(2*TMath::Pi()))) * TMath::Gaus(x,cent,sigm);  


}
//-------------------------------------------------------------------------------------------------------------
double	SinglePeak_step(double *dim, double *par){

	double	x	= dim[0];	// variable x (iterates over, not a parameter per se)?
	double	area	= par[0];	// area of peak
	double	cent	= par[1];	// centroid of peak
	double	sigm	= par[2];	// sigma (crosssection or +/- where 60 ish% of counts lie) of peak
	double	step	= par[3];	// step in bavkground (allowing for a different height in background either side of the peak)
	double	p0	= par[4];	//	}
	double	p1	= par[5];	//	 }-- Parameters from the background fit
	double	p2	= par[6];	//	}

	double	*par_peak	= new double[4];
	double	*par_bg		= new double[3];

	par_peak[0]	= area;		//	Parameters going to the Peak function 
	par_peak[1]	= cent;		//
	par_peak[2]	= sigm;		//
	par_peak[3]	= step;		//

	par_bg[0]	= p0;		//	Parameters going to the BG function
	par_bg[1]	= p1;		//
	par_bg[2]	= p2;		//

	double	result	= 0;	
	result 	+= Peak_step(dim,par_peak);
	result	+= BG(dim,par_bg);

	return	result;

}


//-------------------------------------------------------------------------------------------------------------
// The Voids:
//-------------------------------------------------------------------------------------------------------------

// void PeakFit(TH1D *h, double cent, double sigm, double Low, double High){
void PeakFit(TH1D *h, double cent, double Low, double High){


	// for a proper exclusion zone, both the root forums and the jroot code says to set bin content and error to 0??? 
	// will this fuck up things if I re-draw the histogram? Or potentially make it impossible to fir the peak after fitting the background?
	// -- For now I'll just make the background much larger than normal to ensure that the peak counts are negligible
	
	reject_peak = true;	// boolian to reject the peak from the background fit
	TF1	*fBack	= new TF1("Background_Fit_Test", BG2, Low, High, 4);
	fBack->FixParameter(3,cent);	// fixing a set parameter to the value asked [we want this one to be a "constant" for this background]

	// fBack->SetParameters(3,0,cent);
	// fBack->SetParLimits(0,0,1e4);
	// fBack->SetParLimits(2,cent - 5, cent + 5);





	// Below makes a fit background of above and below the peak separately. -- Works
	/*
	TF1	*fBack1	= new TF1("Background_Fit_Test", BG2, Low, cent-10, 4);		//low BG
	fBack1->FixParameter(3,cent);
	TF1	*fBack2	= new TF1("Background_Fit_Test", BG2, cent+10, High, 4);	//high BG
	fBack2->FixParameter(3,cent);

	// fit 1 and 2 here
	fBack1->SetLineColor(kBlack);
	fBack2->SetLineColor(kBlack);
	h->Fit(fBack1,"RB+");
	h->Fit(fBack2,"RB+");
	*/
	



	// Below doesn't work, dont try it (unless you can fix it)
	/*double AvPar0;
	double AvPar1;
	double AvPar2;
	AvPar0 = (fBack1->GetParameter(0) + fBack2->GetParameter(0));
	AvPar1 = (fBack1->GetParameter(1) + fBack2->GetParameter(1));
	AvPar2 = (fBack1->GetParameter(2) + fBack2->GetParameter(2));

	TF1	*fBack12 = new TF1("Background Fit Test", BG,Low,High,3);		//BG of avg low and high
	fBack12->FixParameter(0,AvPar0); // do same with param errors????
	fBack12->FixParameter(1,AvPar1);
	fBack12->FixParameter(1,AvPar2);
	fBack12->SetLineColor(kBlue);
	h->Fit(fBack12,"RB+");*/



	// Messing about with cloning the histogram then zero weighting the peak bins
	/*
	int *centroid_bin;

	TH1D *h_clone = (TH1D*)h->Clone();	// a cone of histogram "h" to zero weight the bins at the peak
	for(int bin=1; bin<(h_clone->GetNBinsX()) ;bin++){
		centroid_bin = int(cent);	// only meant as a rough check
		if(bin<(cent-10)&&bin>(cent+10)){
			h_clone->SetBinError(bin,0);
		}
	}
	*/





	fBack->SetLineColor(kGreen);
	h->Fit(fBack,"RB+");
	double Par0 = fBack->GetParameter(0);
	double Par1 = fBack->GetParameter(1);
	double Par2 = fBack->GetParameter(2);
	
	
	


	TF1	*fProj	= new TF1("Fit",SinglePeak,Low,High,7);		// ("name", equation of fit, lower bound, upper bound, number of parameters)
	// fProj->SetParameters(100,cent,sigm,0,0,0,0);	// setting parameters start point	(area, centroid, sigma, step, BG par,BG par,BG par)
	fProj->SetParameters(100,cent,1,0,0,0,0);
	// fProj->SetParLimits(0,0,1e4);	// setting upper and lower bound for parameter 0 [area?] (parameter to constrain, lower bound, upper bound)
	fProj->SetParLimits(0,0,1e6);	// setting upper and lower bound for parameter 0 [area?] (parameter to constrain, lower bound, upper bound)
	fProj->SetParLimits(1,cent - 5, cent + 5); // Only allowing the centroid parameter to vary by +/- 5 keV
	fProj->SetParLimits(2,1,10); // limits on sigma, just in case

	// fProj->FixParameter(1,cent);	// fixing a set parameter to the value asked
	// fProj->FixParameter(2,sigm);	// -------------"---------------

	// new TCanvas	// may be useful to add if fitting multiple histograms or looping?
	
	// fixing the parameters for the background fit based on the fit of the background not using the peak counts
	// fProj->FixParameter(4,Par0);	// this is height of the background therefore need to allow it to vary	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//BUT IT LOOKS LIKE LETTING Par0 VARY RESULTS IN A PEAK AT MAX SIGMA IN SOME CASES
	fProj->FixParameter(5,Par1);
	fProj->FixParameter(6,Par2);


	fProj->SetLineColor(kRed);	// it is red by default, just doing this for completion
	h->Fit(fProj,"RB+");
	
	

}

//-------------------------------------------------------------------------------------------------------------

void MakeProj(char const *inhistname, char const *name, double lowbin, double highbin){	// only for use with TH2F/D

	TH2F* h2=(TH2F*)_file0->Get(inhistname);	// THIS ONE WORKS
	TH1D* h1 =(TH1D*) h2->ProjectionY(name,lowbin,highbin);	// makes an y projection of hist between lowbin and highbin

	h1->Draw();
}

//-------------------------------------------------------------------------------------------------------------

void MakeProjHard(char const *name, double lowbin, double highbin){	// only for use with TH2F/D	- changed to TH2D since it was having issues for some reason with Float not Double
	// Got fed up with error messages, so just add the path to the histogram here:
	// TH2D* h2=(TH2D*)_file0->Get("New_Useful_Sr80/AB_dopEn_v_Ring_Comp_Supp_620-700_dopGate_Sr80");	// Comment out these lines if you want to define your own h2 before using this function
	TH2F* h2=(TH2F*)_file0->Get("New_Useful_Pb208/AB_recEn_v_Ring_Comp_Supp_630-700_recGate_Pb208");
	// TH2F* h2=(TH2F*)_file0->Get("New_Useful_UpSt/AB_dopEn_v_Ring_Comp_Supp_620-700_dopGate_UpSt");

	TH1D* h1 =(TH1D*) h2->ProjectionY(name,lowbin,highbin);	// makes an y projection of hist between lowbin and highbin

	h1->Draw();
}

//-------------------------------------------------------------------------------------------------------------

void JustFit(TH1D *h, double cent, double Low, double High){

	// // Commented out block to show what the background fitted on its own looks like
	// reject_peak = true;	// boolian to reject the peak from the background fit
	// TF1	*fBack	= new TF1("Background_Fit_Test", BG2, Low, High, 4);
	// fBack->FixParameter(3,cent);	// fixing a set parameter to the value asked [we want this one to be a "constant" for this background]

	// fBack->SetLineColor(kGreen);
	// h->Fit(fBack,"RB+");

	TF1	*fProj	= new TF1("Fit",SinglePeak,Low,High,7);		// ("name", equation of fit, lower bound, upper bound, number of parameters)
	// fProj->SetParameters(100,cent,sigm,0,0,0,0);	// setting parameters start point	(area, centroid, sigma, step, BG par,BG par,BG par)
	fProj->SetParameters(100,cent,1,0,0,0,0);
	fProj->SetParLimits(0,0,1e4);	// setting upper and lower bound for parameter 0 [area?] (parameter to constrain, lower bound, upper bound)
	fProj->SetParLimits(1,cent - 5, cent + 5); // Only allowing the centroid parameter to vary by +/- 5 keV
	fProj->SetParLimits(2,1,10); // limits on sigma, just in case

	// fProj->FixParameter(1,cent);	// fixing a set parameter to the value asked
	// fProj->FixParameter(2,sigm);	// -------------"---------------

	// new TCanvas	// may be useful to add if fitting multiple histograms or looping?
	

	fProj->SetLineColor(kRed);	// it is red by default, just doing this for completion
	h->Fit(fProj,"RB+");

}

//-------------------------------------------------------------------------------------------------------------

void JustFitStep(TH1D *h, double cent, double Low, double High){

	// // Commented out block to show what the background fitted on its own looks like
	// reject_peak = true;	// boolian to reject the peak from the background fit
	// TF1	*fBack	= new TF1("Background_Fit_Test", BG2, Low, High, 4);
	// fBack->FixParameter(3,cent);	// fixing a set parameter to the value asked [we want this one to be a "constant" for this background]

	// fBack->SetLineColor(kGreen);
	// h->Fit(fBack,"RB+");

	TF1	*fProj	= new TF1("Fit",SinglePeak_step,Low,High,7);		// ("name", equation of fit, lower bound, upper bound, number of parameters)
	// fProj->SetParameters(100,cent,sigm,0,0,0,0);	// setting parameters start point	(area, centroid, sigma, step, BG par,BG par,BG par)
	fProj->SetParameters(100,cent,1,0,0,0,0);
	fProj->SetParLimits(0,0,1e4);	// setting upper and lower bound for parameter 0 [area?] (parameter to constrain, lower bound, upper bound)
	fProj->SetParLimits(1,cent - 5, cent + 5); // Only allowing the centroid parameter to vary by +/- 5 keV
	fProj->SetParLimits(2,1,6); // limits on sigma, just in case

	// fProj->FixParameter(1,cent);	// fixing a set parameter to the value asked
	// fProj->FixParameter(2,sigm);	// -------------"---------------

	// new TCanvas	// may be useful to add if fitting multiple histograms or looping?
	

	fProj->SetLineColor(kRed);	// it is red by default, just doing this for completion
	h->Fit(fProj,"RB+");

}

//-------------------------------------------------------------------------------------------------------------

void PeakFitFix(TH1D *h, double cent, double sigm, double Low, double High){

	// // Commented out block to show what the background fitted on its own looks like
	// reject_peak = true;	// boolian to reject the peak from the background fit
	// TF1	*fBack	= new TF1("Background_Fit_Test", BG2, Low, High, 4);
	// fBack->FixParameter(3,cent);	// fixing a set parameter to the value asked [we want this one to be a "constant" for this background]

	// fBack->SetLineColor(kGreen);
	// h->Fit(fBack,"RB+");

	TF1	*fProj	= new TF1("Fit",SinglePeak_step,Low,High,7);		// ("name", equation of fit, lower bound, upper bound, number of parameters)
	// fProj->SetParameters(100,cent,sigm,0,0,0,0);	// setting parameters start point	(area, centroid, sigma, step, BG par,BG par,BG par)
	fProj->SetParameters(100,cent,1,0,0,0,0);
	fProj->SetParLimits(0,0,1e4);	// setting upper and lower bound for parameter 0 [area?] (parameter to constrain, lower bound, upper bound)
	// fProj->SetParLimits(1,cent - 5, cent + 5); // Only allowing the centroid parameter to vary by +/- 5 keV
	// fProj->SetParLimits(2,1,6); // limits on sigma, just in case

	fProj->FixParameter(1,cent);	// fixing a set parameter to the value asked
	fProj->FixParameter(2,sigm);	// -------------"---------------

	// new TCanvas	// may be useful to add if fitting multiple histograms or looping?
	

	fProj->SetLineColor(kRed);	// it is red by default, just doing this for completion
	h->Fit(fProj,"RB+");

}

//-------------------------------------------------------------------------------------------------------------


// Adding in a "Help" thing since I will forget how to operate this after some time.
void DetritusHelp(){

	std::cout << 
	"\nDetritus.C is a crushed up mix of random codes to make a block of something that works. \n \n"
	"As long as you fit one peak at a time you can fit One, Two, Three, Many or even Lots! All on one histogram (unproven) \n \n"
	"How to work this code: \n \n"
	"Main Functions: \n"
	"\t MakeProj(name of hist you want, name of projection from the TH2D, Low Bin, High Bin) \n \n"
	"\t MakeProjHard(name of projection from the TH2D, Low Bin, High Bin) \n \n"
	// "\t PeakFit( \" name of TH1D \" , centroid, sigma, lower bound for bg fit, upper bound for bg fit) \n \n"
	"\t PeakFit( \" name of TH1D \" , centroid, lower bound for bg fit, upper bound for bg fit) \n \n"
	"\t JustFit( \" name of TH1D \" , centroid, lower bound for bg fit, upper bound for bg fit) \n \n"
	"\t JustFitStep( \" name of TH1D \" , centroid, lower bound for bg fit, upper bound for bg fit) \n \n"
	"\t PeakFitFix( \" name of TH1D \" , centroid (fixed), sigma(fixed), lower bound for bg fit, upper bound for bg fit) \n \n"
	"What each of the functions inside does: \n "
	"\t MakeProj = Makes a 1D projection of a 2D histogram \n"
	"\t MakeProjHard = Same as Makeproj(), but you need to hard code the histogram path \n"
	"\t PeakFit = Fits the Peak at the centroid you define (varies by +/- 5keV),\n \t\t modified background does separate fit without peak and has no step and allowes curvature \n"
	"\t JustFit = Just fits the peak (similar to PeakFit but without a modified background) \n"
	"\t JustFitStep = Same as JustFit but has a step included in the background \n"
	"\t PeakFitFix = If you want to fix the centroid and sigma of a peak \n \n"
	"If you are doing editing, most of the code has comments. Here are the doubles containing formulae defined at the top \n"
	"\t Step = Defines the step allowed in the background either side of the peak \n"
	"\t Peak_step = Equation of the peak allowing for the step \n"
	"\t Peak = Equation of the peak \n"
	"\t BG = Equation of the background (quadratic - if you want, change to y=mx+c style for linear, but this should adapt) \n"
	"\t BG2 = Same as BG but with a chunk cut out relating to the region around the peak itself \n"
	"\t SinglePeak = Full equation of the fit, Peak + Background \n"
	"\t SinglePeak_step = Same as SinglePeak but with a step \n \n"
	"Dont forget that you need to load the .root file when opening root (done by typing \" root name.root \") \n"
	"(This doesn't really matter to mention here as you can't load it up properly without doing this) \n \n"
	"If you are me and using it for the 80Sr, the three TH2D histograms you need are: \n"
	"\t \" New_Useful_Sr80/AB_dopEn_v_Ring_Comp_Supp_620-700_dopGate_Sr80 \" \n"
	"\t \" New_Useful_Pb208/AB_recEn_v_Ring_Comp_Supp_630-700_recGate_Pb208 \" \n"
	"\t \" New_Useful_UpSt/AB_dopEn_v_Ring_Comp_Supp_620-700_dopGate_UpSt \" \n \n"
	"If that is all, kindly fuck off."
	<< std::endl;
}