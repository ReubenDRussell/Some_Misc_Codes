#include "TFitResult.h"
#include "TFitResultPtr.h"

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/// IDEAS FOR IMPROVEMENTS:
//	- Look into "residuals" --> point on curve at energy - point plotted at energy ==> can do for all and work out the standard deviation from this
//	- for S1576: Look into TIGRESS GEANT4 --> add source yeilds to allow for effective extrapolation to higher energies
//
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

Double_t PhotoEfficiency(Double_t* dim, Double_t* par){		// photoefficiency function fit from GRSISort [ https://github.com/GRIFFINCollaboration/GRSISort/blob/master/libraries/TAnalysis/TGRSIFit/TGRSIFunctions.cxx ]
	double sum = 0.0;
	for(int i = 0; i < 9; ++i) {
		sum += par[i] * pow(TMath::Log(dim[0]), i);
	}
	return TMath::Exp(sum);
}

//////////////////////////////////////////////////////////////////////////////////////

Double_t PhotoEfficiency_Radware(Double_t* dim, Double_t* par){		// photoefficiency function fit from JG3 paper (This is the Radware way of doing it I think)

	double a = par[0];
	double b = par[1];
	double c = par[2];
	double d = par[3];
	double e = par[4];
	double f = par[5];
	double g = par[6];

	double x = TMath::Log(dim[0]/100);
	double y = TMath::Log(dim[0]/1000);

	double f1 = a + (b*x) + (c*x*x);
	double f2 = d + (e*y) + (f*y*y);
	double f3 = pow(f1,-g) + pow(f2,-g);
	double f4 = pow(f3,-1/g);
	return TMath::Exp(f4);
}

//////////////////////////////////////////////////////////////////////////////////////

int line_counter(){
	unsigned int number_of_lines = 0;
	FILE *infile = fopen("RelEff_input.txt", "r");
	int ch;

	while (EOF != (ch=getc(infile))){	// while not the end of the file, get a new line in the infile
		if ('\n' == ch){
			++number_of_lines;
		}
	 }
	printf("%u\n", number_of_lines);
	return number_of_lines;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void EffFit(){	// just the standard GRSI photoefficiency fit to data
	int Num_Line = line_counter();
	Double_t Energy[Num_Line], EnergyErr[Num_Line], Eff[Num_Line], EffErr[Num_Line];

	string myText;
	std::stringstream ss;
	float en,enerr, eff, efferr;

	// Read from the text file
	ifstream MyReadFile("RelEff_input.txt");

	int i=0;
	// Use a while loop together with the getline() function to read the file line by line
	while (getline (MyReadFile, myText)) {
		// Energy.ResizeTo(i+1);
		// EnergyErr.ResizeTo(i+1);
		// Eff.ResizeTo(i+1);
		// EffErr.ResizeTo(i+1);

		// Output the text from the file
		ss << myText;
		ss >> en >> enerr >> eff >> efferr;
		ss.clear();
		enerr = 0.01;
		// Energy.push_back(en);
		// EnergyErr.push_back(enerr);
		// Eff.push_back(eff);
		// EffErr.push_back(efferr);
		Energy[i] = en;
		EnergyErr[i] = enerr;
		Eff[i] = eff;
		EffErr[i] = efferr;
		i++;
	} // end while loop


	// WHEN THE NUMBERS ARE NEAR 0, IT SEEMS TO MAKE THE FIT GO TOO FUCKY, THEREFORE THIS IS USEFUL:
	for(int j=0; j<Num_Line; j++){
		Eff[j] *= 10;		// Putting Efficiencies into Percentages (ish)
		EffErr[j] *= 10;		// Putting Efficiencies into Percentages (ish)
	}


	TGraphErrors	*gra	= new TGraphErrors(Num_Line, Energy, Eff, EnergyErr, EffErr);
	TCanvas *canvas = new TCanvas;
	gra->SetTitle("Relative Efficiency");
	gra->SetMarkerStyle(20);
	gra->GetYaxis()->SetRangeUser(0,11);
	gra->Draw("AP");
	gra->GetYaxis()->SetTitle("Relative Efficiency (arb.)");
	gra->GetXaxis()->SetTitle("Energy (keV)");

	TF1	*fit	= new TF1("myfunc",PhotoEfficiency,0,1500,9);

	fit->SetNpx(1000);	// number of points plotted (makes smooth curve/peak)

	gra->Fit(fit);
	fit->Draw("same");


	// saves the picture (you can save the root file separately you lazy shit)
	canvas->SaveAs("Relative_Efficiency_Curve.png");
	canvas->SaveAs("Relative_Efficiency_Curve.pdf");
	// canvas->SaveAs("Relative_Efficiency_Curve.root");

}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void EffFitRad(){	// just using Radware photoefficiency fit to data
	int Num_Line = line_counter();
	Double_t Energy[Num_Line], EnergyErr[Num_Line], Eff[Num_Line], EffErr[Num_Line];

	string myText;
	std::stringstream ss;
	float en,enerr, eff, efferr;

	// Read from the text file
	ifstream MyReadFile("RelEff_input.txt");

	int i=0;
	// Use a while loop together with the getline() function to read the file line by line
	while (getline (MyReadFile, myText)) {
		// Energy.ResizeTo(i+1);
		// EnergyErr.ResizeTo(i+1);
		// Eff.ResizeTo(i+1);
		// EffErr.ResizeTo(i+1);

		// Output the text from the file
		ss << myText;
		ss >> en >> enerr >> eff >> efferr;
		ss.clear();
		enerr = 0.01;
		// Energy.push_back(en);
		// EnergyErr.push_back(enerr);
		// Eff.push_back(eff);
		// EffErr.push_back(efferr);
		Energy[i] = en;
		EnergyErr[i] = enerr;
		Eff[i] = eff;
		EffErr[i] = efferr;
		i++;
	} // end while loop


	// WHEN THE NUMBERS ARE NEAR 0, IT SEEMS TO MAKE THE FIT GO TOO FUCKY, THEREFORE THIS IS USEFUL:
	for(int j=0; j<Num_Line; j++){
		Eff[j] *= 10;		// Putting Efficiencies into Percentages (ish)
		EffErr[j] *= 10;		// Putting Efficiencies into Percentages (ish)
	}


	TGraphErrors	*gra	= new TGraphErrors(Num_Line, Energy, Eff, EnergyErr, EffErr);
	TCanvas *canvas = new TCanvas;
	gra->SetTitle("Relative Efficiency");
	gra->SetMarkerStyle(20);
	gra->GetYaxis()->SetRangeUser(0,11);
	gra->Draw("AP");
	gra->GetYaxis()->SetTitle("Relative Efficiency (arb.)");
	gra->GetXaxis()->SetTitle("Energy (keV)");

	TF1	*fit	= new TF1("myfunc",PhotoEfficiency_Radware,0,1500,7);

	fit->SetNpx(1000);	// number of points plotted (makes smooth curve/peak)

	gra->Fit(fit);
	fit->Draw("same");

	// saves the picture (you can save the root file separately you lazy shit)
	canvas->SaveAs("Relative_Efficiency_Curve.png");
	canvas->SaveAs("Relative_Efficiency_Curve.pdf");
	// canvas->SaveAs("Relative_Efficiency_Curve.root");
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void EffFit_output(float Number){	// uses GRSI photoefficiency curve AND prints out the  value at the asked for Energy
	int Num_Line = line_counter();
	Double_t Energy[Num_Line], EnergyErr[Num_Line], Eff[Num_Line], EffErr[Num_Line];

	string myText;
	std::stringstream ss;
	float en,enerr, eff, efferr;

	// Read from the text file
	ifstream MyReadFile("RelEff_input.txt");

	int i=0;
	// Use a while loop together with the getline() function to read the file line by line
	while (getline (MyReadFile, myText)) {
		// Energy.ResizeTo(i+1);
		// EnergyErr.ResizeTo(i+1);
		// Eff.ResizeTo(i+1);
		// EffErr.ResizeTo(i+1);

		// Output the text from the file
		ss << myText;
		ss >> en >> enerr >> eff >> efferr;
		ss.clear();
		enerr = 0.01;
		// Energy.push_back(en);
		// EnergyErr.push_back(enerr);
		// Eff.push_back(eff);
		// EffErr.push_back(efferr);
		Energy[i] = en;
		EnergyErr[i] = enerr;
		Eff[i] = eff;
		EffErr[i] = efferr;
		i++;
	} // end while loop


	// WHEN THE NUMBERS ARE NEAR 0, IT SEEMS TO MAKE THE FIT GO TOO FUCKY, THEREFORE THIS IS USEFUL:
	for(int j=0; j<Num_Line; j++){
		Eff[j] *= 10;		// Putting Efficiencies into Percentages (ish)
		EffErr[j] *= 10;		// Putting Efficiencies into Percentages (ish)
	}


	TGraphErrors	*gra	= new TGraphErrors(Num_Line, Energy, Eff, EnergyErr, EffErr);
	TCanvas *canvas = new TCanvas;
	gra->SetTitle("Relative Efficiency");
	gra->SetMarkerStyle(20);
	gra->GetYaxis()->SetRangeUser(0,11);
	gra->Draw("AP");
	gra->GetYaxis()->SetTitle("Relative Efficiency (arb.)");
	gra->GetXaxis()->SetTitle("Energy (keV)");

	TF1	*fit	= new TF1("myfunc",PhotoEfficiency,0,1500,9);

	fit->SetNpx(1000);	// number of points plotted (makes smooth curve/peak)

	gra->Fit(fit);
	fit->Draw("same");

	// // saves the picture (you can save the root file separately you lazy shit)
	// canvas->SaveAs("Relative_Efficiency_Curve.png");
	// canvas->SaveAs("Relative_Efficiency_Curve.pdf");
	// // canvas->SaveAs("Relative_Efficiency_Curve.root");



	TFitResultPtr result = gra->Fit("myfunc","S");
	// double x[5] = { 778, 850, 356, 333, 521 };
	double x[1] = {Number};
	double err[1];  // error on the function at point x0
	result->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, false); // (number in array, step in coordinate space, step in dimension space, array of values in x to look at, confidence interval at x (array of errors to output), a number???, normalisation yes/no)
	for(int i=0; i<1; i++){	// i can be changed depending on what other numbers you are putting in??
		cout << " function value at " << x[i] << " = " << fit->Eval(x[i]) << " +/- " << err[i] << endl;
	}

}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

void EffFitRad_output(float Number){	// uses Radware photoefficiency curve AND prints out the value at the asked for Energy
		int Num_Line = line_counter();
	Double_t Energy[Num_Line], EnergyErr[Num_Line], Eff[Num_Line], EffErr[Num_Line];

	string myText;
	std::stringstream ss;
	float en,enerr, eff, efferr;

	// Read from the text file
	ifstream MyReadFile("RelEff_input.txt");

	int i=0;
	// Use a while loop together with the getline() function to read the file line by line
	while (getline (MyReadFile, myText)) {
		// Energy.ResizeTo(i+1);
		// EnergyErr.ResizeTo(i+1);
		// Eff.ResizeTo(i+1);
		// EffErr.ResizeTo(i+1);

		// Output the text from the file
		ss << myText;
		ss >> en >> enerr >> eff >> efferr;
		ss.clear();
		enerr = 0.01;
		// Energy.push_back(en);
		// EnergyErr.push_back(enerr);
		// Eff.push_back(eff);
		// EffErr.push_back(efferr);
		Energy[i] = en;
		EnergyErr[i] = enerr;
		Eff[i] = eff;
		EffErr[i] = efferr;
		i++;
	} // end while loop


	// WHEN THE NUMBERS ARE NEAR 0, IT SEEMS TO MAKE THE FIT GO TOO FUCKY, THEREFORE THIS IS USEFUL:
	for(int j=0; j<Num_Line; j++){
		Eff[j] *= 10;		// Putting Efficiencies into Percentages (ish)
		EffErr[j] *= 10;		// Putting Efficiencies into Percentages (ish)
	}


	TGraphErrors	*gra	= new TGraphErrors(Num_Line, Energy, Eff, EnergyErr, EffErr);
	TCanvas *canvas = new TCanvas;
	gra->SetTitle("Relative Efficiency");
	gra->SetMarkerStyle(20);
	gra->GetYaxis()->SetRangeUser(0,11);
	gra->Draw("AP");
	gra->GetYaxis()->SetTitle("Relative Efficiency (arb.)");
	gra->GetXaxis()->SetTitle("Energy (keV)");

	TF1	*fit	= new TF1("myfunc",PhotoEfficiency_Radware,0,1500,7);

	fit->SetNpx(1000);	// number of points plotted (makes smooth curve/peak)

	gra->Fit(fit);
	fit->Draw("same");

	// // saves the picture (you can save the root file separately you lazy shit)
	// canvas->SaveAs("Relative_Efficiency_Curve.png");
	// canvas->SaveAs("Relative_Efficiency_Curve.pdf");
	// // canvas->SaveAs("Relative_Efficiency_Curve.root");




	TFitResultPtr result = gra->Fit("myfunc","S");
	// double x[5] = { 778, 850, 356, 333, 521 };
	double x[1] = {Number};
	double err[1];  // error on the function at point x0
	result->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, false); // (number in array, step in coordinate space, step in dimension space, array of values in x to look at, confidence interval at x (array of errors to output), a number???, normalisation yes/no)
	for(int i=0; i<1; i++){	// i can be changed depending on what other numbers you are putting in??
		cout << " function value at " << x[i] << " = " << fit->Eval(x[i]) << " +/- " << err[i] << endl;
	}
}

void MASH_help(){	// the help command to tell you what you are doing

	std::cout << 
	"\nMASH is a scrambled together mix of the efficiency codes I have previously made.\n"
	"The code may be a bit janky in places but it should work.\n \n"
	"Use:\n"
	"\tMake an input file called \"RelEff_input.txt\"\n"
	"\tFill this input with:\n"
	"\t\tEnergy  Energy Error  \"Efficiency\"  \"Efficiency\" Error\n"
	"\tSeparate each with a tab (space might work, untested) and no blank lines.\n"
	"\tRun inside ROOT. Get output. Profit.\n\n"
	"The Before Times:\n"
	"\tMeasure the peak areas of the sources\n"
	"\tNote which peaks these are related to and assign an Intensity (see nndc decay charts)\n"
	"\tCalculate Area/(% Intensity/100) (same as Area/Norm_Intensity)\n"
	"\tIF you have the activities you can work out an absolute efficiency\n"
	"\tELSE: Assume this is an efficiency (roughly that is correct) you just need to scale the sources to one another\n"
	"\tYou should have a set of data points which can be plotted\n"
	"\tIf you want to make them look \"nice\" divide the numbers by \n\t a number a little higher than the largest (this assumes the largest is close to max efficiency),\n\t this will make them all under 1, this code has a portion that multiplies by 10 since there is a bug when <<0.\n\n"
	// "After Use:\n"
	// "\tTake the \"Efficiencies\" calculated at the desired points and multiply them \n with the numbers you have measured your peaks at, this is likely to make them wierd looking, but they are all scaled to be the same?"
	"The Functions:\n"
	"\t \"EffFit()\" = GRSI photoefficiency fit to the data (& saves picture).\n"
	"\t \"EffFit_output()\" = same as EffFit(), but outputs Eff at requested number.\n"
	"\t \"EffFitRad(number)\" = RadWare photoefficiency fit to data (SHIT/BROKEN) (& saves picture)\n"
	"\t \"EffFitRad_output(number)\" = same as EffFitRad(), but outputs Eff at requested number.\n\n"
	"I've tested it for most of the common sources used, need to check higher energies (e.g. near 4077 and above)\n"
	"If there are further questions, ask Radar\n \n"
	"That is all, kindly fuck off."
	<< std::endl;

}