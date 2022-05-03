{
	#include "TStyle.h"

	TStyle *vgcstyle;
	//void setvgcstyle() {
	vgcstyle= new TStyle("vgcstyle","VGCStyle");

	//------ define color gradinet
	const Int_t NRGBs = 6;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.2, 0.4, 0.6, 0.8,1 };
    Double_t red[NRGBs]   = { 1, 0.70, 0.5,  0.45, 0.38 , 0.19};
    Double_t green[NRGBs] = { 1, 0.82, 0.61, 0.36, 0.16 , 0.06};
    Double_t blue[NRGBs]  = { 1, 0.9,  0.88, 0.74, 0.43 , 0.148};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    vgcstyle->SetNumberContours(NCont);
	//----------------------------------------------------------------------

	// canvas
	vgcstyle->SetCanvasBorderMode(0);
	vgcstyle->SetCanvasColor(kWhite);
	vgcstyle->SetCanvasDefH(700);//Height of canvas
	vgcstyle->SetCanvasDefW(900); //Width of canvas
	vgcstyle->SetCanvasDefX(0);   //POsition on screen
	vgcstyle->SetCanvasDefY(0);

	// pad
	vgcstyle->SetPadBorderMode(0);
	vgcstyle->SetPadColor(kWhite);
	vgcstyle->SetPadGridX(false);
	vgcstyle->SetPadGridY(false);
	vgcstyle->SetGridColor(0);
	vgcstyle->SetGridStyle(3);
	vgcstyle->SetGridWidth(1);

	// frame
	vgcstyle->SetFrameBorderMode(0);
	vgcstyle->SetFrameBorderSize(1);
	vgcstyle->SetFrameFillColor(0);
	vgcstyle->SetFrameFillStyle(0);
	vgcstyle->SetFrameLineColor(1);
	vgcstyle->SetFrameLineStyle(1);
	vgcstyle->SetFrameLineWidth(1);
  
	// histogram
	vgcstyle->SetHistLineColor(1);
	vgcstyle->SetHistLineStyle(0);
	vgcstyle->SetHistLineWidth(1);
	vgcstyle->SetEndErrorSize(2);
	vgcstyle->SetMarkerStyle(20);

	// fit/function:
	vgcstyle->SetOptFit(0);
	vgcstyle->SetFitFormat("5.4g");
	vgcstyle->SetFuncColor(2);
	vgcstyle->SetFuncStyle(1);
	vgcstyle->SetFuncWidth(1);

	// the date
	vgcstyle->SetOptDate(0);

	// statistics box
	vgcstyle->SetOptFile(0);
	vgcstyle->SetOptStat(0); // To display the mean and RMS: SetOptStat("mr");
	vgcstyle->SetStatColor(kWhite);
	vgcstyle->SetStatFont(42);
	vgcstyle->SetStatFontSize(0.025);
	vgcstyle->SetStatTextColor(1);
	vgcstyle->SetStatFormat("6.4g");
	vgcstyle->SetStatBorderSize(1);
	vgcstyle->SetStatH(0.1);
	vgcstyle->SetStatW(0.15);
	//vgcstyle->SetStatStyle(Style_t style = 1001);
	//vgcstyle->SetStatX(Float_t x = 0);
	//vgcstyle->SetStatY(Float_t y = 0);

	// margins
	vgcstyle->SetPadTopMargin(0.05);
	vgcstyle->SetPadBottomMargin(0.11);
	vgcstyle->SetPadLeftMargin(0.11);
	vgcstyle->SetPadRightMargin(0.05);

	// global title
	vgcstyle->SetOptTitle(0);
	vgcstyle->SetTitleFont(42);
	vgcstyle->SetTitleColor(1);
	vgcstyle->SetTitleTextColor(1);
	vgcstyle->SetTitleFillColor(10);
	vgcstyle->SetTitleFontSize(0.05);
	//vgcstyle->SetTitleH(0); // Set the height of the title box
	//vgcstyle->SetTitleW(0); // Set the width of the title box
	//vgcstyle->SetTitleX(0); // Set the position of the title box
	//vgcstyle->SetTitleY(0.985); // Set the position of the title box
	//vgcstyle->SetTitleStyle(Style_t style = 1001);
	// vgcstyle->SetTitleBorderSize(2);

	// axis titles
	vgcstyle->SetTitleColor(1, "XYZ");
	vgcstyle->SetTitleFont(132, "XYZ");
	vgcstyle->SetTitleSize(0.04, "XYZ");
	//vgcstyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	//vgcstyle->SetTitleYSize(Float_t size = 0.02);
	vgcstyle->SetTitleXOffset(0.9);
	vgcstyle->SetTitleYOffset(1.25);
	//vgcstyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	// axis labels
	vgcstyle->SetLabelColor(1, "XYZ");
	vgcstyle->SetLabelFont(132, "XYZ");
	vgcstyle->SetLabelOffset(0.007, "XYZ");
	vgcstyle->SetLabelSize(0.03, "XYZ");

	// axis
	vgcstyle->SetAxisColor(1, "XYZ");
	vgcstyle->SetStripDecimals(kTRUE);
	vgcstyle->SetTickLength(0.03, "XYZ");
	vgcstyle->SetNdivisions(510, "XYZ");
	vgcstyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	vgcstyle->SetPadTickY(1);

	// log plots (default off)
  	vgcstyle->SetOptLogx(0);
  	vgcstyle->SetOptLogy(0);
  	vgcstyle->SetOptLogz(0);

	// postscript options:
	vgcstyle->SetPaperSize(20.,20.);
	// vgcstyle->SetLineScalePS(Float_t scale = 3);
	// vgcstyle->SetLineStyleString(Int_t i, const char* text);
	// vgcstyle->SetHeaderPS(const char* header);
	// vgcstyle->SetTitlePS(const char* pstitle);
	// vgcstyle->SetBarOffset(Float_t baroff = 0.5);
	// vgcstyle->SetBarWidth(Float_t barwidth = 0.5);
	// vgcstyle->SetPaintTextFormat(const char* format = "g");
	// vgcstyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	// vgcstyle->SetTimeOffset(Double_t toffset);
	// vgcstyle->SetHistMinimumZero(kTRUE);
	vgcstyle->SetHatchesLineWidth(5);
	vgcstyle->SetHatchesSpacing(0.05);

	// legend
	vgcstyle->SetLegendFont(42);
	vgcstyle->SetLegendTextSize(0.03);
	vgcstyle->SetLegendBorderSize(0);
	vgcstyle->SetLegendFillColor(10);
	//void 	SetLegendBorderSize (Width_t size=4)
	//void 	SetLegendFillColor (Color_t color=0)
	//void 	SetLegendFont (Style_t font=62)
	//void 	SetLegendTextSize (Double_t size=0.)


	// colors
	// TRIUMF primary colour palette (+ black && white)
	TColor *tr_cyan      = new TColor(2011,58/255.0, 193/255.0, 227/255.0); // Pantone Process Cyan
	// TRIUMF secondary colour palette
	TColor *tr_cool_gray_10 = new TColor(2012,128/255.0, 130/255.0, 133/255.0);
	TColor *tr_cool_gray_2  = new TColor(2013,209/255.0, 211/255.0, 212/255.0);
	TColor *tr_cool_gray_1  = new TColor(2014,241/255.0, 242/255.0, 242/255.0);
	TColor *tr_yellow    = new TColor(2015,225/255.0, 203/255.0, 5/255.0); // Pantone P 7-8 C
	// TRIUMF tertiary colour palette
	TColor *tr_green     = new TColor(2016,117/255.0, 192/255.0, 67/255.0); //Pantone 376 C
	TColor *tr_navy      = new TColor(2017,35/255.0, 35/255.0, 89/255.0); // Pantone P 105-8 C
	TColor *tr_salmon    = new TColor(2018,243/255.0, 113/255.0, 94/255.0); //Pantone P 55-6 C
	TColor *tr_orange    = new TColor(2019,217/255.0, 132/255.0, 76/255.0); //Pantone P 55-6 C
	// TRIUMF old color palette
	TColor *tr_blue      = new TColor(2020,0.067,0.247,0.486);
	TColor *tr_wine      = new TColor(2021,0.455,0.141,0.122);
	TColor *tr_green2    = new TColor(2022,0.314,0.502,0.255);
	TColor *tr_orange2   = new TColor(2023,0.839,0.38,0.18);
	TColor *tr_gray1     = new TColor(2024,0.353,0.314,0.318);
	TColor *tr_blue2     = new TColor(2025,0.121,0.285,0.488);


	//vgcstyle->SetCanvasPreferGL(kFALSE);
	//vgcstyle->cd();

	// automatic load
	gROOT->SetStyle("vgcstyle"); //uncomment to set this style                                 

}


