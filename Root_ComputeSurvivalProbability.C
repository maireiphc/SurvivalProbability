#if !defined(__CINT__) || defined(__MAKECINT__)

#include "Riostream.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TPaveLabel.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "TColor.h"
#include "TASImage.h"


#endif

// Antonin Maire, 
// Origin : PhD thesis, 2012
//  See CDS reference : CERN-THESIS-2011-263
//   https://cds.cern.ch/record/1490315/, section V.B-5, p.143-145
//    and especially Fig. V.7
// Last modification : 28 mai 2015
// 
// NB : the macro works only if compiled with Root (because of std::vector...) = root -l  ....C++g


// To be launched with ROOT : root -l Root_ComputeSurvivalProbability.C++g

// Int_t Root_ComputeSurvivalProbability(  TString Str_Part1ToDisplay = "kLambdaFromXi", // kLambdaFromOmega or kLambdaFromXi, kLambda, kK0s, kXi, kOmega, kD0, kDplus, kDSplus, kLambdaCplus
//                                         TString Str_Part2ToDisplay = "kLambdaFromOmega",
//                                         Int_t rWrite = 0
//                                     ){



Color_t gColor_1 = kRed+1   ;
Color_t gColor_2 = kOrange+1;
Color_t gColor_3 = kAzure+2 ;
Color_t gColor_4 = kViolet-5;


enum gkColor{
    black = 0,
    red,
    orange,
    yellow,
    green,
    cyan,
    azure,
    blue,
    violet,
    magenta,
    nbColor
};


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,     kRed-10,     kOrange-9,   kYellow-7,     kGreen-8,    kCyan-8,    kAzure-8,    kBlue-9,     kViolet-9,    kMagenta-9,    }; // for syst bands
const Int_t colors[]     = {kBlack,      kRed+1 ,     kOrange+1,   kYellow+2,     kGreen+3,    kCyan+2,    kAzure+2,    kBlue+1,     kViolet-5,    kMagenta+1,    };
const Int_t markers[]    = {kFullCircle, kFullCircle, kOpenSquare, kFullDotMedium,kOpenSquare, kFullCross, kFullCircle, kFullCircle, kOpenCircle,  kOpenDiamond   };



// Masses in Gev/c²
Double_t mOmega         = 1.67245 ;
Double_t mXi            = 1.32171 ;
Double_t mLambda        = 1.115683;
Double_t mK0s           = 0.497614;
Double_t mD0            = 1.86486 ;   
Double_t mDplus         = 1.86962 ;
Double_t mDSplus        = 1.96849 ;
Double_t mLambdaCplus   = 2.28646 ;

// mean lifetime in m
Double_t cTauOmega          = 0.02461;
Double_t cTauXi             = 0.0491 ;
Double_t cTauLambda         = 0.0789 ;
Double_t cTauK0s            = 0.026844;
Double_t cTauD0             = 122.9e-6;
Double_t cTauDplus          = 311.8e-6;
Double_t cTauDSplus         = 149.9e-6;
Double_t cTauLambdaCplus    = 59.9e-6;


enum gkDecayType{
    kLambdaFromXi = 0,
    kLambdaFromOmega,
    kLambda,
    kK0s,
    kXi,
    kOmega,
    kD0,
    kDplus,
    kDSplus,
    kLambdaCplus,
    nPartDecayType
};






void myLegendSetUp(TLegend *currentLegend, float currentTextSize);

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, 
                float currentBottom);

void myGraphSetUp(TGraph *currentGraph, Float_t currentMarkerSize, 
                  int currentMarkerStyle, int currentMarkerColor, 
                  int currentLineStyle,   int currentLineColor, int currentLineWidth);

void myOptions(Int_t lStat);


Double_t dProba_dlxi(Double_t *x, Double_t *par){
    
    
    Double_t Lo  = par[0];
    
    Double_t pXi     = par[1];  // momentum of the mother particle : Omega, Xi, Lambda
    Double_t pOmega  = par[1];  // depending on what you like to consider
    Double_t pPart   = par[1];  // Momentum of the single-decay particle (primary-particle case)
    
    Short_t lDecayType = par[2]; // kLambdaFromXi, kLambdaFromOmega, kLambda, ...
    
    Double_t Lxi            = x[0];
    Double_t Lomega         = x[0];
    Double_t lDecayLength   = x[0];
    
    
    // NOTE 1
    // Pour les cas, kLambdaFromXi et kLambdaFromOmega,
    // voir les hypothèses du chap. 5 thèse, 
    // http://cdsweb.cern.ch/record/1490315
    // NOTE 2
    // Numériquement, attention : il ne faut pas oublier le facteur 1/c qui traine.
    // --> tau devient c.Tau, dans la formule de probabilité de survie du PDG, 
    
    Double_t output = -1;
    
    switch(lDecayType){
        case kLambdaFromXi :        // calcul de l'intégrande
            output = TMath::Exp(-mXi * Lxi / (cTauXi*pXi) ) * mXi/(cTauXi*pXi) * TMath::Exp( -mLambda *(Lo-Lxi) / (cTauLambda * 0.85 * pXi) );
            break;
        case kLambdaFromOmega :     // calcul de l'intégrande
            output = TMath::Exp(-mOmega * Lomega / (cTauOmega*pOmega) ) * mOmega/(cTauOmega*pOmega) * TMath::Exp( -mLambda * (Lo-Lomega) / (cTauLambda * 0.66 * pOmega) );
            break;
        case kLambda :
            output = TMath::Exp(-mLambda * lDecayLength / (cTauLambda * pPart) );  
            break;
        case kK0s :
            output = TMath::Exp(-mK0s * lDecayLength / (cTauK0s * pPart) );  // FIXME
            break;
        case kXi :
            output = TMath::Exp(-mXi * lDecayLength / (cTauXi * pPart) );  // FIXME
            break;            
        case kOmega :
            output = TMath::Exp(-mOmega * lDecayLength / (cTauOmega * pPart) );  // FIXME
            break;
        case kD0 :
            output = TMath::Exp(-mD0 * lDecayLength / (cTauD0 * pPart) );  
            break;
        case kDplus :
            output = TMath::Exp(-mDplus * lDecayLength / (cTauDplus * pPart) );  
            break;
        case kDSplus :
            output = TMath::Exp(-mDSplus * lDecayLength / (cTauDSplus * pPart) );  
            break;
        case kLambdaCplus :
            output = TMath::Exp(-mLambdaCplus * lDecayLength / (cTauLambdaCplus * pPart) );  
            break;            
               
    } // end switch 
    
    return output;
  
    
}

Int_t ComputeProbability(TF1 *lProbaFunc, Double_t *lArLo, Double_t *lArSurvivalProba, 
                         Double_t LoMax,    Int_t NbPoint, 
                         const Char_t *ch_DecayType, Double_t pMother){
    
    TString Str_DecayType(ch_DecayType);
    
    for(Int_t iPoint = 0; iPoint < NbPoint; iPoint++){
        // Lo values
        lArLo[iPoint] = LoMax/NbPoint * iPoint;
        
        // Proba
        lProbaFunc ->SetParameter(0, lArLo[iPoint] ); // Lo in meter

        if(     Str_DecayType.EqualTo("LambdaFromXi"))      lArSurvivalProba [iPoint] = lProbaFunc->Integral(0., lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("LambdaFromOmega"))   lArSurvivalProba [iPoint] = lProbaFunc->Integral(0., lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("Lambda"))            lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("K0s"))               lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("Xi"))                lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("Omega"))             lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("D0"))                lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("Dplus"))             lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("DSplus"))            lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else if(Str_DecayType.EqualTo("LambdaCplus"))       lArSurvivalProba [iPoint] = lProbaFunc->Eval( lArLo[iPoint] );
        else {
            Printf("ComputeProbability : wrong particle type... exit !");
            return 0;
        }

        
        if(! (iPoint%40) ) 
           Printf( "(%s - p(mother) = %.2f GeV/c) / [value %03d] Lo = %.6f m : P(L>Lo) = %.3f", Str_DecayType.Data(), pMother,  iPoint,  lArLo[iPoint], lArSurvivalProba [iPoint] );
    }// end loop

    return 1;
}


Int_t Root_ComputeSurvivalProbability(  TString Str_Part1ToDisplay = "kLambdaFromXi", // kLambdaFromOmega or kLambdaFromXi, kLambda, kK0s, kXi, kOmega, kD0, kDplus, kDSplus, kLambdaCplus
                                        TString Str_Part2ToDisplay = "kLambdaFromOmega",
                                        Int_t rWrite = 0
                                    ){
    
    
    if( Str_Part1ToDisplay.EqualTo("kLambdaFromXi")    || 
        Str_Part1ToDisplay.EqualTo("kLambdaFromOmega") ||
        Str_Part1ToDisplay.EqualTo("kLambda")          ||
        Str_Part1ToDisplay.EqualTo("kK0s")             ||
        Str_Part1ToDisplay.EqualTo("kXi")              ||
        Str_Part1ToDisplay.EqualTo("kOmega")           ||
        Str_Part1ToDisplay.EqualTo("kD0")              ||
        Str_Part1ToDisplay.EqualTo("kDplus")           ||
        Str_Part1ToDisplay.EqualTo("kDSplus")          ||
        Str_Part1ToDisplay.EqualTo("kLambdaCplus")
    )
        Printf("Ok, choice of 1st decaying particle, valid... one can proceed");
    else{
        Printf("Issue with 1st display choice... exit !");
        return -3;
    }
  
    if( Str_Part2ToDisplay.EqualTo("kLambdaFromXi")    || 
        Str_Part2ToDisplay.EqualTo("kLambdaFromOmega") ||
        Str_Part2ToDisplay.EqualTo("kLambda")          ||
        Str_Part2ToDisplay.EqualTo("kK0s")             ||
        Str_Part2ToDisplay.EqualTo("kXi")              ||
        Str_Part2ToDisplay.EqualTo("kOmega")           ||
        Str_Part2ToDisplay.EqualTo("kD0")              ||
        Str_Part2ToDisplay.EqualTo("kDplus")           ||
        Str_Part2ToDisplay.EqualTo("kDSplus")          ||
        Str_Part2ToDisplay.EqualTo("kLambdaCplus")     ||
        Str_Part2ToDisplay.EqualTo("")
    )
        Printf("Ok, choice of 2nd decaying particle, valid... one can proceed");
    else{
        Printf("Issue with 2nd display choice... exit !");
        return -3;
    }
  
  
  
  
    
    
    // I. - Compute the different probabilities, for the various hypotheses on pMother

    
    Int_t lNPartType = nPartDecayType;
    Int_t lNbPtCases = 3;
    Double_t LoMax  =  2.5; // in meter
  
    //-------  
    vector<vector<Double_t> > pPart;

    // NOTE : set up the size of the 2D momentum array  (lNPartType x lNbPtCases)
    pPart.resize( lNPartType );
    for (Int_t iPart = 0; iPart < lNPartType ; ++iPart) {
        pPart[iPart].resize( lNbPtCases );
    }// end set-up size of the 2D array
    
    
    pPart[kLambdaFromXi     ][0] = 1.0;     pPart[kLambdaFromXi     ][1] = 3.0;     pPart[kLambdaFromXi     ][2] =  8.0;   // NOTE Xi momentum in GeV/c, not Lambda momentum !
    pPart[kLambdaFromOmega  ][0] = 1.0;     pPart[kLambdaFromOmega  ][1] = 3.0;     pPart[kLambdaFromOmega  ][2] =  8.0;   // NOTE Omega momentum in GeV/c, not Lambda momentum !
    pPart[kLambda           ][0] = 0.5;     pPart[kLambda           ][1] = 3.0;     pPart[kLambda           ][2] =  5.0;   // primary Lambda momentum in GeV/c, diff with "Lambda in cascade" case : 0.85*pXi[0]
    pPart[kK0s              ][0] = 2.0;     pPart[kK0s              ][1] = 5.0;     pPart[kK0s              ][2] = 10.0;
    pPart[kXi               ][0] = 0.5;     pPart[kXi               ][1] = 0.8;     pPart[kXi               ][2] =  1.0;
    pPart[kOmega            ][0] = 2.0;     pPart[kOmega            ][1] = 5.0;     pPart[kOmega            ][2] = 10.0;
    pPart[kD0               ][0] = 2.0;     pPart[kD0               ][1] = 5.0;     pPart[kD0               ][2] = 10.0;
    pPart[kDplus            ][0] = 2.0;     pPart[kDplus            ][1] = 5.0;     pPart[kDplus            ][2] = 10.0;
    pPart[kDSplus           ][0] = 2.0;     pPart[kDSplus           ][1] = 5.0;     pPart[kDSplus           ][2] = 10.0;
    pPart[kLambdaCplus      ][0] = 2.0;     pPart[kLambdaCplus      ][1] = 5.0;     pPart[kLambdaCplus      ][2] = 10.0;
    
   
    
    
    //-------
    vector<vector<TGraph*> > grProba;

    // NOTE : set up the size of the 2D graph array  (lNPartType x lNbPtCases)
    grProba.resize( lNPartType );
    for (Int_t iPart = 0; iPart < lNPartType ; ++iPart) {
        grProba[iPart].resize( lNbPtCases );
    }// end set-up size of the 2D array    
          
    
    TF1 *lProbaFunc = new TF1("lProbaFunc", dProba_dlxi, 0, 1000, 3);
        lProbaFunc ->SetParNames("L0", "p(Mother)");

    Double_t lArLo[500] = {0.};
    Double_t lArSurvivalProba[500] = {0.};
    TString Str_PartType("");


    
    

    
for( Int_t iPart = 0; iPart < lNPartType ; iPart++){    
    
  
    switch(iPart){
        case kLambdaFromXi :
            LoMax = 2.5;
            Str_PartType = "LambdaFromXi";
            break;
        case kLambdaFromOmega :    
            LoMax = 2.5;
            Str_PartType = "LambdaFromOmega";
            break;
        case kLambda :
            LoMax = 2.5;
            Str_PartType = "Lambda";
            break;
        case kK0s :
            LoMax = 2.0;
            Str_PartType = "K0s";
            break;
        case kXi :
            LoMax = 2.0;
            Str_PartType = "Xi";
            break;
        case kOmega :
            LoMax = 2.0;
            Str_PartType = "Omega";
            break;            
        case kD0 :
            LoMax = 0.003;
            Str_PartType = "D0";
            break;
        case kDplus :
            LoMax = 0.006;
            Str_PartType = "Dplus";
            break;
        case kDSplus :
            LoMax = 0.003;
            Str_PartType = "DSplus";
            break;
        case kLambdaCplus :
            LoMax = 0.0012;
            Str_PartType = "LambdaCplus";
            break;
        default:
            Printf("I- Unknown particle case [%d]... exit !", iPart);
            return -11;
            
                
    } // end switch 
    
    
    
    // - 3 momenta tested 
    lProbaFunc ->SetParameter(2, iPart);

    for(Int_t ipTCase =0; ipTCase < lNbPtCases; ipTCase++){
      
      Printf(" ");
      Printf("-- Particle [%d] / Momentum case [%d] : pT = %.2f GeV/c --", iPart, ipTCase, pPart[iPart][ipTCase] ) ;
      lProbaFunc ->SetParameter(1, pPart[iPart][ipTCase]); // Particle momentum in GeV/c
      
        if( !ComputeProbability(lProbaFunc, lArLo, lArSurvivalProba, LoMax, 500, Str_PartType.Data(), lProbaFunc->GetParameter(1)   ) ){
                Printf("Sthg wrong with the ComputeProbability for %s in pT = %.2f GeV/c ", Str_PartType.Data(), lProbaFunc->GetParameter(1) ); 
                return -10;        
        }
        
        
        // NOTE Conversion des abscisses des m vers mm, pour les particules à courte distance de vol (mésons D)
        //  La conversion a lieu APRES les calculs de proba avec des cTau en mètre !
        if(iPart == kD0 || iPart == kDplus || iPart == kDSplus || iPart == kLambdaCplus)
            for(Int_t iBin = 0; iBin < 500; iBin++) lArLo[iBin] = lArLo[iBin]*1000.;
            
        grProba[iPart][ipTCase] = new TGraph(500, lArLo, lArSurvivalProba);
        grProba[iPart][ipTCase]->SetName( Form("grProba_%s_%d", Str_PartType.Data(), ipTCase) );
        Printf("-- Particle [%d] / Momentum case [%d] / graph name : %s --", iPart, ipTCase, grProba[iPart][ipTCase]->GetName() ) ;
        
    }// end for loop pTCase

    Printf("//------------------------------------------------------------------------------------------end %s",  Str_PartType.Data() );    
    Printf(" ");
 
} // end loop over particle, ~l.260













  
  // II. - Display

    myOptions(0);
    gROOT->ForceStyle();


    TCanvas *myCan = new TCanvas("myCan", "Survival Probability");
    myCan->Draw();
    myCan->cd();
    myCan->ToggleEventStatus();
    myCan->ToggleEditor();
    myCan->ToggleToolBar();


    TPad *myPad = new TPad("myPad", "The pad",0,0,1,1);
    myPadSetUp(myPad,0.15,0.04,0.04,0.15);
    
    myPad->Draw();
    myPad->cd();
    
    
    TH1F *myBlankHisto = new TH1F("myBlankHisto","Blank Histogram",100, 0, 250);
    myBlankHisto->SetNdivisions(510,"x");
    myBlankHisto->SetNdivisions(510,"y");
    myBlankHisto->Draw();
    
    TLegend *myLegend = 0x0;
  
    TLatex *system = new TLatex();
        system->SetNDC();
        //  system->SetTextFont(42);
        system->SetTextSize(0.04);
        system->SetTextColor(kGray+2);  
  
  
    Double_t xMax           =  2.5;
    Double_t xMin           = -0.1; // in meter
    Double_t yMax           = 0.75;
    Double_t lConvertFactor = 1;    // conversion from m to mm, where needed.
    TString Str_LatexPartType("");
    TString Str_LatexPartInfo("");
    TString Str_Xtitle("");
    TString Str_Ytitle("");
    Int_t lUseCaseDisplay1  = -1;
    Int_t lUseCaseDisplay2  = -1;
    Int_t lDecayType        = -1;
    Int_t lColor            = -1;  
    
  
       
    if( Str_Part1ToDisplay.EqualTo("kLambdaFromXi"     ) )  lUseCaseDisplay1 = kLambdaFromXi    ;   
    if( Str_Part1ToDisplay.EqualTo("kLambdaFromOmega"  ) )  lUseCaseDisplay1 = kLambdaFromOmega ;
    if( Str_Part1ToDisplay.EqualTo("kLambda"           ) )  lUseCaseDisplay1 = kLambda          ;
    if( Str_Part1ToDisplay.EqualTo("kK0s"              ) )  lUseCaseDisplay1 = kK0s             ;
    if( Str_Part1ToDisplay.EqualTo("kXi"               ) )  lUseCaseDisplay1 = kXi              ;
    if( Str_Part1ToDisplay.EqualTo("kOmega"            ) )  lUseCaseDisplay1 = kOmega           ;    
    if( Str_Part1ToDisplay.EqualTo("kD0"               ) )  lUseCaseDisplay1 = kD0              ;
    if( Str_Part1ToDisplay.EqualTo("kDplus"            ) )  lUseCaseDisplay1 = kDplus           ;
    if( Str_Part1ToDisplay.EqualTo("kDSplus"           ) )  lUseCaseDisplay1 = kDSplus          ;
    if( Str_Part1ToDisplay.EqualTo("kLambdaCplus"      ) )  lUseCaseDisplay1 = kLambdaCplus     ;
 
    
    if( Str_Part2ToDisplay.EqualTo("kLambdaFromXi"     ) )  lUseCaseDisplay2 = kLambdaFromXi    ;   
    if( Str_Part2ToDisplay.EqualTo("kLambdaFromOmega"  ) )  lUseCaseDisplay2 = kLambdaFromOmega ;
    if( Str_Part2ToDisplay.EqualTo("kLambda"           ) )  lUseCaseDisplay2 = kLambda          ;
    if( Str_Part2ToDisplay.EqualTo("kK0s"              ) )  lUseCaseDisplay2 = kK0s             ;
    if( Str_Part2ToDisplay.EqualTo("kXi"               ) )  lUseCaseDisplay2 = kXi              ;
    if( Str_Part2ToDisplay.EqualTo("kOmega"            ) )  lUseCaseDisplay2 = kOmega           ;        
    if( Str_Part2ToDisplay.EqualTo("kD0"               ) )  lUseCaseDisplay2 = kD0              ;
    if( Str_Part2ToDisplay.EqualTo("kDplus"            ) )  lUseCaseDisplay2 = kDplus           ;
    if( Str_Part2ToDisplay.EqualTo("kDSplus"           ) )  lUseCaseDisplay2 = kDSplus          ;
    if( Str_Part2ToDisplay.EqualTo("kLambdaCplus"      ) )  lUseCaseDisplay2 = kLambdaCplus     ;
    if( Str_Part2ToDisplay.EqualTo(""                  ) )  lUseCaseDisplay2 = -1;
    
    


for(Int_t iDraw = 0; iDraw < 2; iDraw++){  
    
    Int_t lUseCaseDisplay = -1;
    if(iDraw == 0) lUseCaseDisplay = lUseCaseDisplay1;
    if(iDraw == 1) {
        lUseCaseDisplay = lUseCaseDisplay2;
        if(lUseCaseDisplay < 0) continue;   
    }
    
    switch(lUseCaseDisplay){
        case kLambdaFromXi :
            lDecayType = kLambdaFromXi;
            lColor = colors[ red ];
            xMax =  2.5; // en mètre
            xMin = -0.1;
            yMax = 0.8;
            if(iDraw == 0)  lConvertFactor = 1;
            Str_PartType = "LambdaFromXi";
            Str_LatexPartType = "#Xi^{-}";
            Str_LatexPartInfo = Form("#Lambda from #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (m) #scale[0.7]{at y = 0} ";
            Str_Ytitle = Form(" #color[%d]{P_{#Lambda} [ L >L_{0} , #font[52]{p}_{T}(#Xi) ]}", lColor);
            break;
        case kLambdaFromOmega : 
            lDecayType = kLambdaFromOmega;
            lColor =  colors[ azure ];
            xMax =  2.5; // en mètre
            xMin = -0.1;
            yMax = 0.8;
            if(iDraw == 0)  lConvertFactor = 1;
            Str_PartType = "LambdaFromOmega";
            Str_LatexPartType = "#Omega^{-}";
            Str_LatexPartInfo = Form("#Lambda from #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (m) #scale[0.7]{at y = 0} ";
            Str_Ytitle = Form("#color[%d]{P_{#Lambda} [ L >L_{0} , #font[52]{p}_{T}(#Omega) ]}", lColor);
            break;
        case kLambda :
            lDecayType = kLambda;
            lColor = colors[ orange ]; 
            xMax =  2.5; // en mètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1;
            Str_PartType = "Lambda";
            Str_LatexPartType = "#Lambda";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (m) #scale[0.7]{at y = 0} ";
            Str_Ytitle = Form("#color[%d]{P_{#Lambda} [ L >L_{0} , #font[52]{p}_{T}(#Lambda) ]}", lColor);       
            break;
        case kK0s :
            lDecayType = kK0s;
            lColor = colors[ yellow ];
            xMax =  2.0; // en mètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1;
            Str_PartType = "K0s";
            Str_LatexPartType = "K^{0}_{S}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (m) #scale[0.7]{at y = 0} ";
            Str_Ytitle = Form("#color[%d]{P_{K^{0}_{S}} [ L >L_{0} , #font[52]{p}_{T}(K^{0}_{S}) ]}", lColor);
            break;
        case kXi :
            lDecayType = kXi;
            lColor = colors[ blue ]; 
            xMax =  2.0; // en mètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1;
            Str_PartType = "Xi";
            Str_LatexPartType = "#Xi^{-}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (m) #scale[0.7]{at y = 0} ";
            Str_Ytitle = Form("#color[%d]{P_{#Xi} [ L >L_{0} , #font[52]{p}_{T}(#Xi) ]}", lColor);       
            break;  
        case kOmega :
            lDecayType = kOmega ;
            lColor = colors[ violet ]; 
            xMax =  2.0; // en mètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1;
            Str_PartType = "Omega";
            Str_LatexPartType = "#Omega^{-}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (m) #scale[0.7]{at y = 0} ";
            Str_Ytitle = Form("#color[%d]{P_{#Omega} [ L >L_{0} , #font[52]{p}_{T}(#Omega) ]}", lColor);       
            break;                         
        case kD0 :
            lDecayType = kD0;
            lColor = colors[ blue ];
            xMax =  3;  // en millimètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0) lConvertFactor = 1000;
            Str_PartType = "D0";
            Str_LatexPartType = "D^{0}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (mm)";
            Str_Ytitle = Form( "#color[%d]{P_{D^{0}} [ L >L_{0} , #font[52]{p}_{T}(D^{0}) ]}", lColor);
            break;
        case kDplus :
            lDecayType = kDplus;
            lColor = colors[ violet ];
            xMax =  6; // en millimètre
            xMin = -0.2;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1000;
            Str_PartType = "Dplus";
            Str_LatexPartType = "D^{+}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (mm)";
            Str_Ytitle = Form("#color[%d]{P_{D^{+}} [ L >L_{0} , #font[52]{p}_{T}(D^{+}) ]}", lColor);
            break;
        case kDSplus :
            lDecayType = kDSplus;
            lColor = colors[ orange ] ;
            xMax =  3; // en millimètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1000;
            Str_PartType = "DSplus";
            Str_LatexPartType = "D^{+}_{S}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (mm)";
            Str_Ytitle = Form("#color[%d]{P_{D^{+}_{S}} [ L >L_{0} , #font[52]{p}_{T}(D^{+}_{S}) ]}", lColor);
            break;
        case kLambdaCplus :
            lDecayType = kLambdaCplus;
            lColor = colors[ green ];
            xMax =  1.2; // en millimètre
            xMin = -0.1;
            yMax = 1.0;
            if(iDraw == 0)  lConvertFactor = 1000;
            Str_PartType = "LambdaCplus";
            Str_LatexPartType = "#Lambda^{+}_{C}";
            Str_LatexPartInfo = Form("primary #color[%d]{%s}", lColor, Str_LatexPartType.Data() );
            if(iDraw == 0)  Str_Xtitle = "L_{0} (mm)";
            Str_Ytitle = Form("#color[%d]{P_{#Lambda^{+}_{C}} [ L >L_{0} , #font[52]{p}_{T}(#Lambda^{+}_{C}) ]}", lColor);
            break;
        default:
            Printf("II - Unknown particle case [%d]... exit !", lUseCaseDisplay1);
            return -11;
            
                
    } // end switch 
    
  
  if(iDraw == 0){
    myBlankHisto->SetBins(100, xMin, xMax);
    myBlankHisto->SetMaximum( yMax );
    myBlankHisto->SetXTitle( Str_Xtitle.Data() );
    myBlankHisto->SetYTitle( Str_Ytitle.Data() );
    myBlankHisto->Draw();
  }
  
  if(iDraw == 1) {
      Str_Ytitle.Append( Form(" or %s", myBlankHisto->GetYaxis()->GetTitle()  ) );
      myBlankHisto->SetYTitle( Str_Ytitle.Data() );
  }
  
  
   
  
  
    //myGraphSetUp(grProba[kLambdaFromXi][0], Float_t currentMarkerSize, 
    //                            int currentMarkerStyle, int currentMarkerColor, 
    //                            int currentLineStyle,   int currentLineColor, int currentLineWidth);
  
        myGraphSetUp(grProba[lDecayType][0], 0.2, kOpenSquare, lColor, 7,  lColor, 2);
        myGraphSetUp(grProba[lDecayType][1], 0.2, kOpenSquare, lColor, 2,  lColor, 2);
        myGraphSetUp(grProba[lDecayType][2], 0.2, kOpenSquare, lColor, 1,  lColor, 2);
        
        grProba[lDecayType][0]->Draw("C");
        grProba[lDecayType][1]->Draw("C");
        grProba[lDecayType][2]->Draw("C");
         
        if(iDraw == 0){
            if(lUseCaseDisplay2 < 0)    myLegend = new TLegend(0.62,0.6, 0.95,0.85);
            else                        myLegend = new TLegend(0.62,0.49, 0.95,0.92);
        }
        
        if(iDraw == 0) {
            system->DrawLatex(0.35, 0.89, "Survival probability of " );
            system->DrawLatex(0.39, 0.83, Form("%s", Str_LatexPartInfo.Data() ) );
        }
        if(iDraw == 1) system->DrawLatex(0.39, 0.78, Form("and %s", Str_LatexPartInfo.Data() ) );
        
        if(iDraw == 1) myLegend->AddEntry(myBlankHisto, " ", "");
        myLegend->AddEntry(grProba[lDecayType][0], Form("#font[52]{p}_{T}(%s) = %.1f GeV/#it{c}", Str_LatexPartType.Data(), pPart[lDecayType][0] ), "l");
        myLegend->AddEntry(grProba[lDecayType][1], Form("#font[52]{p}_{T}(%s) = %.1f GeV/#it{c}", Str_LatexPartType.Data(), pPart[lDecayType][1] ), "l");
        myLegend->AddEntry(grProba[lDecayType][2], Form("#font[52]{p}_{T}(%s) = %.1f GeV/#it{c}", Str_LatexPartType.Data(), pPart[lDecayType][2] ), "l");
         

  
   
  


}// end loop iDraw

    myLegendSetUp(myLegend,0.04);
    myLegend->Draw();
    
    
    TLatex *lComment2 = new TLatex();
        //lComment2->SetNDC();     
        lComment2->SetTextColor(kGray+3);
        lComment2->SetTextSize(0.027);
        lComment2->SetTextFont(52);
        lComment2->SetTextAngle(90);
    
    TLine *lineDetector = new TLine(1.00*lConvertFactor, 0.02, 1.00*lConvertFactor, 0.60);
        lineDetector->SetLineStyle( 1 );
        lineDetector->SetLineColor( kGray+1 );
        lineDetector->SetLineWidth( 1 );
        // lineDetector->Draw("same");
        
        
        
       

    lineDetector ->DrawLine(0.029*lConvertFactor, 0.02, 0.029*lConvertFactor, 0.60);
    lComment2->DrawLatex(0.01*lConvertFactor,0.52,"Beam pipe");
    
    lineDetector ->DrawLine(0.039*lConvertFactor, 0.02, 0.039*lConvertFactor, 0.65);
    lineDetector ->DrawLine(0.076*lConvertFactor, 0.02, 0.076*lConvertFactor, 0.65);
    lComment2->DrawLatex(0.09*lConvertFactor,0.66,"SPD_{1 + 2}");
    
    lineDetector ->DrawLine(0.150*lConvertFactor, 0.02, 0.150*lConvertFactor, 0.62);
    lineDetector ->DrawLine(0.239*lConvertFactor, 0.02, 0.239*lConvertFactor, 0.62);
    lComment2->DrawLatex(0.22*lConvertFactor,0.63,"SDD_{1 + 2}");
    
    lineDetector ->DrawLine(0.38*lConvertFactor, 0.02, 0.38*lConvertFactor, 0.63);
    lineDetector ->DrawLine(0.43*lConvertFactor, 0.02, 0.43*lConvertFactor, 0.63);
    lComment2->DrawLatex(0.43*lConvertFactor,0.64,"SSD_{1 + 2}");
    
    lineDetector ->DrawLine(0.848*lConvertFactor, 0.02, 0.848*lConvertFactor, 0.41);           
    lComment2->DrawLatex(0.91*lConvertFactor,0.35,"TPC");




  if (rWrite == 1)  myCan->SaveAs( Form("SurvivalProba-%s-%s.eps", Str_Part1ToDisplay.Data(), Str_Part2ToDisplay.Data() ) );
  if (rWrite == 2)  myCan->SaveAs( Form("SurvivalProba-%s-%s.png", Str_Part1ToDisplay.Data(), Str_Part2ToDisplay.Data() ) );
 
    return 1;

  
}




void myLegendSetUp(TLegend *currentLegend, float currentTextSize){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  //currentLegend->SetFillStyle(0);
  currentLegend->SetFillStyle(1001);
  currentLegend->SetFillColor(0);
  currentLegend->SetLineWidth(1);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  //currentLegend->SetNColumns(2);  
  
  return;
}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, 
                float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  currentPad->SetFillColor(0);
  currentPad->SetBorderMode(0);
  currentPad->SetBorderSize(2);
   //currentPad->SetLogx();
  currentPad->SetTicks();
  currentPad->SetFrameBorderMode(0);
  currentPad->SetFrameBorderMode(0);    
  currentPad->SetGridx();
  currentPad->SetGridy();
  
  return;
}

void myGraphSetUp(TGraph *currentGraph, Float_t currentMarkerSize, 
                  int currentMarkerStyle, int currentMarkerColor, 
                  int currentLineStyle,   int currentLineColor, int currentLineWidth){
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  currentGraph->SetLineWidth(currentLineWidth);
  return;
}

void myOptions(Int_t lStat){
  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.1,"xyz");  
  gStyle->SetTitleSize(0.05,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}
