// emacs -*- C++ -*-

/**
   Script to save the post-fit histograms and yields.

   This are the main steps in the script:
   0. read in workspace
   1. perform a standard fit and save a snapshot
   2. for category in [SR, CR, ...]:
        for component in [background components]:
           compute integral and error varying only these component
           collect components that should be reported together (eg. multiple signal production processes)
        3. compute integral and error for merged components
        4. perform two global fits: one with muhat, and one with mu=0


   How to run:
   root -l -b
   root [0] .L DumpPostFitHistos.C+
   root [1] DumpPostFitHistos("workspaces/taumu_lh_ll/ws_LFV_combined_AllSYS_model.root", "", "taumu_lh_ll.root", "", true, true, false, false)

   Alternatively, as standalone:
   compile with
   g++ -Wall `root-config  --cflags` -o DumpPostFitHistos DumpPostFitHistos.C `root-config --libs` -lRooFitCore -lRooFit -lRooStats -lHistFactory
   run with
   ./DumpPostFitHistos -h

   original version:
   ruthmann@cern.ch for the htautau group, 2013
   hlfv version:
   davide.gerbaudo@gmail.com,  Dec 2015
*/

// C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <map>
#include <sstream>

// Root
#include "TFile.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TMarker.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"

// RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "Roo1DTable.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooRandom.h"
#include "TStopwatch.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooBinning.h"
#include "RooAddition.h"
// RooStat
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
//#include "RooStats/MinNLLTestStat.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;


TH1D *TGraphAsymmErrorsToTH1D(TGraphAsymmErrors *g, TH1D *hmodel);
TH1D *TGraphToTH1D(TGraph *g, TH1D *hmodel);
void  ExtractHistoFromCanvas(TPad *can, vector<TH1*> &vec_histo, TGraphAsymmErrors* &gdataAsymErr);

RooFitResult* FitPDF( ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata,
                      int &MinuitStatus, int &HessStatus, double &Edm,
                      TString minimType = "Minuit2", bool useMinos = false );
void TransportNPVals(RooArgSet* inSet, RooArgSet* outSet);
void DumpHistograms(RooWorkspace* w,ModelConfig*mc,  RooSimultaneous* pdfVar, RooAbsData* data, TString outFile, RooFitResult* fitresGlobal=0, TString tag="", bool convertaxis=false);
TString workspaceName = "combined";
TString modelConfigName = "ModelConfig";
TString ObsDataName = "obsData";//"asimov125_1";
RooArgSet* params=0;

void DumpPostFitHistos(TString workspaceIn, TString workspaceVar, TString outFile, TString snapshotName="", bool doFit=true,bool transportCovariance=false, bool JustFit=false, bool convertaxis=false, bool fixNuisanceParameters=false  );
TH1* ConvertBDTAxis(TH1 *hist);

//======================================================
// ================= Main function =====================
//======================================================

void DumpPostFitHistos(TString workspaceIn, TString workspaceVar, TString outFile, TString snapshotName, bool doFit, bool transportCovariance, bool JustFit, bool convertaxis, bool fixNuisanceParameters ){
    /*
      workspaceIn: Workspace to load
      workspaceVar: This is the workspace where the NPs should be transferred to
      if "" the 1st workspace will be dumped
      outFile: Name the outputfile
      snapshotName:  snpshot which will be tried to load
      doFit: in all cases do a fit
    */
    cout<<"DumpPostFitHistos() starting at "<<TTimeStamp().AsString()<<endl
        <<" running options:"<<endl
        <<" workspaceIn '"<<workspaceIn<<"'"
        <<" workspaceVar '"<<workspaceVar<<"'"
        <<" outFile '"<<outFile<<"'"
        <<" snapshotName '"<<snapshotName<<"'"
        <<" doFit '"<<doFit<<"'"
        <<" transportCovariance '"<<transportCovariance<<"'"
        <<" JustFit '"<<JustFit<<"'"
        <<" convertaxis '"<<convertaxis<<"'"
        <<endl;
    // Load workspaceIn
    TFile*file_wsIn = TFile::Open(workspaceIn);
    if (!file_wsIn) {
        cout << "The file " << workspaceIn << " is not found/created, will stop here." << endl;
        return;
    }
    if(!(RooWorkspace*) file_wsIn->Get(workspaceName)){
        cout <<"workspace not found" << endl;
        return;
    }
    RooWorkspace* wIn      = (RooWorkspace*) file_wsIn->Get(workspaceName);
    ModelConfig* mcIn     = (ModelConfig*) wIn->obj(modelConfigName);
    RooAbsData* dataIn   = wIn->data(ObsDataName);
    RooSimultaneous* pdfIn = (RooSimultaneous*) wIn->pdf("simPdf");
    RooArgSet* paramsIn = (RooArgSet*) pdfIn->getParameters(*dataIn) ;

    RooWorkspace* wVar =0;
    ModelConfig* mcVar   =0;
    RooAbsData* dataVar =0;
    RooSimultaneous* pdfVar=0;
    RooArgSet* paramsVar=0;
    if(workspaceVar!=""){
        // Load workspaceVar
        TFile*file_wsVar = TFile::Open(workspaceVar);
        if (!file_wsVar) {
            cout << "The file " << workspaceVar << " is not found/created, will stop here." << endl;
            return;
        }
        if(!(RooWorkspace*) file_wsVar->Get(workspaceName)){
            cout <<"workspace not found" << endl;
            return;
        }
        wVar      = (RooWorkspace*) file_wsVar->Get(workspaceName);
        mcVar     = (ModelConfig*) wVar->obj(modelConfigName);
        dataVar   = wVar->data(ObsDataName);
        pdfVar = (RooSimultaneous*) wVar->pdf("simPdf");
        paramsVar = (RooArgSet*) pdfVar->getParameters(*dataIn) ;
        wVar->saveSnapshot("snapshot_paramsVals_initial",*paramsVar);
    }


    RooFitResult* r=0;
    // Try to load snapshot
    // Need to do extra work because Swagato delivered this in the a non-initial state
    // if(!wIn->loadSnapshot("nominalGlobs")){cout<<"Swagatos advertised snapshot nominalNuis is not available"<<endl; };
    // if(!wIn->loadSnapshot("nominalNuis")){cout<<"Swagatos advertised snapshot nominalNuis is not available"<<endl; };
    RooRealVar * poi = (RooRealVar*) mcIn->GetParametersOfInterest()->first();
    //  poi->setVal(1);
    poi->setRange(-100,100);


    TIterator* InIter = paramsIn->createIterator();
    RooRealVar* var;
    while ((var = (RooRealVar*)InIter->Next())){
        TString name =  var->GetName();
        if( name.Contains("SigXsecOverSM") ){
            if( name.Contains("WW") ){
                var->setConstant(1);
                cout<<"Setting WW mu to 1 "<< name<<endl;

            }
            var->setRange(-200,200);
            cout<<"Adjusting range of "<< name<<endl;
        }

        // if( name.Contains("mu_XS8_VBF") || name.Contains("mu_XS8_ggH") ){
        //   var->setRange(-200,200);
        //   var->setConstant(0);

        // }
    }
    //Set mu constant to 1:


    //  poi->setConstant();


    InIter = paramsIn->createIterator();
    while ((var = (RooRealVar*)InIter->Next())){
        TString name =  var->GetName();
        if( name.Contains("ATLAS_norm_") ){
            cout<<" Recover initial state: Setting "<<name<<" to 1."<<endl;
            var->setVal(1.);
        }
    }

    cout<<"Initial parameters?"<<endl;
    paramsIn->Print("v");
    if(!wIn->loadSnapshot("snapshot_paramsVals_initial"))
        wIn->saveSnapshot("snapshot_paramsVals_initial",*paramsIn);

    if (!doFit && !transportCovariance){
        if(!  wIn->loadSnapshot(snapshotName) ){
            cout<<" ERROR: Was unable to load snapshot "<<snapshotName<<" from wsIn"<<endl;
        }
    }
    else{ // Do a fit and save the snapshot afterwards
        int MinuitStat,HessStat;
        double Edm;
        r =  FitPDF(  mcIn,  pdfIn, dataIn,
                      MinuitStat, HessStat, Edm,
                      "Minuit2", false);
        wIn->saveSnapshot(snapshotName, *paramsIn);
        if( MinuitStat !=0 && MinuitStat !=1) cout<<" ERROR: Fit Failed"<<endl;

        cout<<"Fit result:"<<endl;
        r->Print();
        cout<<"---"<<endl;

        TString tmpFileName= outFile;
        tmpFileName.ReplaceAll(".root","_fitresult.txt");
        ofstream myfile;
        myfile.open( tmpFileName );
        myfile<<"muval="<<poi->getVal()<<" + "<<poi->getErrorHi()<<" - "<<poi->getErrorLo()<<endl;
        myfile.close();
        if(JustFit) return;
    }

    struct IsNuisanceParameter{
        bool operator()(const TString &n) {
            return (n.BeginsWith("alpha_ATLAS") or n.BeginsWith("alpha_Fakes") or
                    n.BeginsWith("alpha_QCDscale") or n.BeginsWith("alpha_pdf") or
                    n.BeginsWith("fl1pt_l1pt"));
        }
    };
    struct {
        void operator()(TIterator* it) {
            while (RooRealVar *v = static_cast<RooRealVar*>(it->Next())){
                if(IsNuisanceParameter()(v->GetName())) {
                    v->setConstant(true);
                    cout<<"Fixing nuisance parameter '"<<v->GetName()<<"' to constant"<<endl;
                }
            }
        }
    } fixSystNuisanceParameters;
    if(fixNuisanceParameters) {
        cout<<"---------------------------------------"<<endl
            <<"Fixing the systematic nuisance parameters and then running another fit"<<endl
            <<"---------------------------------------"<<endl;
        fixSystNuisanceParameters(paramsIn->createIterator());
        int MinuitStat,HessStat;
        double Edm;
        r =  FitPDF(  mcIn,  pdfIn, dataIn,
                      MinuitStat, HessStat, Edm,
                      "Minuit2", false);
        wIn->saveSnapshot(snapshotName, *paramsIn);
        if( MinuitStat !=0 && MinuitStat !=1) cout<<" ERROR: Fit Failed"<<endl;

        cout<<"---------------------------------------"<<endl
            <<"Fit result:"<<endl;
        r->Print();
        cout<<"---------------------------------------"<<endl;
        DumpHistograms(wIn, mcIn, pdfIn, dataIn, outFile, r, "_noSyst" , convertaxis);
    }

    if(workspaceVar==""){
        cout<<"---------------------------------------"<<endl
            <<"Done with the preliminary fits; now going to dump the histograms"<<endl
            <<"  first DumpHistograms for the unconstrained fit..."<<endl
            <<"---------------------------------------"<<endl;
        wIn->saveSnapshot("THESnapShot",*paramsIn);
        params=paramsIn;
        DumpHistograms(wIn, mcIn, pdfIn, dataIn, outFile, r ,"", convertaxis);
        // Crazy idea: do another fit:
        cout<<"---------------------------------------"<<endl
            <<"  then DumpHistograms for the mu=0 fit..."<<endl
            <<"---------------------------------------"<<endl;
        RooRealVar * poi = (RooRealVar*) mcIn->GetParametersOfInterest()->first();
        poi->setVal(0);
        poi->setConstant(true);
        if(fixNuisanceParameters) fixSystNuisanceParameters(paramsIn->createIterator());
        int MinuitStat,HessStat;
        double Edm;
        r =  FitPDF(  mcIn,  pdfIn, dataIn,
                      MinuitStat, HessStat, Edm,
                      "Minuit2", false);
        if (!  wIn->saveSnapshot("THESnapShot",*paramsIn)){cout<<"Failed to overwrite the snapshot"<<endl;}
        DumpHistograms(wIn, mcIn, pdfIn, dataIn, outFile, r, "_condMu0" , convertaxis);
        // one more fit, now with mu=1 (to get the exp signal)
        wIn->loadSnapshot("THESnapShot"); // needed?
        poi->setConstant(false);
        poi->setVal(1);
        poi->setConstant(true);
        if(fixNuisanceParameters) fixSystNuisanceParameters(paramsIn->createIterator());
        r =  FitPDF(  mcIn,  pdfIn, dataIn,
                      MinuitStat, HessStat, Edm,
                      "Minuit2", false);
        if (!  wIn->saveSnapshot("THESnapShot",*paramsIn)){cout<<"Failed to overwrite the snapshot"<<endl;}
        DumpHistograms(wIn, mcIn, pdfIn, dataIn, outFile, r, "_condMu1" , convertaxis);


        return;
    }

    //  paramsIn->Print("v");
    // Now we have a set of parameters in the state we want to transport to the second set
    TransportNPVals(paramsIn, paramsVar);
    cout<<" Second WS Parameters: "<<endl;
    paramsVar->Print("v");
    wVar->saveSnapshot("THESnapShot",*paramsIn);
    if(!transportCovariance) r=0;
    params=paramsVar;
    DumpHistograms(wVar, mcVar, pdfVar, dataVar, outFile, r, "", convertaxis );

}

void DumpHistograms(RooWorkspace* w,ModelConfig*mc,  RooSimultaneous* pdfVar, RooAbsData* data, TString outFile, RooFitResult* fitresGlobal, TString tag , bool convertaxis){
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<" DumpHistograms........."<<endl;
    string linebreak="--------------------------------------------------\n";
    RooCategory* channelCat = (RooCategory*) (&pdfVar->indexCat());
    TIterator *iter = channelCat->typeIterator() ;
    RooCatType *tt  = NULL;

    const RooArgSet* nuis=  mc->GetNuisanceParameters();
    int channelCounter = -1;
    while((tt=(RooCatType*) iter->Next()) ){
        channelCounter++;
        TString tmpFileName= outFile;
        tmpFileName.ReplaceAll(".root",TString(tt->GetName())+".root");
        TFile* file=0;
        if (tag!=""){
            file=new TFile(tmpFileName,"UPDATE");
        }
        else
            file=new TFile(tmpFileName,"RECREATE");
        TString modelName(tt->GetName());
        modelName.Append("_model");

        cout<<"\n\n"<<linebreak+linebreak<<"   On Category["<<channelCounter<<"]: "<<modelName<<"\n\n"<<linebreak+linebreak<<endl;
        if( modelName.Contains("btop") || modelName.Contains("vtop") || modelName.Contains("vzll") || modelName.Contains("bzll") ){
            cout<<" ... skip"<<endl;
            continue;
        }

        RooAbsPdf  *pdftmp  = pdfVar->getPdf( tt->GetName() );
        if (!pdftmp) cout<<"ERROR: was not able to getPDF"<<endl;
        RooAbsData *datatmp = data->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),tt->GetName()));
        if (!datatmp) cout<<"ERROR: was not able to get data"<<endl;
        RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
        if (!obstmp) cout<<"ERROR: was not able to get observables"<<endl;
        RooRealVar *obs     = ((RooRealVar*) obstmp->first());
        if (!obs) cout<<"ERROR: was not able to get observable"<<endl;
        RooRealVar* binWidth = ((RooRealVar*) pdftmp->getVariables()->find(Form("binWidth_obs_x_%s_0",tt->GetName()))) ;
        if(!binWidth) { cout << "No bin width " << tt->GetName() << endl; }


        RooRealSumPdf *pdfmodel = (RooRealSumPdf*) (pdftmp->getComponents())->find(modelName);
        if(!pdfmodel){
            cout<<"ERROR: could not get pdf"<<endl;
            pdftmp->getComponents()->Print("v");
        }
        RooRealVar* binWidth_model = ((RooRealVar*) pdfmodel->getVariables()->find(Form("binWidth_obs_x_%s_0",tt->GetName()))) ;
        if(!binWidth_model) { cout << "No model binwidth " << tt->GetName() << endl; }

        RooArgList funcList =  pdfmodel->funcList();
        // RooArgSet* components = pdfmodel->getComponents();
        //components->Print("v");
        // cout<<"------"<<endl;
        RooProduct* comp = 0;
        RooLinkedListIter funcIter = funcList.iterator() ;

        // cout<<"Starting wouters trick"<<endl;
        // It would be great to somehow get the pdfs directly
        // // Given RooAbsPdf* pdf and RooRealVar* obs
        // cout<<"ObsMin="<<obs->getMin()<<" ObsMax="<<obs->getMax()<<endl;
        // list<Double_t>* bl = pdfmodel->binBoundaries(*obs,obs->getMin(),obs->getMax()) ;
        // Double_t* ba = new Double_t[bl->size()] ; int i=0 ;
        // for (list<double>::iterator it=bl->begin() ; it!=bl->end() ; ++it) { ba[i++] = *it ; }
        // RooBinning binning(bl->size()-1,ba) ;
        // cout<<"End wouter trick"<<endl;

        TCanvas* c2 = new TCanvas( "dummy", "dummy", 500,500 );
        RooPlot* frame = obs->frame();

        float postFitIntegral = pdftmp->expectedEvents(*obs);

        datatmp->plotOn(frame,MarkerSize(1),Name("Data"),DataError(RooAbsData::Poisson));

        double chi2 = frame->chiSquare();



        // cout << "     Post Fit Yields" << endl;
        w->loadSnapshot("THESnapShot");
        ofstream myfile;
        myfile.open(tmpFileName.ReplaceAll(".root",tag+".txt") );


        funcIter = funcList.iterator() ;

        //Hacky and name dependent way to get correct uncertainties on sums of components
        RooArgList signalComps;
        RooArgList bkdComps;
        RooArgList zttComps;

        float nsig=0;

        int categoryCounter = -1;
        while( (comp = (RooProduct*) funcIter.Next()) ) {
            categoryCounter++;
            TString compname(comp->GetName());
            cout<<"\n"<<linebreak<<"On category["<<categoryCounter<<"]: "<<compname<<"\n"<<linebreak<<endl;

            compname.ReplaceAll("L_x_","");
            compname.ReplaceAll(tt->GetName(),"");
            compname.ReplaceAll("_overallSyst_x_StatUncert","");
            compname.ReplaceAll("_overallSyst_x_HistSyst","");
            compname.ReplaceAll("_overallSyst_x_Exp","");

            // This here is ugly and hacky.. depends on the naming. I twill work though for Swagatos combination in November 2013
            struct {
                bool operator()(const TString &n) {
                    return ((n.Contains("WH") || n.Contains("ZH") || n.Contains("ggH") || n.Contains("VBF"))
                            and not n.Contains("HWW"));
                }
            } isSignalHtautau;
            struct {
                bool operator()(const TString &n) {
                    return (n.Contains("signal") and (n.Contains("ggH") || n.Contains("vbf") || n.Contains("VH")));
                }
            } isSignalHlfv;

            if(isSignalHlfv(compname)) {
                signalComps.add(*comp);
                nsig+= (comp->createIntegral(*obs))->getVal() * binWidth->getVal() ;
            } else {
                bkdComps.add(*comp);
            }
            if ( compname.Contains("Ztt") ) { zttComps.add(*comp); }

            RooAbsReal* integral = comp->createIntegral(*obs);
            double nom= integral->getVal() * binWidth->getVal();
            double 	err_correct=0;
            //      cout << "\t" << comp->GetName() << "\t" << (comp->createIntegral(*obs))->getVal() * binWidth->getVal() << endl;

            // Get the correct event yield:
            double Ntemp=(comp->createIntegral(*obs))->getVal() * binWidth->getVal();
            pdfmodel->plotOn(frame,LineWidth(0),Components(*comp),LineColor(0), LineStyle(0), Normalization(Ntemp,RooAbsReal::NumEvent),Name("NoStacked_"+compname));
            if(fitresGlobal){
                // This here does not work: gives the wrong uncertainty because the normalisation to NTemp is forced - even for varied NPs..
                //pdfmodel->plotOn(frame,LineWidth(0),Components(*comp),LineColor(0), LineStyle(0), Normalization(Ntemp,RooAbsReal::NumEvent),VisualizeError(*fitresGlobal,1),Name("FitError_AfterFit_"+compname));

                err_correct = integral->getPropagatedError(*fitresGlobal )* binWidth->getVal();

                // That was Nicolas idea and works:
                comp->plotOn(frame,LineWidth(0),LineColor(0), LineStyle(0), Normalization(1,RooAbsReal::RelativeExpected ),VisualizeError(*fitresGlobal,1),Name("FitError_AfterFit_"+compname));


            }
            myfile<<compname<<" = "<<nom <<" +- "<<err_correct<<endl;
            TH1F tmph=TH1F(TString("Integral_")+compname+tag,TString("Integral_")+compname+tag,1,0,1);
            tmph.SetBinContent(1,nom);
            tmph.SetBinError(1,err_correct);
            tmph.Write();
        } // Loop through components


        RooRealVar * poi = (RooRealVar*) mc->GetParametersOfInterest()->first();

        // Get total expected signal and its error:
        // This here is ugly and hacky.. depends on the naming. I twill work though for Swagatos combination in November 2013
        if(fitresGlobal){


            // Try Nicolas hack:
            if( signalComps.getSize()>1 ){
                RooAddition* tmp_pdfSum = new RooAddition( "pdf_Sum", "pdf_Sum", signalComps);
                RooAbsReal* integral = tmp_pdfSum->createIntegral(*obs);
                double nom= integral->getVal() * binWidth->getVal();
                double err_correct = integral->getPropagatedError(*fitresGlobal )* binWidth->getVal();;
                myfile<<"Signal = "<<nom <<" +- "<<err_correct<<endl;
                tmp_pdfSum->plotOn(frame, VisualizeError(*fitresGlobal,1),  Normalization(1,RooAbsReal::RelativeExpected),Name("FitError_AfterFit_FullSignal") );

                TH1F tmph=TH1F(TString("Integral_TotSignal")+tag,TString("Integral_TotSignal")+tag, 1,0,1);
                tmph.SetBinContent(1,nom);
                tmph.SetBinError(1,err_correct);
                tmph.Write();
            }
            if( bkdComps.getSize()>1 ){
                RooAddition* tmp_pdfSum = new RooAddition( "pdf_Sum_bkd", "pdf_Sum_bkd", bkdComps);
                RooAbsReal* integral = tmp_pdfSum->createIntegral(*obs);
                double nom= integral->getVal() * binWidth->getVal();
                double err_correct = integral->getPropagatedError(*fitresGlobal )* binWidth->getVal();;
                myfile<<"Bkd = "<<nom <<" +- "<<err_correct<<endl;

                TH1F tmph=TH1F(TString("Integral_TotBkd")+tag,TString("Integral_TotBkd")+tag, 1,0,1);
                tmph.SetBinContent(1,nom);
                tmph.SetBinError(1,err_correct);
                tmph.Write();
            }

            if( zttComps.getSize()>1 ){
                RooAddition* tmp_pdfSum_ztt = new RooAddition( "pdf_Sum_ztt", "pdf_Sum_ztt", zttComps);
                tmp_pdfSum_ztt->plotOn(frame, VisualizeError(*fitresGlobal,1),  Normalization(1,RooAbsReal::RelativeExpected),Name("FitError_AfterFit_Ztt") );

                RooAbsReal* integral = tmp_pdfSum_ztt->createIntegral(*obs);
                double nom= integral->getVal() * binWidth->getVal();
                double err_correct = integral->getPropagatedError(*fitresGlobal )* binWidth->getVal();;
                myfile<<"Ztt = "<<nom <<" +- "<<err_correct<<endl;
                TH1F tmph=TH1F(TString("Integral_TotZtt")+tag,TString("Integral_TotZtt")+tag, 1,0,1);
                tmph.SetBinContent(1,nom);
                tmph.SetBinError(1,err_correct);
                tmph.Write();
            }
        }
        // //The combined model
        if(fitresGlobal){
            pdftmp->plotOn(frame,FillColor(kOrange),LineWidth(2),LineColor(kBlue),VisualizeError(*fitresGlobal,1),
                           Normalization(1,RooAbsReal::RelativeExpected) ,Name("FitError_AfterFit"));
            RooAbsReal*  integral= pdfmodel->createIntegral(*obs);
            //      double  nom= integral->getVal() * binWidth_model->getVal();
            //      double  err_correct = integral->getPropagatedError(*fitresGlobal )* binWidth->getVal();
            double  nom= integral->getVal() ;
            double  err_correct = integral->getPropagatedError(*fitresGlobal );
            myfile<<"Total Model muhat (with stat NP) = "<<nom <<" +- "<<err_correct<<endl;
            TH1F tmph=TH1F(TString("Integral_TotModelMuHat")+tag,TString("Integral_TotModelMuHat")+tag, 1,0,1);
            tmph.SetBinContent(1,nom);
            tmph.SetBinError(1,err_correct);
            tmph.Write();

            // Now also get the error for bkd only
            if (poi)
                poi->setVal(0.);
            integral= pdfmodel->createIntegral(*obs);
            // nom= integral->getVal() * binWidth->getVal();
            // err_correct = integral->getPropagatedError(*fitresGlobal )* binWidth->getVal();
            nom= integral->getVal();
            err_correct = integral->getPropagatedError(*fitresGlobal);
            myfile<<"Total Model mu0 (with stat NP) = "<<nom <<" +- "<<err_correct<<endl;

            tmph=TH1F(TString("Integral_TotModelMu0_uncondfit")+tag,TString("Integral_TotModelMu0_uncondfit")+tag, 1,0,1);
            tmph.SetBinContent(1,nom);
            tmph.SetBinError(1,err_correct);
            tmph.Write();

            pdftmp->plotOn(frame,FillColor(kOrange),LineWidth(2),LineColor(kBlue),VisualizeError(*fitresGlobal,1),
                           Normalization(1,RooAbsReal::RelativeExpected) ,Name("FitError_AfterFit_Mu0"));


        }

        //Create histo of data
        TString histName = modelName+"_postfit";
        histName = modelName+"_data";
        TH1* h_data = datatmp->createHistogram(histName,*obs);
        if(convertaxis)
            h_data= ConvertBDTAxis(h_data);

        myfile<<"Data = "<<h_data->Integral() <<endl;
        myfile.close();

        // Now prefit:
        w->loadSnapshot("snapshot_paramsVals_initial");
        if (poi)
            poi->setVal(1.);
        float preFitIntegral=0;
        float preFitIntegral_comp=0;

        funcIter = funcList.iterator() ;
        while( (comp = (RooProduct*) funcIter.Next()) ) {
            cout << "\t" << comp->GetName() << "\t" << (comp->createIntegral(*obs))->getVal() * binWidth->getVal() << endl;
            preFitIntegral_comp+=( comp->createIntegral(*obs))->getVal() * binWidth->getVal() ;
            if (poi)
                poi->setVal(1.0);

            // Plot all components also before fit!
            double Ntemp=(comp->createIntegral(*obs))->getVal() * binWidth->getVal();
            TString compname(comp->GetName());
            compname.ReplaceAll("L_x_","");
            compname.ReplaceAll(tt->GetName(),"");
            compname.ReplaceAll("_overallSyst_x_StatUncert","");
            compname.ReplaceAll("_overallSyst_x_HistSyst","");
            compname.ReplaceAll("_overallSyst_x_Exp","");

            pdfmodel->plotOn(frame,LineWidth(0),Components(*comp),LineColor(0), LineStyle(0), Normalization(Ntemp,RooAbsReal::NumEvent),Name("NoStacked_BeforeFitMu1_"+compname));

        }
        // Set mu to 0 for beforefitbkd
        if (poi)
            poi->setVal(0.0);
        preFitIntegral = pdftmp->expectedEvents(*obs);
        pdftmp->plotOn(frame,LineWidth(2),Name("BeforeFit_BkgTot"),LineStyle(kDashed), Normalization(preFitIntegral,RooAbsReal::NumEvent));
        //Draw the frame
        c2->cd();
        frame->Draw();

        vector<TH1*> vec_histo;
        TGraphAsymmErrors *gdata;
        ExtractHistoFromCanvas(c2,vec_histo,gdata);

        file->cd();
        h_data->Write(TString("DATA")+=tag);
        for(int iH=0;iH<vec_histo.size();iH++){

            if(convertaxis)
                vec_histo[iH]= ConvertBDTAxis(vec_histo[iH]);

            if (tag!=""){
                vec_histo[iH]->Write( vec_histo[iH]->GetName() + tag);
            }
            else
                vec_histo[iH]->Write();
        }
        file->Close();

    } // end while(tt) (loop on categories)
    cout<<"DumpPostFitHistos() finishing at "<<TTimeStamp().AsString()<<endl;
    return;
}


void TransportNPVals(RooArgSet* inSet, RooArgSet* outSet){
    /*
      Loops trough parameters in inSet and applies the values to vars in outSet
    */

    TIterator* InIter = inSet->createIterator();
    RooRealVar* var;
    while ((var = (RooRealVar*)InIter->Next())){
        TString name =  var->GetName();
        //    cout<<name<<endl;
        if(name.Contains("binWidth_obs_x_"))continue;
        if(name.Contains("gamma_stat"))continue;
        if(name.Contains("nom_"))continue;

        RooRealVar* outVar = (RooRealVar*) outSet->find(name);
        if (!outVar){
            cout<<"Warning: Was not able to find "<<name<<" in outSet."<<endl;
            continue;
        }

        outVar->setVal(var->getVal());
        outVar->setAsymError(var->getErrorLo(),var->getErrorHi() );
        //    cout<<"Setting "<<name<<" to "<< outVar->getVal()<<endl;
    }

    return;
}



RooFitResult* FitPDF( ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata,
                      int &MinuitStatus, int &HessStatus, double &Edm,
                      TString minimType, bool useMinos ) {
    cout<<model<<endl;
    model->Print();

    RooArgSet* constrainedParams = fitpdf->getParameters(*fitdata);
    RemoveConstantParameters(constrainedParams);
    Constrain(*constrainedParams);

    const RooArgSet* glbObs = model->GetGlobalObservables();

    RooRealVar * poi = (RooRealVar*) model->GetParametersOfInterest()->first();
    cout << "Constatnt POI " << poi->isConstant() << endl;
    cout << "Value of POI  " << poi->getVal() << endl;

    RooAbsReal * nll = fitpdf->createNLL(*fitdata, Constrain(*constrainedParams), GlobalObservables(*glbObs), Offset(1) );
    nll->enableOffsetting(true);

    double nllval = nll->getVal();

    std::cout << "initial parameters" << std::endl;
    constrainedParams->Print("v");

    std::cout << "INITIAL NLL = " << nllval << std::endl;

    static int nrItr = 0;
    int maxRetries = 3;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimType);
    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    int save_strat = strat;
    RooMinimizer minim(*nll);
    minim.setStrategy(strat);
    minim.setPrintLevel(1);
    minim.setEps(1);

    TStopwatch sw; sw.Start();

    int status=-99;
    HessStatus=-99;
    Edm = -99;
    RooFitResult * r;
    while (nrItr<maxRetries && status!=0 && status!=1){

        cout << endl;
        cout << endl;
        cout << endl;
        cout << "Fit try n°" << nrItr+1 << endl;
        cout << "======================" << endl;
        cout << endl;


        ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(),ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        HessStatus= minim.hesse();
        r = minim.save();
        Edm = r->edm();

        //up the strategy
        bool FitIsNotGood = ((status!=0 && status!=1) || (HessStatus!=0 && HessStatus!=1) || Edm>1.0);
        if (FitIsNotGood && strat<2){
            cout << endl;
            cout << "   *******************************" << endl;
            cout << "   * Increasing Minuit strategy (was " << strat << ")" << endl;
            strat++;
            cout << "   * Fit failed with : " << endl;
            cout << "      - minuit status " << status << endl;
            cout << "      - hess status " << HessStatus << endl;
            cout << "      - Edm = " << Edm << endl;
            cout << "   * Retrying with strategy " << strat << endl;
            cout << "   ********************************" << endl;
            cout << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            HessStatus= minim.hesse();
            r = minim.save();
            Edm = r->edm();
        }

        FitIsNotGood = ((status!=0 && status!=1) || (HessStatus!=0 && HessStatus!=1) || Edm>1.0);
        if (FitIsNotGood && strat < 2){
            cout << endl;
            cout << "   ********************************" << endl;
            cout << "   * Increasing Minuit strategy (was " << strat << ")" << endl;
            strat++;
            cout << "   * Fit failed with : " << endl;
            cout << "      - minuit status " << status << endl;
            cout << "      - hess status " << HessStatus << endl;
            cout << "      - Edm = " << Edm << endl;
            cout << "   * Retrying with strategy " << strat << endl;
            cout << "   ********************************" << endl;
            cout << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            HessStatus= minim.hesse();
            r = minim.save();
            Edm = r->edm();
        }

        FitIsNotGood = ((status!=0 && status!=1) || (HessStatus!=0 && HessStatus!=1) || Edm>1.0);
        if (FitIsNotGood && strat < 2){
            cout << endl;
            cout << "   *******************************" << endl;
            cout << "   * Increasing Minuit strategy (was " << strat << ")" << endl;
            strat++;
            cout << "   * Fit failed with : " << endl;
            cout << "      - minuit status " << status << endl;
            cout << "      - hess status " << HessStatus << endl;
            cout << "      - Edm = " << Edm << endl;
            cout << "   * Retrying with strategy " << strat << endl;
            cout << "   ********************************" << endl;
            cout << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
            HessStatus= minim.hesse();
            r = minim.save();
            Edm = r->edm();
        }

        if(useMinos) { minim.minos(); }

        // To run Minos only on mu
        RooArgSet *MyArgSet = (RooArgSet*) model->GetParametersOfInterest();
        // Add all mus:
        RooArgSet* paramsIn = (RooArgSet*) fitpdf->getParameters(*fitdata) ;
        TIterator* InIter = paramsIn->createIterator();
        RooRealVar* var;
        while ((var = (RooRealVar*)InIter->Next())){
            TString name =  var->GetName();
            if( name.Contains("SigXsecOverSM") || name.Contains("mu_XS8_VBF")|| name.Contains("mu_XS8_ggH") ){
                MyArgSet->add(*var);
                cout<<"Adding to minos set: "<< name<<endl;
            }
        }
        minim.minos( *MyArgSet );

        FitIsNotGood = ((status!=0 && status!=1) || (HessStatus!=0 && HessStatus!=1) || Edm>1.0);
        if ( FitIsNotGood) nrItr++;
        if (nrItr == maxRetries) {
            cout << endl;
            cout << endl;
            cout << endl;
            cout << "***********************************************************" << endl;
            cout << "WARNING::Fit failure unresolved with status " << status << endl;
            cout << "   Please investigate your workspace" << endl;
            cout << "   Find a wall : you will need it to crash your head on it" << endl;
            cout << "***********************************************************" << endl;
            cout << endl;
            cout << endl;
            cout << endl;
            MinuitStatus = status;

            return r;
        }

    }

    r = minim.save();
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "***********************************************************" << endl;
    cout << "         FIT FINALIZED SUCCESSFULLY : " << endl;
    cout << "            - minuit status " << status << endl;
    cout << "            - hess status " << HessStatus << endl;
    cout << "            - Edm = " << Edm << endl;
    cout << "            - NLL = " << setprecision(6)<<setiosflags(ios::fixed)<< nll->getVal() <<endl;
    cout << "            - MU = "<< ( (RooRealVar*) model->GetParametersOfInterest()->first())->getVal()<<" + "<< ( (RooRealVar*) model->GetParametersOfInterest()->first())->getErrorHi()<<" - " << ( (RooRealVar*) model->GetParametersOfInterest()->first())->getErrorLo()<<endl;
    cout << " -- " ; sw.Print();
    cout << "***********************************************************" << endl;
    cout << endl;
    cout << endl;
    cout << endl;

    MinuitStatus = status;
    sw.Print();


    return r;
} // FitPDF



void ExtractHistoFromCanvas(TPad *can, vector<TH1*> &vec_histo, TGraphAsymmErrors* &gdataAsymErr){

    vec_histo.clear();

    vector<TObject*> vec_MyObjects;
    TList *ListOfObject = can->GetListOfPrimitives();
    if (ListOfObject) {
        TIter next(ListOfObject);
        TObject *ThisObject;
        vec_MyObjects.clear();
        while ( (ThisObject=(TObject*)next()) ) vec_MyObjects.push_back( ThisObject->Clone() );
    }

    TH1D *hmodel(0);
    for (unsigned i=0 ; i<vec_MyObjects.size() ; i++){
        TString ObjectName = vec_MyObjects[i]->GetName();
        TString ClassName  = vec_MyObjects[i]->ClassName();

        if (ObjectName.Contains("NotAppears")) continue;
        if (ObjectName.Contains("Stacked_") && !ObjectName.Contains("NoStacked_")  ) continue;

        if(ClassName=="TH1D" && !hmodel){
            hmodel = (TH1D*) can->FindObject(ObjectName);
            if(hmodel){
                cout<<"Found a histogram: name="<<ObjectName<<endl;
                hmodel->Print("ranges");
            }
        }
        if (!hmodel) continue;

        TH1D *htemp(0);
        if(ClassName=="RooHist"){
            TGraphAsymmErrors *gRooHist = (TGraphAsymmErrors*) can->FindObject(ObjectName);
            if (!gRooHist) cout << "TGraphAsymmErrors object for " << ObjectName << " (class : " << ClassName << ") was not found !" << endl;
            if ( ObjectName=="Data" ){
                gdataAsymErr = (TGraphAsymmErrors*) gRooHist->Clone("DataAsymError");
                gdataAsymErr->SetTitle("DataAsymError");
                gdataAsymErr->SetLineWidth(2);
                gdataAsymErr->SetLineColor(1);
                gdataAsymErr->SetLineStyle(1);
            }
            htemp = TGraphAsymmErrorsToTH1D(gRooHist,hmodel);
            can->cd();
            gRooHist->Draw("LP");

        }

        if(ClassName=="RooCurve"){
            TGraph *gRooCurve = (TGraph*) can->FindObject(ObjectName);
            if (!gRooCurve) cout << "TGraph object for " << ObjectName << " (class : " << ClassName << ") was not found !" << endl;
            htemp = TGraphToTH1D(gRooCurve,hmodel);
        }

        if (htemp) vec_histo.push_back(htemp);

    }

    return;
}

TH1D* TGraphToTH1D(TGraph *g, TH1D *hmodel){

    TString gname = (TString)g->GetName();
    bool IsForError=false;
    if (gname.Contains("FitError_AfterFit")) IsForError=true;
    //  cout<<"What is in this hmodel?"<<endl;
    //hmodel->Print();
    // cout<<"call a clone?"<<endl;
    TH1D* hres = (TH1D*) hmodel->Clone();
    hres->Reset();
    TString hnameOld = (TString)g->GetName();
    cout<<"calling TGrapgToTH1D on"<<hnameOld<<endl;
    //  if (hnameOld.Contains("BeforeFit") && !hnameOld.Contains("BkgBeforeFit")) hnameOld.ReplaceAll("BeforeFit","BeforeFit_BkgTot");
    hnameOld.ReplaceAll("NoStacked_","");
    hnameOld.ReplaceAll("BkgBeforeFit_","BeforeFit_");
    hres->SetName("My_"+hnameOld);
    hres->SetTitle(hnameOld);
    hres->SetLineWidth( 2 );
    hres->SetLineColor( 1 );
    hres->SetLineStyle( 1 );

    // Get TH1F of the errors
    if (IsForError){
        cout<<"Is for Error , N="<<g->GetN() <<endl;
        //    g->Print();
        double xl,yl,xh,yh;
        xl=yl=xh=yh=0;
        for(int i=0;i<g->GetN();i++) {
            g->GetPoint(i,xl,yl);
            //cout<<" "<<i<<":"<<xl<<","<<yl<<endl;
        }


        int NN=g->GetN() / 2;

        TGraph* gUp=new TGraph(NN);
        TGraph* gLow=new TGraph(NN);

        vector<double> ratios;
        for(int i=0;i<NN;i++) {
            g->GetPoint(i,xl,yl);
            g->GetPoint(i+NN,xh,yh);
            gUp->SetPoint(NN-1-i,xh,yh);
            gLow->SetPoint(i,xl,yl);

            // cout<<"setting gup ("<<NN-1-i<<") :  "<<xh<<" "<<yh<<endl;
            // cout<<"setting glow("<<i<<") :  "<<xl<<" "<<yl<<endl;

        }

        int Nbins = hmodel->GetNbinsX();
        double xc,yc,xu,yu,xl2,yl2;
        for(int i=0;i<NN;i++) {
            gUp->GetPoint(i,xu,yu);
            gLow->GetPoint(i,xl2,yl2);
            //  cout<<" xu="<<xu<<" xl2="<<xl2<<endl;
            //xc = (xu - xl2)/2.;
            for (int ib=0 ; ib<Nbins+1 ; ib++){
                double dx=hmodel->GetBinWidth(ib);
                double xbin=hmodel->GetBinCenter(ib);
                xc=xl2+dx/2.0;
                //	cout<<"  ..My guess for the bincenter of the tgraph is:"<<xc<<endl;
                //cout<<"  xl="<<xl2<<" xu="<<xu<<endl;

                yc=(yu+yl2)/2.0;
                if (fabs(xc-xbin)<dx/10.){
                    // cout<<"Decided on bin. xc="<<xc<<" xbin="<<xbin<<endl;
                    //hres->SetBinContent(ib+1,yc);
                    //hres->SetBinError(ib+1,fabs((yu-yl2)/2.0));
                    hres->SetBinContent(ib,yc);
                    hres->SetBinError(ib,fabs((yu-yl2)/2.0));
                }
            }
        }
        // cout<<"will return something"<<endl;

        return hres;
    }


    // Get the TH1F for central value
    int Npts = g->GetN();
    double xprevious=0.;
    for (int ip=0 ; ip<Npts ; ip++){

        double x=0.0;
        double y=0.0;
        g->GetPoint(ip,x,y);
        if (ip==0) xprevious=x;
        // if( ip==0 && y==0){ cout<<"I suspect black magic. go away"<<endl; continue;}
        cout<<" point:"<<ip<<" x="<<x<<" y="<<y<<" xprevious="<<xprevious<<endl;
        //Nils just commented out to get an overview..
        if ( fabs(xprevious-x)>hres->GetBinWidth(1)/100. ){
            double xbin = (x+xprevious)/2. ;
            double ycontent = y;

            int ib = hres->FindBin(xbin);
            cout<<"...decided on xposition:"<<xbin<<" corr.binnr="<<ib<<endl;
            if (ib<=hres->GetNbinsX()) hres->SetBinContent(ib,ycontent);
        }
        xprevious = x;
    }

    return hres;
}


TH1D* TGraphAsymmErrorsToTH1D(TGraphAsymmErrors *g, TH1D *hmodel){
    TH1D* hres = (TH1D*) hmodel->Clone();
    hres->Reset();
    TString hnameOld = (TString)g->GetName();
    cout<<"Calling TGraphAsymmErrorsToTH1D on"<<hnameOld<<endl;
    if (hnameOld.Contains("BeforeFit") && !hnameOld.Contains("BkgBeforeFit")) hnameOld.ReplaceAll("BeforeFit","BeforeFit_BkgTot");
    hnameOld.ReplaceAll("NoStacked_","");
    hnameOld.ReplaceAll("BkgBeforeFit_","BeforeFit_");
    hres->SetName("My_"+hnameOld);
    hres->SetTitle(hnameOld);
    hres->SetLineWidth( 2 );
    hres->SetLineColor( 1 );
    hres->SetLineStyle( 1 );


    int Npts = g->GetN();
    double xprevious = 0.0;
    for (int ip=0 ; ip<Npts ; ip++){
        double x,y,erry;
        g->GetPoint(ip,x,y);
        erry = g->GetErrorY(ip);
        if (ip==0) xprevious=x;
        if (  fabs(xprevious-x)>hres->GetBinWidth(1)/100.  ){
            double xbin = (x+xprevious)/2. ;
            //      if (ip==0) xbin=x;
            //xbin = x;
            cout<<"xbin="<<xbin<<" x="<<x<<" xprevious="<<xprevious<<endl;
            double ycontent = y;
            int ib = hres->FindBin(xbin);
            if (ib<=hres->GetNbinsX()) {
                hres->SetBinContent(ib,ycontent);
                hres->SetBinError(ib,erry);
            }
        }
        xprevious = x;
    }

    TCanvas *can = new TCanvas(hres->GetName(),hres->GetName(),600,600);
    can->cd();
    hres->Draw();
    return hres;
}



/* Converts the histogram to have a axis between -1 and 1*/

TH1* ConvertBDTAxis(TH1 *hist){
    TH1D*  newh=new TH1D(TString(hist->GetName())+"newaxis",TString(hist->GetName())+"newaxis",hist->GetNbinsX(),-1,1);
    for(int b=0; b<=hist->GetNbinsX()+1; b++ ){
        newh->SetBinContent(b, hist->GetBinContent(b));
        newh->SetBinError(b, hist->GetBinError(b));
    }
    return newh;
}
//----------------------------------------------------------
int main(int argc, char** argv) {
    string inputWorkspace;
    string workspaceVariable; //only used if you want to apply the fit results to another distribution
    string outputFile;
    // keep the same defaults as function def (might need to update)
    bool doFit = true;
    bool transportCovariance = false;
    bool justFit = false;
    bool convertAxis = false;
    bool fixSystematicParameters = false;

    int optind(1);
    while ((optind < argc)) {
        std::string sw = argv[optind];
        if(sw[0]!='-') { /*if(dbg)*/ cout<<"skip "<<sw<<endl; optind++; continue; }
        if     (sw=="-i"||sw=="--input") { inputWorkspace = argv[++optind]; }
        else if(sw=="-V"||sw=="--variable" ) { workspaceVariable = argv[++optind]; }
        else if(sw=="-o"||sw=="--output") { outputFile = argv[++optind]; }
        else if(sw=="--skip-fit") { doFit = false; }
        else if(sw=="--transport-cov") { transportCovariance = true; }
        else if(sw=="--just-fit") { justFit = true; }
        else if(sw=="--convert-axis") { convertAxis = true; }
        else if(sw=="--constant-syst") { fixSystematicParameters = true; }
        // \todo else if(sw=="-h"||sw=="--help" ) { usage(argv[0], matrixFile.c_str()); return 0; }
        else cout<<"Unknown switch "<<sw<<endl;
        optind++;
    } // end while(optind<argc)

    string snapshotname = "";
    DumpPostFitHistos(inputWorkspace.c_str(), workspaceVariable.c_str(), outputFile.c_str(), snapshotname,
                      doFit, transportCovariance, justFit, convertAxis, fixSystematicParameters);
                      // true, true, false, false);
    return 0;
}
