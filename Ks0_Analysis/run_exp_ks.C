#wrapper to run cuts.C in ACLIC mode. (much faster)
{
   gROOT->Reset();
   gROOT->LoadMacro("v0dst.C+");
   gROOT->LoadMacro("cuts_exp_ks.C+");
   cuts_exp_ks(0,100);
}

