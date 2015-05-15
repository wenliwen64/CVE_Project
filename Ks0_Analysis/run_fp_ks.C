#wrapper to run cuts.C in ACLIC mode. (much faster)
{
   gROOT->Reset();
   gROOT->LoadMacro("v0dst.C+");
   gROOT->LoadMacro("cuts_fp_ks.C+");
   cuts_fp_ks(0,100);
}

