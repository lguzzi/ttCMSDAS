import ROOT

#Define function to calc signal strength: XSBR(data)/XSBR(expected)
def mu_calc( N_data, N_bkg, N_sig ): return (N_data-N_bkg)/N_sig

#Function to calculate uncertainty on mu
# N_sig = XSBR(theory)*eff*Acc*Lumi
# unc_sig = 
def mu_unc_calc( 
