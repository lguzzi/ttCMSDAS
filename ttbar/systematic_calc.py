import ROOT

#Channels: ee, mumu, emu
channelString = "temp"

#Open file to read yield histograms
f_input = ROOT.TFile("/gpfs/ddn/cms/user/jlangfor/top_quark_exercise/CMSSW_10_2_10/src/ttCMSDAS/ttbar/%s/TT.root"%channelString)

#Define systematic dictionary
syst_dict = {"TT_muonID":"Y_SFmuon", "TT_electronID":"Y_SFelectron"}

#Extract nominal events from first histogram
h_initialize = f_input.Get( syst_dict["TT_muonID"] )
n_nominal = h_initialize.GetBinContent(1)
n_nominal_statUnc = h_initialize.GetBinError(1)

print "N_signal = %4.2f +- %4.2f"%(n_nominal,n_nominal_statUnc)

#Iterate over keys, values in dictionary and extract the yield variations
for syst, hSyst_string in syst_dict.iteritems():

  #Initiate string to hold yield variations for systematic
  syst_string = "%s "%syst_dict

  #Get the corresponding histogram
  h_syst = f_input.Get("%s"%hSyst_string)

  #For "specific" systematics: treat differently
  if( syst = "TT_QCDScale" ): continue
  else:

    #Get number of events with up and down variation
    n_up = h_syst.GetBinContent(2)
    n_down = h_syst.GetBinContent(3)
  
    #Yield variations (percentage)
    yieldVar_up = n_up/n_nominal
    yieldVar_down = n_down/n_nominal
    yieldVar_cental = (0.5*abs(n_up-n_down))/n_nominal

    #Add information to string
    syst_string += "%4.2f/%4.2f (%4.2f)"%(yieldVar_up,yieldVar_down,yieldVar_central)

    print syst_string
    
