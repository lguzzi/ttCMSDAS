import ROOT

# Global variables
intLumi = 296.08 # pb^-1, uncertainty = 3.5 %

#Channels: ee, mumu, emu
channelString = "temp"

#Open file to read yield histograms
f_input = ROOT.TFile("/gpfs/ddn/cms/user/jlangfor/top_quark_exercise/CMSSW_10_2_10/src/ttCMSDAS/ttbar/%s/TT.root"%channelString)

#Define systematic dictionary
syst_dict = {"TT_muonID":"Y_SFmuon", "TT_electronID":"Y_SFelec"}

#Extract nominal events from first histogram
h_initialize = f_input.Get( syst_dict["TT_muonID"] )
n_nominal = h_initialize.GetBinContent(1)
n_nominal_statUnc = h_initialize.GetBinError(1)

print ""
print "N_signal = %5.4f +- %5.4f"%(n_nominal*intLumi,n_nominal_statUnc*intLumi)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Iterate over keys, values in dictionary and extract the yield variations
for syst, hSyst_string in syst_dict.iteritems():

  #Initiate string to hold yield variations for systematic
  syst_string = "%s "%syst

  #Get the corresponding histogram
  h_syst = f_input.Get("%s"%hSyst_string)

  #For "specific" systematics: treat differently
  if( syst == "TT_QCDScale" ): continue
  else:

    #Get number of events with up and down variation
    n_up = h_syst.GetBinContent(2)
    n_down = h_syst.GetBinContent(3)
  
    #Yield variations (percentage)
    yieldVar_up = n_up/n_nominal
    yieldVar_down = n_down/n_nominal
    yieldVar_central = (0.5*abs(n_up-n_down))/n_nominal

    #Add information to string
    syst_string += "%5.4f/%5.4f (%4.2f%%/%4.2f%%) (central = %4.2f%%)"%(yieldVar_up,yieldVar_down,abs(yieldVar_up-1)*100,abs(yieldVar_down-1)*100,(yieldVar_central*100))

    print syst_string
    
print ""
