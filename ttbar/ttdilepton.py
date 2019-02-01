import os,sys
basepath = os.path.abspath(__file__).rsplit('/ttCMSDAS/',1)[0]+'/ttCMSDAS/'
   
   
sys.path.append(basepath)

from framework.analysis import analysis
from framework.functions import DeltaPhi, DiPt, InvMass, lepton, jet
from ROOT.TMath import Sqrt as sqrt
from ROOT import *

################ Analysis
class ttdilepton(analysis):
  def init(self):
    # Load SF files
    if not self.isData:
      self.LoadHisto('MuonIsoSF', basepath+'./inputs/MuonISO.root', 'NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta') # pt, abseta
      self.LoadHisto('MuonIdSF',  basepath+'./inputs/MuonID.root',  'NUM_TightID_DEN_genTracks_pt_abseta') # pt, abseta
      self.LoadHisto('ElecSF',    basepath+'./inputs/ElecTightCBid94X.root',  'EGamma_SF2D') # eta, pt

    # Objects for the analysis
    self.selLeptons = []
    self.selJets = []
    
    self.selJets_jes_up = []
    self.selJets_jes_do = []
    self.selJets_jer_up = []
    self.selJets_jer_do = []
   
    self.pmet = TLorentzVector()

    # Create output histograms
    self.CreateTH1F("Lep0Pt",   "", 24, 0, 120)
    self.CreateTH1F("Lep1Pt",   "", 24, 0, 120)
    self.CreateTH1F("Lep0Eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("Lep1Eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("InvMass",  "", 60, 0, 300)
    self.CreateTH1F("DilepPt",  "", 40, 0, 200)
    self.CreateTH1F("DeltaPhi", "", 20, 0, 1)
 
    self.CreateTH1F('Ymm_yield', "", 1, 0, 1)
    self.CreateTH1F('Yme_yield', "", 1, 0, 1)
    self.CreateTH1F('Yee_yield', "", 1, 0, 1)
    
    self.CreateTH1F("Yee_MatrixEl", "", 9, 0, 9) 
    self.CreateTH1F("Yee_SFmuon", "", 3, 0, 3)
    self.CreateTH1F("Yee_SFelec", "", 3, 0, 3)
    self.CreateTH1F('Yee_SFPU'  , "", 3, 0, 3)
    self.CreateTH1F("Yee_Pdf", "", 31, 0, 31) 
    self.CreateTH1F("Yee_Alphas", "", 3, 0, 3) 
 
    self.CreateTH1F("Ymm_MatrixEl", "", 9, 0, 9) 
    self.CreateTH1F("Ymm_SFmuon", "", 3, 0, 3)
    self.CreateTH1F("Ymm_SFelec", "", 3, 0, 3)
    self.CreateTH1F('Ymm_SFPU'  , "", 3, 0, 3)
    self.CreateTH1F("Ymm_Pdf", "", 31, 0, 31) 
    self.CreateTH1F("Ymm_Alphas", "", 3, 0, 3) 
  
    self.CreateTH1F("Yme_MatrixEl", "", 9, 0, 9) 
    self.CreateTH1F("Yme_SFmuon", "", 3, 0, 3)
    self.CreateTH1F("Yme_SFelec", "", 3, 0, 3)
    self.CreateTH1F('Yme_SFPU'  , "", 3, 0, 3)
    self.CreateTH1F("Yme_Pdf", "", 31, 0, 31) 
    self.CreateTH1F("Yme_Alphas", "", 3, 0, 3)
    
    self.CreateTH1F("Yme_jetS", "", 3, 0, 3)
    self.CreateTH1F("Yee_jetS", "", 3, 0, 3)
    self.CreateTH1F("Ymm_jetS", "", 3, 0, 3)
    
    
    self.CreateTH1F("Yee_jetR", "", 3, 0, 3)
    self.CreateTH1F("Ymm_jetR", "", 3, 0, 3)
    self.CreateTH1F("Yme_jetR", "", 3, 0, 3)
    
    
  

  def resetObjects(self):
    ''' Reset the list where the objects are stored '''
    self.selLeptons = []
    self.selJets = []
    self.selJets_jes_do = []
    self.selJets_jes_up = []
    self.selJets_jer_up = []
    self.selJets_jer_do = []
    self.pmet = TLorentzVector()

  def FillHistograms(self, leptons, jets, pmet, flavour):
    ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
    
    if flavour == 0: flav = 'mm'
    if flavour == 1: flav = 'ee'
    if flavour == 2: flav = 'me'


    if not len(leptons) >= 2: return # Just in case

    # nominal
    self.weight = self.EventWeight * self.SFmuon * self.SFelec * self.PUSF
    
    # muon Up and Down
    self.weightSFmuonUp = self.EventWeight * (self.SFmuon + self.SFmuonErr) * self.SFelec * self.PUSF
    self.weightSFmuonDown = self.EventWeight * (self.SFmuon - self.SFmuonErr) * self.SFelec * self.PUSF
    
    # electrons Up and Down
    self.weightSFelecUp = self.EventWeight * self.SFmuon * (self.SFelec + self.SFelecErr) * self.PUSF
    self.weightSFelecDown = self.EventWeight * self.SFmuon * (self.SFelec - self.SFelecErr) * self.PUSF

    ## PU Up and Down
    self.weightSFPUUp   =  self.weight * self.PUUpSF / self.PUSF
    self.weightSFPUDown =  self.weight * self.PUDoSF / self.PUSF

    # Re-calculate the observables
    lep0  = leptons[0]; lep1 = leptons[1]
    l0pt  = lep0.Pt();  l1pt  = lep1.Pt()
    l0eta = lep0.Eta(); l1eta = lep1.Eta()
    dphi  = DeltaPhi(lep0, lep1)
    mll   = InvMass(lep0, lep1)
    dipt  = DiPt(lep0, lep1)
    
    ## JES and JER on top of everything
    if len(self.selJets       ) >= 2: self.obj["Y%s_jetS" %flav].Fill(0.5, self.weight)
    if len(self.selJets       ) >= 2: self.obj["Y%s_jetR" %flav].Fill(0.5, self.weight)

    if len(self.selJets_jes_up) >= 2: self.obj["Y%s_jetS" %flav].Fill(1.5, self.weight)
    if len(self.selJets_jes_do) >= 2: self.obj["Y%s_jetS" %flav].Fill(2.5, self.weight)
    
    if len(self.selJets_jer_up) >= 2: self.obj["Y%s_jetR" %flav].Fill(1.5, self.weight)
    if len(self.selJets_jer_do) >= 2: self.obj["Y%s_jetR" %flav].Fill(2.5, self.weight)
    
    ## From now on selections should be the nominal ones, and only ket number was missing
    ## NOTE should we use jet id?
    if len(self.selJets) < 2: return
    ### Fill the histograms
    self.obj['Lep0Pt'].Fill(l0pt, self.weight)
    self.obj['Lep1Pt'].Fill(l1pt, self.weight)
    self.obj['Lep0Eta'].Fill(l0eta, self.weight)
    self.obj['Lep1Eta'].Fill(l1eta, self.weight)
    self.obj["InvMass"].Fill(mll, self.weight)
    self.obj['DilepPt'].Fill(dipt, self.weight)
    self.obj['DeltaPhi'].Fill(dphi/3.141592, self.weight)
    
    self.obj['Y%s_SFmuon' %flav].Fill(0.5, self.weight)
    self.obj['Y%s_SFmuon' %flav].Fill(1.5, self.weightSFmuonUp)
    self.obj['Y%s_SFmuon' %flav].Fill(2.5, self.weightSFmuonDown)
    
    self.obj['Y%s_SFelec' %flav].Fill(0.5, self.weight)
    self.obj['Y%s_SFelec' %flav].Fill(1.5, self.weightSFelecUp)
    self.obj['Y%s_SFelec' %flav].Fill(2.5, self.weightSFelecDown)

    self.obj['Y%s_SFPU' %flav ].Fill(0.5, self.weight        )
    self.obj['Y%s_SFPU' %flav ].Fill(1.5, self.weightSFPUUp  )
    self.obj['Y%s_SFPU' %flav ].Fill(2.5, self.weightSFPUDown)

    #self.obj['Y%s_yield' %flav].Fill(0.5, self.weight)
 
    for ii in range(9):
        if ii == 6 or ii == 2: continue     ## skip the unphysical values (2-0.5 and 0.5, 2)
        self.obj['Y%s_MatrixEl' %flav].Fill(ii + 0.5, self.weight * self.SFlhe[ii])
    
    for ii in range(31):
        self.obj['Y%s_Pdf' %flav].Fill(ii + 0.5, self.weight * self.LHEPdf[ii])
    
    for ii in range(3):
        self.obj['Y%s_Alphas' %flav].Fill(ii + 0.5, self.weight * self.LHEAlphas[ii])
    
  def insideLoop(self, t):
    self.resetObjects()

    ### Lepton selection
    ###########################################
    if not self.isData: nGenLep = t.nGenDressedLepton 
    
    ##### Muons
    for i in range(t.nMuon):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Muon_pt[i], t.Muon_eta[i], t.Muon_phi[i], t.Muon_mass[i])
      charge = t.Muon_charge[i]

      # Tight ID, tight ISO, RelIso04 < 0.15, tight IP
      if not t.Muon_tightId[i]: continue # Tight ID
      if not t.Muon_pfRelIso04_all[i] < 0.15: continue
      dxy = abs(t.Muon_dxy[i])
      dz  = abs(t.Muon_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 12 GeV, |eta| < 2.4
      if p.Pt() < 12 or abs(p.Eta()) > 2.4: continue
      self.selLeptons.append(lepton(p, charge, 13))
       
    ##### Electrons
    for i in range(t.nElectron):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Electron_pt[i], t.Electron_eta[i], t.Electron_phi[i], t.Electron_mass[i])
      charge = t.Electron_charge[i]
      etaSC    = abs(p.Eta());
      convVeto = t.Electron_convVeto[i]

      # Tight cut-based Id, convVeto, RelIso03 tight, tight IP
      if not t.Electron_cutBased[i] >= 4: continue
      if not convVeto: continue
      relIso03 = t.Electron_pfRelIso03_all[i]
      if   etaSC <= 1.479 and relIso03 > 0.0361: continue
      elif etaSC >  1.479 and relIso03 > 0.094:  continue
      dxy = abs(t.Electron_dxy[i])
      dz  = abs(t.Electron_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 12 GeV, |eta| < 2.4
      if p.Pt() < 12 or abs(p.Eta()) > 2.4: continue
      self.selLeptons.append(lepton(p, charge, 11))

    leps = self.selLeptons
    pts  = [lep.Pt() for lep in leps]
    self.selLeptons = [lep for _,lep in sorted(zip(pts,leps))]
    
    ###Jets
    ###########################################
    ntightJet=0
    nBJet=0
    MC=1
    for i in range(t.nJet):
      p = TLorentzVector()
      jes_pup = TLorentzVector()
      jes_pdo = TLorentzVector()
      jer_pup = TLorentzVector()
      jer_pdo = TLorentzVector()
      #p.SetPtEtaPhiM(t.Jet_pt[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass[i])
      p.SetPtEtaPhiM(t.Jet_pt[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass[i])
     
      jes_pup.SetPtEtaPhiM(t.Jet_pt_jesTotalUp[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass_jesTotalUp[i])
      jes_pdo.SetPtEtaPhiM(t.Jet_pt_jesTotalDown[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass_jesTotalDown[i])
      
      jer_pup.SetPtEtaPhiM(t.Jet_pt_jerUp[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass_jerUp[i])
      jer_pdo.SetPtEtaPhiM(t.Jet_pt_jerDown[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass_jerDown[i])
      
      if(self.isData==1):
        MC=-1
      else:
        MC=1
      if t.Jet_jetId[i] < 2: 
        continue #continue if it is not tight 
      if abs(p.Eta()) > 2.4: continue 
      if p.Pt() > 25             :       self.selJets.append(jet(p, t.Jet_btagCSVV2[i],MC,i,t.Jet_btagDeepB[i]))
      
      if jes_pup.Pt() > 25       :       self.selJets_jes_up.append(jet(jes_pup, t.Jet_btagCSVV2[i],MC,i,t.Jet_btagDeepB[i]))
      if jes_pdo.Pt() > 25       :       self.selJets_jes_do.append(jet(jes_pdo, t.Jet_btagCSVV2[i],MC,i,t.Jet_btagDeepB[i]))

      if jer_pup.Pt() > 25       :       self.selJets_jer_up.append(jet(jer_pup, t.Jet_btagCSVV2[i],MC,i,t.Jet_btagDeepB[i]))
      if jer_pdo.Pt() > 25       :       self.selJets_jer_do.append(jet(jer_pdo, t.Jet_btagCSVV2[i],MC,i,t.Jet_btagDeepB[i]))
   
     ###P Met
    self.pmet.SetPtEtaPhiM(t.MET_pt, 0, t.MET_phi, 0)
    
    ### Calculate the weights
    self.SFelec = 1; self.SFmuon = 1; self.SFelecErr = 0; self. SFmuonErr = 0
    self.weightSFmuonUp = 1 
    self.weightSFmuonDown = 1 
    self.weightSFelecUp = 1 
    self.weightSFelecDown = 1
    if not self.isData:
      for lep in self.selLeptons:
        if lep.IsMuon():
          sf, err = self.GetSFandErr('MuonIsoSF, MuonIdSF', lep.Pt(), TMath.Abs(lep.Eta()))
          self.SFmuon*=sf
          self.SFmuonErr+=err*err
        else:
          sf, err = self.GetSFandErr('ElecSF', lep.Eta(), lep.Pt())
          self.SFelec*=sf
          self.SFelecErr+=err*err
      self.SFelecErr = sqrt(self.SFelecErr)
      self.SFmuonErr = sqrt(self.SFmuonErr)

    # PU SF --> PLEASE CHECK THAT THE WEIGHTS ARE IN THE TREES THAT YOU'RE USING!
    if not self.isData:
      self.PUSF   = t.puWeight
      self.PUUpSF = t.puWeightUp
      self.PUDoSF = t.puWeightDown
    else:
      self.PUSF   = 1; self.PUUpSF = 1; self.PUDoSF = 1
    
    self.SFlhe = [1]*9
    #LHE weights
    if not self.isData:
        self.SFlhe = t.LHEScaleWeight

    self.LHEPdfAlphaWeight = [1]*33
    self.LHEPdf = [1]*31
    self.LHEAlphas = [1]*3
    #LHE PDF and alpha weights
    if not self.isData:
       self.LHEPdfAlphaWeight = t.LHEPdfWeight
       self.LHEPdf = [el for eli, el in enumerate(t.LHEPdfWeight) if eli < 31]
       self.LHEAlphas = [self.LHEPdfAlphaWeight[0]] + [el for eli, el in enumerate(t.LHEPdfWeight) if eli >= 31]

    ### Event selection
    ###########################################
    
    ### Dilepton pair: 2 leptons, opposite sign, mll > 20 GeV, leading lep pT > 20 GeV
    if not len(leps) >= 2:      return 
    l0 = leps[0]; l1 = leps[1]
    if l0.charge*l1.charge > 0: return 
    if l0.Pt() < 20:            return 
    if InvMass(l0,l1) < 20:     return  
    
    ##if t.nJet < 2             :   return
    ##if t.jet_jetId & 011 == 1 :   return  ## should be 0?   

    Z_MASS_PDG = 91.1876
    inv_mass = InvMass(l0, l1)

    ## HLT selection
    ##      simple HLT selection for doubleLep samples
    if    l0.IsMuon() and l1.IsMuon():
        flavour = 0
        ## DY offline cut
        if abs(inv_mass - Z_MASS_PDG) < 15:
            return
        if self.pmet.Pt() < 35:
            return

        ## online cut
        if not t.HLT_HIL3DoubleMu0: 
            return
    
    elif  l0.IsElec() and l1.IsElec():
        flavour = 1
        ## DY offline cut
        if abs(inv_mass - Z_MASS_PDG) < 15:
             return
        if self.pmet.Pt() < 35:
             return

        ## online cut
        if not t.HLT_HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ: 
            return     
    
    ##      depends on the sample for emu samples
    elif (l0.IsMuon() and l1.IsElec()) or \
         (l0.IsElec() and l1.IsMuon()) :
        
        flavour = 2
        ## online only
        if self.sampleName == 'HighEGJet':
            if not (t.HLT_HIEle20_WPLoose_Gsf and not t.HLT_HIL3Mu20):
                return
        if self.sampleName == 'SingleMuon':
            if not (t.HLT_HIL3Mu20):
                return
    
    else:
        print 'ERROR: 0 muons and 0 electrons were found'
        import pdb; pdb.set_trace()
 
    ### Fill the histograms
    self.FillHistograms(self.selLeptons, self.selJets, self.pmet, flavour)
