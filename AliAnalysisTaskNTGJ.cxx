#include <climits>

// This is purely to suppress warning inside ROOT once instanciating
// std::vector<bool>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <TCollectionProxyInfo.h>
#pragma GCC diagnostic pop

#include <Compression.h>
#include <TFile.h>
#include <TSystem.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <AliVEvent.h>
#include <AliAODEvent.h>
#pragma GCC diagnostic pop
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <AliVCaloCells.h>
#include <AliOADBContainer.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>
#include <AliAODMCParticle.h>
#include <AliESDMuonTrack.h>

#include <AliMultSelection.h>
#include <AliEventplane.h>

#include <AliMagF.h>

#include <AliAODTrack.h>
#include <AliAODVertex.h>
#include <AliAODTracklets.h>
#include <AliAODMCHeader.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#include <AliAnalysisTaskEmcalEmbeddingHelper.h>
#pragma GCC diagnostic pop

#include "AliAnalysisTaskNTGJ.h"

#ifndef __CINT__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wdeprecated-register"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>

#include "eLut.cpp"
#include "half.cpp"

#include "keras_model.cc"

// #include "tiny_dnn/tiny_dnn.h"

#pragma GCC diagnostic pop
#include "special_function.h"
#include "mc_truth.h"
#include "emcal.h"
#include "jet.h"
#include "track_cuts.h"
#ifdef WITH_EFP7
#include "einstein_sum.h"
#include "efp7.cc"
#else // WITH_EFP7
#define FILL_EFP7 {}
#endif // WITH_EFP7
#include "isolation.h"
#endif // __CINT__

/////////////////////////////////////////////////////////////////////
// Replacement dlfcn.h for CINT

#ifndef RTLD_NOW
#define RTLD_NOW      0x00002
#endif // RTLD_NOW
#ifndef RTLD_GLOBAL
#define RTLD_GLOBAL   0x00100
#endif // RTLD_GLOBAL
#ifndef RTLD_NODELETE
#define RTLD_NODELETE 0x01000
#endif // RTLD_NODELETE

extern "C" {
extern void *dlopen(const char *, int);
extern void *dlsym(void *, const char *);
extern char *dlerror(void);
}

/////////////////////////////////////////////////////////////////////

// FIXME: This needs to be moved to somewhere else, later on

namespace {

    class alice_jec_t {
    protected:
        const char *_version;
    public:
        alice_jec_t(void)
            : _version("lhc16c2-2017-06-pre1")
        {
        }
        TLorentzVector jec_p(TLorentzVector p) const
        {
            const double pt = p.Pt();

            // Dummy value
            const double pt_calib = 2.0 * pt;

            const double eta = p.Eta();

            // Mean value for jet pT raw within 2.5-18 GeV
            static const double par_eta[] = { 0.015, 7.7, 0.65 };

            const double eta_calib = par_eta[0] *
                (TMath::Erf(par_eta[1] * (eta - par_eta[2])) +
                 TMath::Erf(par_eta[1] * (eta + par_eta[2])));

            TLorentzVector calib;

            calib.SetPtEtaPhiM(pt_calib, eta_calib, p.Phi(), p.M());

            return calib;
        }
        const char *version(void) const
        {
            return _version;
        }
    };

}

ClassImp(AliAnalysisTaskNTGJ);

#define EMCAL_GEOMETRY_NAME "EMCAL_COMPLETE12SMV1_DCAL_8SM"

#define ntrigger_class NTRIGGER_CLASS_MAX

#define B CHAR_MIN
#define b 0
#define S SHRT_MIN
#define s 0
#define I INT_MIN
#define i 0
#define F NAN
#define D NAN
#define L LONG_MIN
#define l 0
#define O false
#define BRANCH(b, t)                            \
    _branch_ ## b((t)),
#define BRANCH_ARRAY(b, d, t)
#define BRANCH_ARRAY2(b, d, e, t)
#define BRANCH_STR(b)
#define BRANCH_STR_ARRAY(b, d)                          \
    _branch_ ## b(std::vector<std::string>()),

#define CLASS_INITIALIZATION                                \
    _emcal_geometry_name(EMCAL_GEOMETRY_NAME),              \
    _tree_event(NULL),                                      \
    MEMBER_BRANCH                                           \
    _run_number_current(INT_MIN),                           \
    _f1_ncluster_tpc_linear_pt_dep(NULL),                   \
    _track_cut(std::vector<AliESDtrackCuts>()),             \
    _reco_util(new AliEMCALRecoUtils),                      \
    _emcal_geometry(NULL),                                  \
    _muon_track_cut(new AliMuonTrackCuts),                  \
    _ncell(EMCAL_NCELL),                                    \
    _emcal_geometry_filename(""),                           \
    _emcal_local2master_filename(""),                       \
    _force_ue_subtraction(false),                           \
    _skim_cluster_min_e(-INFINITY),                         \
    _skim_track_min_pt(-INFINITY),                          \
    _skim_muon_track_min_pt(-INFINITY),                     \
    _skim_jet_min_pt(std::vector<double>(3, -INFINITY)),    \
    _skim_jet_average_pt(-INFINITY),                        \
    _skim_multiplicity_tracklet_min_n(INT_MIN),             \
    _skim_sum_eg_ntrial(0),                                 \
    _stored_track_min_pt(-INFINITY),                        \
    _stored_jet_min_pt_raw(-INFINITY),                      \
    _nrandom_isolation(0),                                  \
    _is_embed(false),                                       \
    _emcal_mask(std::vector<bool>()),                       \
    _emcal_cell_position(NULL),                             \
    _keras_model_photon_discrimination(NULL),               \
    _alien_plugin(NULL),                                    \
    _metadata_filled(false)

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(const char *name)
    : AliAnalysisTaskEmcal(name), CLASS_INITIALIZATION
{
    DefineOutput(1, TTree::Class());
}

AliAnalysisTaskNTGJ::AliAnalysisTaskNTGJ(
    const AliAnalysisTaskNTGJ &x)
    : AliAnalysisTaskEmcal(), CLASS_INITIALIZATION
{
    _emcal_geometry_filename = x._emcal_geometry_filename;
    _emcal_local2master_filename = x._emcal_local2master_filename;

    _force_ue_subtraction = x._force_ue_subtraction;
    _skim_cluster_min_e = x._skim_cluster_min_e;
    _skim_track_min_pt = x._skim_track_min_pt;
    _skim_muon_track_min_pt = x._skim_muon_track_min_pt;
    _skim_jet_min_pt = x._skim_jet_min_pt;
    _skim_jet_average_pt = x._skim_jet_average_pt;
    _skim_multiplicity_tracklet_min_n =
        x._skim_multiplicity_tracklet_min_n;
    _stored_track_min_pt = x._stored_track_min_pt;
    _stored_jet_min_pt_raw = x._stored_jet_min_pt_raw;
    _nrandom_isolation = x._nrandom_isolation;
    _is_embed = x._is_embed;
}

#undef ntrigger_class

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY
#undef B
#undef b
#undef S
#undef s
#undef I
#undef i
#undef F
#undef D
#undef L
#undef l
#undef O

AliAnalysisTaskNTGJ &AliAnalysisTaskNTGJ::operator=(
    const AliAnalysisTaskNTGJ &x)
{
    if (this == &x) {
        return *this;
    }

    _emcal_geometry_filename = x._emcal_geometry_filename;
    _emcal_local2master_filename = x._emcal_local2master_filename;

    _force_ue_subtraction = x._force_ue_subtraction;
    _skim_cluster_min_e = x._skim_cluster_min_e;
    _skim_track_min_pt = x._skim_track_min_pt;
    _skim_muon_track_min_pt = x._skim_muon_track_min_pt;
    _skim_jet_min_pt = x._skim_jet_min_pt;
    _skim_jet_average_pt = x._skim_jet_average_pt;
    _skim_multiplicity_tracklet_min_n =
        x._skim_multiplicity_tracklet_min_n;
    _stored_track_min_pt = x._stored_track_min_pt;
    _stored_jet_min_pt_raw = x._stored_jet_min_pt_raw;
    _nrandom_isolation = x._nrandom_isolation;
    _is_embed = x._is_embed;

    return *this;
}

AliAnalysisTaskNTGJ::~AliAnalysisTaskNTGJ(void)
{
    if (_tree_event != NULL) {
        delete _tree_event;
    }
    if (_emcal_cell_position != NULL) {
        delete reinterpret_cast<std::vector<point_2d_t> *>
            (_emcal_cell_position);
    }
    if (_keras_model_photon_discrimination != NULL) {
        delete reinterpret_cast<KerasModel *>
            (_keras_model_photon_discrimination);
    }

    delete _reco_util;
}

void AliAnalysisTaskNTGJ::UserCreateOutputObjects(void)
{
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    // fAliEventCuts.fGreenLight = true;

    TFile *file = OpenFile(1);

    if (file != NULL) {
        file->SetCompressionSettings(ROOT::CompressionSettings(
                                         ROOT::kZLIB, 9));
    }

    /////////////////////////////////////////////////////////////////

    _tree_event = new TTree("_tree_event", "");

#define BRANCH(b, t)                    \
    _tree_event->Branch(                \
        #b, &_branch_ ## b, #b "/" #t);
#define BRANCH_ARRAY(b, d, t)                   \
    _tree_event->Branch(                        \
        #b, _branch_ ## b, #b "[" #d "]/" #t);
#define BRANCH_ARRAY2(b, d, e, t)                       \
    _tree_event->Branch(                                \
        #b, _branch_ ## b, #b "[" #d "][" #e "]/" #t);
#define BRANCH_STR(b)                           \
    _tree_event->Branch(                        \
        #b, _branch_ ## b, #b "/C");
#define BRANCH_STR_ARRAY(b, d)                  \
    _tree_event->Branch(                        \
        #b, &_branch_ ## b);

    MEMBER_BRANCH;

#undef BRANCH
#undef BRANCH_ARRAY
#undef BRANCH_ARRAY2
#undef BRANCH_STR
#undef BRANCH_STR_ARRAY

    PostData(1, _tree_event);

    /////////////////////////////////////////////////////////////////
}

#undef MEMBER_BRANCH

Bool_t AliAnalysisTaskNTGJ::Run()
{
    AliClusterContainer *cluster_container = NULL;
    std::vector<AliTrackContainer*> *track_containers = new std::vector<AliTrackContainer*>;
    AliMCParticleContainer *mc_container = NULL;

    loadEmcalGeometry(); // Fernando
    getEvent(); // Alwina
    getContainers(cluster_container, track_containers, mc_container);
    setTrackCuts(); // Dhruv
    getMultiplicityCentralityEventPlane(); // Alwina
    loadPhotonNNModel(); // Fernando
    loadMC(); // Preeti
    getBeamProperties(); // Alwina
    if(skimClusterE()) {  // Fernando
        return false;
    }
    getMetadata(); // Alwina
    getEmcalCellInfo(); // Fernando
    getPrimaryMCParticles(); // Rey
    initializeFastjetVectors();
    doTrackLoop(); // Dhruv
    doClusterLoopForAreaDetermination();
    computeVoronoiAreas();
    initializeMoreFastjetVectors();
    doMCParticleLoop(); // Rey
    fastjetTruthJets();
    doClusterLoopForJets();
    doManyFastjetThings();
    getUEEstimate();
    doClusterLoop(); // Fernando, except isolation stuff
    fillJetBranches();
    skimJets();
    fillCellBranches(); // Fernando
    fillMuonBranches();
    fillEgNtrial(); // Alwina

    _tree_event->Fill();
    return true;
}

void AliAnalysisTaskNTGJ::loadEmcalGeometry()
{

}

void AliAnalysisTaskNTGJ::getEvent()
{

}

void AliAnalysisTaskNTGJ::getContainers(AliClusterContainer *cluster_container,
                                        std::vector<AliTrackContainer*> *track_containers,
                                        AliMCParticleContainer *mc_container)
{
    cluster_container = GetClusterContainer(0);
    track_containers->push_back(GetTrackContainer(0));
    if (_is_embed) {
        track_containers->push_back(GetTrackContainer(1));
    }
    mc_container = GetMCParticleContainer("mcparticles");

    AliDebugStream(2) << "Event has " << cluster_container->GetNAcceptEntries() << " clusters" << std::endl;
    int ntrack = 0;
    for (auto track_container : *track_containers) {
        ntrack += track_container->GetNAcceptEntries();
    }
    AliDebugStream(2) << "Event has " << ntrack << " tracks" << std::endl;
    if (mc_container) {
        AliDebugStream(2) << "Event has " << mc_container->GetNParticles() << " MC particles" << std::endl;
    }
}

void AliAnalysisTaskNTGJ::setTrackCuts()
{

}

void AliAnalysisTaskNTGJ::getMultiplicityCentralityEventPlane()
{

}

void AliAnalysisTaskNTGJ::loadPhotonNNModel()
{

}

void AliAnalysisTaskNTGJ::loadMC()
{

}

void AliAnalysisTaskNTGJ::getBeamProperties()
{

}

bool AliAnalysisTaskNTGJ::skimClusterE()
{
    return false;
}

void AliAnalysisTaskNTGJ::getMetadata()
{

}

void AliAnalysisTaskNTGJ::getEmcalCellInfo()
{

}

void AliAnalysisTaskNTGJ::getPrimaryMCParticles()
{

}

void AliAnalysisTaskNTGJ::initializeFastjetVectors()
{

}

void AliAnalysisTaskNTGJ::doTrackLoop()
{

}

void AliAnalysisTaskNTGJ::doClusterLoopForAreaDetermination()
{

}

void AliAnalysisTaskNTGJ::computeVoronoiAreas()
{

}

void AliAnalysisTaskNTGJ::initializeMoreFastjetVectors()
{

}

void AliAnalysisTaskNTGJ::doMCParticleLoop()
{

}

void AliAnalysisTaskNTGJ::fastjetTruthJets()
{

}

void AliAnalysisTaskNTGJ::doClusterLoopForJets()
{

}

void AliAnalysisTaskNTGJ::doManyFastjetThings()
{

}

void AliAnalysisTaskNTGJ::getUEEstimate()
{

}

void AliAnalysisTaskNTGJ::doClusterLoop()
{
    fillClusterBranches();
    fillIsolationBranches();
    fillPhotonNNBranches();
}

void AliAnalysisTaskNTGJ::fillClusterBranches()
{

}

void AliAnalysisTaskNTGJ::fillIsolationBranches()
{

}

void AliAnalysisTaskNTGJ::fillPhotonNNBranches()
{

}

void AliAnalysisTaskNTGJ::fillJetBranches()
{

}

void AliAnalysisTaskNTGJ::skimJets()
{

}

void AliAnalysisTaskNTGJ::fillCellBranches()
{

}

void AliAnalysisTaskNTGJ::fillMuonBranches()
{

}

void AliAnalysisTaskNTGJ::fillEgNtrial()
{

}

AliEMCALRecoUtils *AliAnalysisTaskNTGJ::GetEMCALRecoUtils(void)
{
    return _reco_util;
}

void AliAnalysisTaskNTGJ::SetAliROOTVersion(const char *version)
{
    strncpy(_branch_version_aliroot, version, BUFSIZ);
    _branch_version_aliroot[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetAliPhysicsVersion(const char *version)
{
    strncpy(_branch_version_aliphysics, version, BUFSIZ);
    _branch_version_aliphysics[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetGridDataDir(const char *dir)
{
    strncpy(_branch_grid_data_dir, dir, BUFSIZ);
    _branch_grid_data_dir[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::SetGridDataPattern(const char *pattern)
{
    strncpy(_branch_grid_data_pattern, pattern, BUFSIZ);
    _branch_grid_data_pattern[BUFSIZ - 1] = '\0';
}

void AliAnalysisTaskNTGJ::
SetEMCALGeometryFilename(const char *emcal_geometry_filename)
{
    _emcal_geometry_filename = emcal_geometry_filename;
}

void AliAnalysisTaskNTGJ::
SetEMCALLocal2MasterFilename(const char *emcal_local2master_filename)
{
    _emcal_local2master_filename = emcal_local2master_filename;
}

void AliAnalysisTaskNTGJ::
SetForceUESubtraction(bool force_ue_subtraction)
{
    _force_ue_subtraction = force_ue_subtraction;
}

void AliAnalysisTaskNTGJ::SetSkimClusterMinE(double min_e)
{
    _skim_cluster_min_e = min_e;
}

void AliAnalysisTaskNTGJ::SetSkimTrackMinPt(double min_pt)
{
    _skim_track_min_pt = min_pt;
}

void AliAnalysisTaskNTGJ::SetSkimMuonTrackMinPt(double min_pt)
{
    _skim_muon_track_min_pt = min_pt;
}

void AliAnalysisTaskNTGJ::SetSkimJetMinPt(double min_pt_1,
                                          double min_pt_2,
                                          double min_pt_3)
{
    _skim_jet_min_pt.clear();

    if (min_pt_1 != -INFINITY) {
        _skim_jet_min_pt.push_back(min_pt_1);
        if (min_pt_2 != -INFINITY) {
            _skim_jet_min_pt.push_back(min_pt_2);
            if (min_pt_3 != -INFINITY) {
                _skim_jet_min_pt.push_back(min_pt_3);
            }
        }
    }
}

void AliAnalysisTaskNTGJ::SetSkimJetAveragePt(double average_pt)
{
    _skim_jet_average_pt = average_pt;
}

void AliAnalysisTaskNTGJ::SetSkimMultiplicityTrackletMinN(int min_n)
{
    _skim_multiplicity_tracklet_min_n = min_n;
}

void AliAnalysisTaskNTGJ::SetStoredTrackMinPt(double min_pt)
{
    _stored_track_min_pt = min_pt;
}

// This is primarily intended for JEC derivation, and therefore
// defined in pT raw

void AliAnalysisTaskNTGJ::SetStoredJetMinPtRaw(double min_pt_raw)
{
    _stored_jet_min_pt_raw = min_pt_raw;
}

void AliAnalysisTaskNTGJ::
SetNRandomIsolation(unsigned int nrandom_isolation)
{
    _nrandom_isolation = nrandom_isolation;
}

void AliAnalysisTaskNTGJ::SetIsEmbed(bool is_embed)
{
    _is_embed = is_embed;
}
