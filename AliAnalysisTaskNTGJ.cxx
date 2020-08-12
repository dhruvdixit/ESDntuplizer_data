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
	AliVEvent *event = NULL;
	AliESDEvent *esd_event = NULL;
	AliAODEvent *aod_event = NULL;

	AliClusterContainer *cluster_container = NULL;
	std::vector<AliTrackContainer*> *track_containers = new std::vector<AliTrackContainer*>;
	AliMCParticleContainer *mc_container = NULL;

	loadEmcalGeometry(); // Fernando [testing]
	if (!getEvent(event, esd_event, aod_event)) { // getEvent returns false if the event is null
		return false;
	}
	getContainers(cluster_container, track_containers, mc_container);
	setTrackCuts(); // Dhruv
	getMultiplicityCentralityEventPlane(event);
	loadPhotonNNModel(); // Fernando [testing]

	if (mc_container) {
		loadMC(aod_event);
	}

	getBeamProperties(event, esd_event, aod_event);
	if (!skimMultiplicityTracklet(event)) { // skimMultiplicityTracklet returns true if we should keep the event
		return false;
	}
    if (!skimClusterE(cluster_container)) {  //Fernando: skimClusterE returns true if we should keep event [testing]
      return false; 
    }
	getMetadata(esd_event, aod_event);
	getEmcalCellInfo(); // Fernando

	std::vector<size_t> stored_mc_truth_index;
	std::vector<Int_t> reverse_stored_mc_truth_index;
	std::vector<Int_t> reverse_stored_parton_algorithmic_index;

	if (mc_container) {
		getPrimaryMCParticles(mc_container,
				              &stored_mc_truth_index,
                              &reverse_stored_mc_truth_index,
                              &reverse_stored_parton_algorithmic_index); // Rey
	}

	initializeFastjetVectors();
	doTrackLoop(); // Dhruv
	doClusterLoopForAreaDetermination();
	computeVoronoiAreas();
	initializeMoreFastjetVectors();

	if (mc_container) {
		doMCParticleLoop(mc_container,
				         esd_event,
				         reverse_stored_mc_truth_index); // Rey
	}

	fastjetTruthJets();
	doClusterLoopForJets();
	doManyFastjetThings();
	getUEEstimate();
    doClusterLoop(); // Fernando, except isolation stuff [lol]
    fillJetBranches();
    skimJets();
    fillCellBranches(); // Fernando [Testing]
    fillMuonBranches();
    fillEgNtrial();

    _tree_event->Fill();
    return true;
}

void AliAnalysisTaskNTGJ::loadEmcalGeometry()
{

  if (!_emcal_geometry_filename.empty() &&
      !_emcal_local2master_filename.empty() &&
      _emcal_geometry == NULL) {
    if (!gSystem->
	AccessPathName(gSystem->ExpandPathName(
		       _emcal_geometry_filename.c_str()))) 
      {
	TGeoManager::Import(_emcal_geometry_filename.c_str());
	_emcal_geometry = AliEMCALGeometry::
	  GetInstance(_emcal_geometry_name);
      }
    
    AliOADBContainer emcal_geometry_container("emcal");

    if (!gSystem->AccessPathName(gSystem->ExpandPathName(
	_emcal_local2master_filename.c_str()))) 
      {

	emcal_geometry_container.
	  InitFromFile(_emcal_local2master_filename.c_str(),"AliEMCALgeo");
      }

    TObjArray *geometry_matrix = dynamic_cast<TObjArray *>(
				 emcal_geometry_container.GetObject(
			         _branch_run_number, "EmcalMatrices"));

    if (_emcal_geometry != NULL && geometry_matrix != NULL) {
      const Int_t nsm = _emcal_geometry->GetEMCGeometry()->
	GetNumberOfSuperModules();

      for (Int_t sm = 0; sm < nsm; sm++) {

	_emcal_geometry->SetMisalMatrix(
			 dynamic_cast<TGeoHMatrix *>(
			 geometry_matrix->At(sm)),sm);
      }
      _branch_has_misalignment_matrix = true;
    }
  }

  if (_emcal_mask.size() != EMCAL_NCELL) {
    _emcal_mask.resize(EMCAL_NCELL);
#if 0 // Keep = 1 to for an actual EMCAL mask (and not all channels
      // turned on)
    for (unsigned int i = 0; i < EMCAL_NCELL; i++) {
      _emcal_mask[i] = inside_edge(i, 1);
    }
    // #include "bad_channel.h"
    for (unsigned int i = 0; bad_channel_emcal[i] != -1; i++) {
      if (inside_edge(bad_channel_emcal[i], 1)) {
	unsigned int bad_cell_3_3[9];

	cell_3_3(bad_cell_3_3, bad_channel_emcal[i]);
	for (size_t j = 0; j < 9; j++) {
	  _emcal_mask[bad_cell_3_3[j]] = false;
	}
      }
    }
#else
    for (unsigned int i = 0; i < EMCAL_NCELL; i++) {
      _emcal_mask[i] = true;
    }
#endif
  }
}

bool AliAnalysisTaskNTGJ::getEvent(AliVEvent *&event,
		AliESDEvent *&esd_event,
		AliAODEvent *&aod_event)
{
	// lines 406-432, 436-449
	event = InputEvent();

	if (event == NULL) {
		return false;
	}

	if (event->GetRunNumber() != _run_number_current) {
		_metadata_filled = false;
		_run_number_current = event->GetRunNumber();
		if (_muon_track_cut != NULL) {
			_muon_track_cut->SetAllowDefaultParams(kTRUE);
			_muon_track_cut->SetRun(
					(AliInputEventHandler *)
					((AliAnalysisManager::GetAnalysisManager())->
					 GetInputEventHandler()));
		}
	}

	esd_event = dynamic_cast<AliESDEvent *>(event);

	if (esd_event == NULL) {
		aod_event = dynamic_cast<AliAODEvent *>(event);
		if (aod_event == NULL) {
			return false;
		}
	}

	_branch_run_number = event->GetRunNumber();
	_branch_period_number = event->GetPeriodNumber();
	_branch_orbit_number = event->GetOrbitNumber();
	_branch_bunch_crossing_number = event->GetBunchCrossNumber();
	_branch_time_stamp = esd_event != NULL ?
		esd_event->GetTimeStamp() :
		aod_event->GetTimeStamp();
	_branch_trigger_mask[0] = esd_event != NULL ?
		esd_event->GetTriggerMask() :
		aod_event->GetTriggerMask();
	_branch_trigger_mask[1] = esd_event != NULL ?
		esd_event->GetTriggerMaskNext50() :
		aod_event->GetTriggerMaskNext50();
	_branch_has_misalignment_matrix = false;

	return true;
}

void AliAnalysisTaskNTGJ::getContainers(AliClusterContainer *&cluster_container,
		std::vector<AliTrackContainer*> *&track_containers,
		AliMCParticleContainer *&mc_container)
{
	// set cluster_container, track_containers, and mc_container
	// by taking in references to pointers
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

void AliAnalysisTaskNTGJ::getMultiplicityCentralityEventPlane(AliVEvent *event)
{
	// lines 547-604
	AliVVZERO *v0 = event->GetVZEROData();

	if (v0 != NULL) {
		for (size_t i = 0; i < 64; i++) {
			_branch_multiplicity_v0[i] = v0->GetMultiplicity(i);
		}
	}

	static const char *centrality_method[9] = {
		"V0M", "CL0", "CL1",
		"V0Mplus05", "V0Mplus10", "V0Mminus05", "V0Mminus10",
		"SPDClustersCorr", "SPDTracklets"
	};

	AliMultSelection *mult_selection = static_cast<AliMultSelection *>
		(event->FindListObject("MultSelection"));

	std::fill(_branch_centrality,
			_branch_centrality + sizeof(_branch_centrality) /
			sizeof(*_branch_centrality), NAN);
	if (mult_selection != NULL) {
		for (size_t i = 0; i < 9; i++) {
			_branch_centrality[i] =
				mult_selection->GetMultiplicityPercentile(
						centrality_method[i]);
		}
	}
	else {
		AliCentrality *centrality = event->GetCentrality();

		if (centrality != NULL) {
			for (size_t i = 0; i < 9; i++) {
				_branch_centrality[i] =
					centrality->GetCentralityPercentile(
							centrality_method[i]);
			}
		}
	}
	// Copy for easy cutting (where duplicated global variable are not
	// very wasteful)
	_branch_centrality_v0m = _branch_centrality[0];

	std::fill(_branch_event_plane_psi_v0,
			_branch_event_plane_psi_v0 +
			sizeof(_branch_event_plane_psi_v0) /
			sizeof(*_branch_event_plane_psi_v0), NAN);

	AliEventplane *event_plane = event->GetEventplane();

	if (event_plane != NULL) {
		for (size_t i = 1; i < 4; i++) {
			_branch_event_plane_psi_v0[i - 1] =
				event->GetEventplane()->CalculateVZEROEventPlane(
						event, 10, i,
						_branch_event_plane_q_v0[i - 1][0],
						_branch_event_plane_q_v0[i - 1][1]);
		}
	}
}

void AliAnalysisTaskNTGJ::loadPhotonNNModel()
{
  if (_keras_model_photon_discrimination == NULL) {
    _keras_model_photon_discrimination = new KerasModel;
    reinterpret_cast<KerasModel *>(
      _keras_model_photon_discrimination)->
      LoadModel("photon_discr.model");
  }
}

void AliAnalysisTaskNTGJ::loadMC(AliAODEvent *aod_event)
{
	// lines 616-683
	_branch_eg_signal_process_id = INT_MIN;
	_branch_eg_mpi = INT_MIN;
	_branch_eg_pt_hat = NAN;
	_branch_eg_cross_section = NAN;
	_branch_eg_weight = NAN;
	// This should be default to zero to avoid counting mishap
	_branch_eg_ntrial = 0;

	_branch_eg_scale_pdf = NAN;
	_branch_eg_alpha_qcd = NAN;
	_branch_eg_alpha_qed = NAN;
	std::fill(_branch_eg_pdf_id,
			_branch_eg_pdf_id + sizeof(_branch_eg_pdf_id) /
			sizeof(*_branch_eg_pdf_id), INT_MIN);
	std::fill(_branch_eg_pdf_x,
			_branch_eg_pdf_x + sizeof(_branch_eg_pdf_x) /
			sizeof(*_branch_eg_pdf_x), NAN);
	std::fill(_branch_eg_pdf_x_pdf,
			_branch_eg_pdf_x_pdf + sizeof(_branch_eg_pdf_x_pdf) /
			sizeof(*_branch_eg_pdf_x_pdf), NAN);

	AliGenPythiaEventHeader *mc_truth_pythia_header = NULL;

	// embedding MC header is done separately for now
	if (_is_embed) {
		mc_truth_pythia_header = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetPythiaHeader();
	}
	else {
		AliAODMCHeader *aod_mc_header = dynamic_cast<AliAODMCHeader*>
			(aod_event->GetList()->FindObject(AliAODMCHeader::StdBranchName()));

		Int_t nGenerators = aod_mc_header->GetNCocktailHeaders();
		for (Int_t igen = 0; igen < nGenerators; igen++) {
			AliGenEventHeader *eventHeaderGen = aod_mc_header->GetCocktailHeader(igen);
			TString name = eventHeaderGen->GetName();

			if (name.CompareTo("Pythia") == 0) {
				mc_truth_pythia_header = dynamic_cast<AliGenPythiaEventHeader*>(eventHeaderGen);
			}
		}
	}

	if (mc_truth_pythia_header != NULL) {
		_branch_eg_weight = mc_truth_pythia_header->EventWeight();

		TArrayF eg_primary_vertex(3);

		mc_truth_pythia_header->PrimaryVertex(eg_primary_vertex);
		for (Int_t i = 0; i < 3; i++) {
			_branch_eg_primary_vertex[i] = eg_primary_vertex.At(i);
		}

		_branch_eg_signal_process_id = mc_truth_pythia_header->ProcessType();
		_branch_eg_mpi = mc_truth_pythia_header->GetNMPI();
		_branch_eg_pt_hat = mc_truth_pythia_header->GetPtHard();

		_branch_eg_cross_section =  mc_truth_pythia_header->GetXsection();

		// Count ntrial, because the event might get skimmed away
		// by the ntuplizer
		_skim_sum_eg_ntrial += mc_truth_pythia_header->Trials();
	}
}

void AliAnalysisTaskNTGJ::getBeamProperties(AliVEvent *event,
		AliESDEvent *esd_event,
		AliAODEvent *aod_event)
{
	// lines 685-792
	if (esd_event != NULL) {
		esd_event->InitMagneticField();
	}
	else if (aod_event != NULL) {
		aod_event->InitMagneticField();
	}

	getPrimaryVertex(event, esd_event);
	getPrimaryVertexSPD(event, esd_event, aod_event);
	getPileup(esd_event, aod_event);

	_branch_event_selected = fInputHandler->IsEventSelected();
}

void AliAnalysisTaskNTGJ::getPrimaryVertex(AliVEvent *event,
		AliESDEvent *esd_event)
{
	// lines 692-717
	std::fill(_branch_primary_vertex,
			_branch_primary_vertex +
			sizeof(_branch_primary_vertex) /
			sizeof(*_branch_primary_vertex), NAN);
	std::fill(_branch_primary_vertex_sigma,
			_branch_primary_vertex_sigma +
			sizeof(_branch_primary_vertex_sigma) /
			sizeof(*_branch_primary_vertex_sigma), NAN);
	_branch_primary_vertex_ncontributor = INT_MIN;

	const AliVVertex *primary_vertex = event->GetPrimaryVertex();

	if (primary_vertex != NULL) {
		primary_vertex->GetXYZ(_branch_primary_vertex);
	}
	if (esd_event != NULL) {
		const AliESDVertex *esd_primary_vertex =
			esd_event->GetPrimaryVertex();

		if (esd_primary_vertex != NULL) {
			esd_primary_vertex->GetSigmaXYZ(
					_branch_primary_vertex_sigma);
			_branch_primary_vertex_ncontributor =
				esd_primary_vertex->GetNContributors();
		}
	}
}

void AliAnalysisTaskNTGJ::getPrimaryVertexSPD(AliVEvent *event,
		AliESDEvent *esd_event,
		AliAODEvent *aod_event)
{
	// lines 719-745
	std::fill(_branch_primary_vertex_spd,
			_branch_primary_vertex_spd +
			sizeof(_branch_primary_vertex_spd) /
			sizeof(*_branch_primary_vertex_spd), NAN);
	std::fill(_branch_primary_vertex_spd_sigma,
			_branch_primary_vertex_spd_sigma +
			sizeof(_branch_primary_vertex_spd_sigma) /
			sizeof(*_branch_primary_vertex_spd_sigma), NAN);
	_branch_primary_vertex_spd_ncontributor = INT_MIN;

	const AliVVertex *primary_vertex_spd =
		event->GetPrimaryVertexSPD();

	if (primary_vertex_spd != NULL) {
		primary_vertex_spd->GetXYZ(_branch_primary_vertex_spd);
	}
	if (esd_event != NULL) {
		const AliESDVertex *esd_primary_vertex_spd =
			esd_event->GetPrimaryVertexSPD();

		if (esd_primary_vertex_spd != NULL) {
			esd_primary_vertex_spd->GetSigmaXYZ(
					_branch_primary_vertex_spd_sigma);
			_branch_primary_vertex_spd_ncontributor =
				esd_primary_vertex_spd->GetNContributors();
		}
	}
	else if (aod_event != NULL) {
		// lines 762-770
		const AliAODVertex *aod_primary_vertex_spd =
			aod_event->GetPrimaryVertexSPD();

		/*if (aod_primary_vertex_spd != NULL) {
		  aod_primary_vertex_spd->GetSigmaXYZ(
		  _branch_primary_vertex_spd_sigma);
		  _branch_primary_vertex_spd_ncontributor =
		  aod_primary_vertex_spd->GetNContributors();
		  }//*/
	}
}

void AliAnalysisTaskNTGJ::getPileup(AliESDEvent *esd_event,
		AliAODEvent *aod_event)
{
	// lines 747-761, 772-780
	_branch_is_pileup_from_spd_3_08 = false;
	_branch_is_pileup_from_spd_5_08 = false;
	_branch_npileup_vertex_spd = INT_MIN;
	_branch_ncluster_tpc = INT_MIN;
	if (esd_event != NULL) {
		_branch_is_pileup_from_spd_3_08 =
			esd_event->IsPileupFromSPD(3, 0.8);
		_branch_is_pileup_from_spd_5_08 =
			esd_event->IsPileupFromSPD(5, 0.8);
		_branch_npileup_vertex_spd =
			esd_event->GetNumberOfPileupVerticesSPD();
		_branch_ncluster_tpc =
			esd_event->GetNumberOfTPCClusters();
	}
	else if (aod_event != NULL) {
		_branch_is_pileup_from_spd_3_08 =
			aod_event->IsPileupFromSPD(3, 0.8);
		_branch_is_pileup_from_spd_5_08 =
			aod_event->IsPileupFromSPD(5, 0.8);
		_branch_npileup_vertex_spd =
			aod_event->GetNumberOfPileupVerticesSPD();
		_branch_ncluster_tpc =
			aod_event->GetNumberOfTPCClusters();
	}
}

bool AliAnalysisTaskNTGJ::skimMultiplicityTracklet(AliVEvent *event)
{
	// lines 784-792
	// returns true if we should KEEP the event
	if (_skim_multiplicity_tracklet_min_n > INT_MIN) {
		AliVMultiplicity *multiplicity = event->GetMultiplicity();

		if (multiplicity == NULL ||
				!(multiplicity->GetNumberOfTracklets() >=
					_skim_multiplicity_tracklet_min_n)) {
			return false;
		}
	}

	return true;
}

bool AliAnalysisTaskNTGJ::skimClusterE(AliClusterContainer *calo_cluster)
{
  if (_skim_cluster_min_e > -INFINITY) {
    double cluster_e_max = -INFINITY;

    for (Int_t i = 0; i < calo_cluster->GetNEntries(); i++) {
            AliVCluster *c =
	      static_cast<AliVCluster *>(calo_cluster->GetCluster(i));

            if (!(c->GetNCells() > 1) ||
                cell_masked(c, _emcal_mask)) {
	      continue;
            }

            TLorentzVector p;

            c->GetMomentum(p, _branch_primary_vertex);
            cluster_e_max = std::max(cluster_e_max, p.E());
    }
    if (!(cluster_e_max >= _skim_cluster_min_e)) {
      // Discard this event
      return false;
    }
  }
  return true; //if no min_e set, return true to keep event
}

void AliAnalysisTaskNTGJ::getMetadata(AliESDEvent *esd_event,
		AliAODEvent *aod_event)
{
	// line 434, moved here because it doesn't seem to be used anywhere else
	alice_jec_t jec;

	// lines 831-889
	if (!_metadata_filled) {
		// Use gitattributes ident mechanism to track the blob object
		// name
		strncpy(_branch_id_git, "$Id$", BUFSIZ);
		_branch_id_git[BUFSIZ - 1] = '\0';
		strncpy(_branch_version_jec, jec.version(), BUFSIZ);
		_branch_version_jec[BUFSIZ - 1] = '\0';
		if (esd_event != NULL) {
			_branch_beam_energy = esd_event->GetBeamEnergy();
			for (size_t i = 0; i < 2; i++) {
				_branch_beam_particle[i] =
					esd_event->GetBeamParticle(i);
			}

			const AliESDRun *esd_run = esd_event->GetESDRun();

			if (esd_run != NULL) {
				_branch_trigger_class.clear();
				for (size_t i = 0; i < NTRIGGER_CLASS_MAX; i++) {
					_branch_trigger_class.push_back(
							std::string(esd_run->GetTriggerClass(i)));
				}
			}
		}
		else if (aod_event != NULL) {
			// FIXME: AOD not handled
			_branch_beam_energy = NAN;
			std::fill(_branch_beam_particle,
					_branch_beam_particle +
					sizeof(_branch_beam_particle) /
					sizeof(*_branch_beam_particle), 0);
		}

#ifdef WITH_EFP7
		strncpy(_branch_debug_blas_version,
				cblas.version_str().c_str(), BUFSIZ);
#else // WITH_EFP7
		_branch_debug_blas_version[0] = '\0';
#endif // WITH_EFP7

		_metadata_filled = true;
	}
	else {
		// To make sure no space is wasted, set metadata in all
		// subsequent events to empty
		_branch_id_git[0] = '\0';
		_branch_version_aliroot[0] = '\0';
		_branch_version_aliphysics[0] = '\0';
		_branch_version_jec[0] = '\0';
		_branch_grid_data_dir[0] = '\0';
		_branch_grid_data_pattern[0] = '\0';
		std::fill(_branch_beam_particle,
				_branch_beam_particle +
				sizeof(_branch_beam_particle) /
				sizeof(*_branch_beam_particle), 0);
		_branch_trigger_class.clear();

		_branch_debug_blas_version[0] = '\0';
	}
}

void AliAnalysisTaskNTGJ::getEmcalCellInfo()
{

}
// =================================================================================================================
void AliAnalysisTaskNTGJ::getPrimaryMCParticles(AliMCParticleContainer *mc_container,
                                                std::vector<size_t> *stored_mc_truth_index,
                                                std::vector<Int_t> *reverse_stored_mc_truth_index,
                                                std::vector<Int_t> *reverse_stored_parton_algorithmic_index)
{
	// lines 948-998
    AliDebugStream(3) << "loops 1 and 2 through MC container" << std::endl;

	stored_mc_truth_index->resize(mc_container->GetNParticles(), ULONG_MAX);
	size_t nmc_truth = 0;
	for (Int_t i = 0; i < mc_container->GetNParticles(); i++) {
		// Bookkeeping for primary final state particles
		if (final_state_primary(mc_container, i)) {
			stored_mc_truth_index->at(i) = nmc_truth;
			reverse_stored_mc_truth_index->push_back(i);
			nmc_truth++;
		}


		// Bookkeeping for partons
		if (parton_cms_algorithmic(mc_container, i)) {
			reverse_stored_parton_algorithmic_index->push_back(i);
		}
	}

	// Assign secondaries to the primaries
	for (Int_t i = 0;
			i < mc_container->GetNParticles(); i++) {
		// Skip primaries
		if (final_state_primary(mc_container, i)) {
			continue;
		}

		Int_t j = i;
		bool has_physical_primary_ancestor = false;

		// Ensure termination even if there is a loop
		for (size_t k = 0; k < 1000; k++) {
			const AliAODMCParticle *p = mc_container->GetMCParticle(j);
			if (p == NULL) {
				break;
			}
			j = p->GetMother();
			if (!(j >= 0 &&
						j < mc_container->GetNParticles())) {
				break;
			}
			if (final_state_primary(mc_container, j)) {
				has_physical_primary_ancestor = true;
				break;
			}
		}
		if (has_physical_primary_ancestor) {
			stored_mc_truth_index->at(i) = stored_mc_truth_index->at(j);
		}
	}	
}
// ================================================================================================================= 
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
// ================================================================================================================= 
void AliAnalysisTaskNTGJ::doMCParticleLoop(AliMCParticleContainer *mc_container,
                                           AliESDEvent *esd_event,
                                           std::vector<Int_t> reverse_stored_mc_truth_index)
{
    // lines 1336-1376, 1406-1474
    AliDebugStream(3) << "loop 3 through MC container" << std::endl;

    enum {
        BEAM_PARTICLE_P = 1001
    };

    const bool subtract_ue =
        _force_ue_subtraction ? true :
        esd_event != NULL ?
        !(esd_event->GetBeamParticle(0) == BEAM_PARTICLE_P &&
                esd_event->GetBeamParticle(1) == BEAM_PARTICLE_P) :
        false;

    // may need to move this outside of this function
    _branch_nmc_truth = 0;

	enum {
		PDG_CODE_ELECTRON_MINUS             = 11,
		PDG_CODE_ELECTRON_NEUTRINO          = 12,
		PDG_CODE_MUON_MINUS                 = 13,
		PDG_CODE_MUON_NEUTRINO              = 14,
		PDG_CODE_TAU_MINUS                  = 15,
		PDG_CODE_TAU_NEUTRINO               = 16,
		PDG_CODE_PHOTON                     = 22,
	};

	double met_truth_kahan_error[2] = { 0, 0 };

	for (std::vector<Int_t>::const_iterator iterator =
			reverse_stored_mc_truth_index.begin();
			iterator != reverse_stored_mc_truth_index.end();
			iterator++) {
		const AliAODMCParticle *p = mc_container->GetMCParticle(*iterator);

		if (p == NULL) {
			// Keep consistent indexing, though this should never
			// happen
			_branch_mc_truth_e[_branch_nmc_truth] = NAN;
			_branch_mc_truth_pt[_branch_nmc_truth] = NAN;
			_branch_mc_truth_eta[_branch_nmc_truth] = NAN;
			_branch_mc_truth_phi[_branch_nmc_truth] = NAN;
			_branch_mc_truth_charge[_branch_nmc_truth] =
				CHAR_MIN;
			_branch_mc_truth_pdg_code[_branch_nmc_truth] =
				SHRT_MIN;
			_branch_mc_truth_status[_branch_nmc_truth] =
				UCHAR_MAX;
			_branch_mc_truth_generator_index[_branch_nmc_truth] =
				UCHAR_MAX;
			_branch_nmc_truth++;
			continue;
		}

		if (p->GetGeneratorIndex() != 0 || !subtract_ue) {
			const unsigned int abs_pdg_code =
				std::abs(p->PdgCode());

			switch (abs_pdg_code) {
				case PDG_CODE_ELECTRON_NEUTRINO:
				case PDG_CODE_MUON_NEUTRINO:
				case PDG_CODE_TAU_NEUTRINO:
					// Remove all (stable) neutrinos from the truth
					// jet definition
					break;
				default:
					kahan_sum(_branch_met_truth[0],
							met_truth_kahan_error[0], p->Px());
					kahan_sum(_branch_met_truth[1],
							met_truth_kahan_error[1], p->Py());
			}
		}

		_branch_mc_truth_e[_branch_nmc_truth] = half(p->E());
		_branch_mc_truth_pt[_branch_nmc_truth] = half(p->Pt());
		_branch_mc_truth_eta[_branch_nmc_truth] = half(p->Eta());
		_branch_mc_truth_phi[_branch_nmc_truth] =
			half(angular_range_reduce(p->Phi()));
		_branch_mc_truth_charge[_branch_nmc_truth] =
			std::min(static_cast<Short_t>(CHAR_MAX),
					std::max(static_cast<Short_t>(CHAR_MIN),
						p->Charge()));
		_branch_mc_truth_pdg_code[_branch_nmc_truth] =
			std::min(SHRT_MAX,
					std::max(SHRT_MIN, p->PdgCode()));
		_branch_mc_truth_status[_branch_nmc_truth] =
			std::min(static_cast<UInt_t>(UCHAR_MAX),
					std::max(static_cast<UInt_t>(0),
						mc_container->
						GetMCParticle(*iterator)->
						MCStatusCode()));
		// _branch_mc_truth_final_state_primary[_branch_nmc_truth] =
		//     final_state_primary(mc_container, *iterator);
		// _branch_mc_truth_first_parent[_branch_nmc_truth] =
		//     p->GetMother();
		// _branch_mc_truth_first_child[_branch_nmc_truth] =
		//     p->GetDaughterFirst();
		// _branch_mc_truth_first_child[_branch_nmc_truth] =
		//     p->GetDaughterLast();
		_branch_mc_truth_generator_index[_branch_nmc_truth] =
			std::min(static_cast<Short_t>(UCHAR_MAX),
					std::max(static_cast<Short_t>(0),
						p->GetGeneratorIndex()));

		if (p->GetMother() >= 0 &&
				p->GetMother() <
				mc_container->GetNParticles()) {
			const AliAODMCParticle *pp =
				mc_container->GetMCParticle(
						p->GetMother());

			_branch_mc_truth_first_parent_pdg_code
				[_branch_nmc_truth] = pp->PdgCode();
			_branch_mc_truth_first_parent_e
				[_branch_nmc_truth] = half(pp->E());
			_branch_mc_truth_first_parent_pt
				[_branch_nmc_truth] = half(pp->Pt());
			_branch_mc_truth_first_parent_eta
				[_branch_nmc_truth] = half(pp->Eta());
			_branch_mc_truth_first_parent_phi
				[_branch_nmc_truth] = half(pp->Phi());
			_branch_mc_truth_sibling_index[_branch_nmc_truth] =
				*iterator ==
				pp->GetDaughterFirst() ?
				pp->GetDaughterLast() :
				*iterator ==
				pp->GetDaughterLast() ?
				pp->GetDaughterFirst() :
				USHRT_MAX;
		}

		_branch_nmc_truth++;
		if (_branch_nmc_truth >= NMC_TRUTH_MAX) {
			break;
		}
	} 

}
// ================================================================================================================= 
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
  std::fill(&_branch_cell_position[0][0],
              &_branch_cell_position[0][0] +
	    sizeof(_branch_cell_position) /
	    sizeof(_branch_cell_position[0][0]), NAN);
  std::fill(_branch_cell_voronoi_area,
              _branch_cell_voronoi_area +
	    sizeof(_branch_cell_voronoi_area) /
	    sizeof(*_branch_cell_voronoi_area), NAN);
  if (_emcal_geometry != NULL) {
    _emcal_cell_position = new std::vector<point_2d_t>();
    for (unsigned int cell_id = 0; cell_id < EMCAL_NCELL;
	 cell_id++) {
      TVector3 v;

      _emcal_geometry->GetGlobal(cell_id, v);
      v -= TVector3(_branch_primary_vertex);

      _branch_cell_position[cell_id][0] = v.X();
      _branch_cell_position[cell_id][1] = v.Y();
      _branch_cell_position[cell_id][2] = v.Z();
      reinterpret_cast<std::vector<point_2d_t> *>
	(_emcal_cell_position)->push_back(
	    point_2d_t(v.Eta(), v.Phi()));
    }
    //FIXME: check to see if "computeVoronoiAreas" is just for jet area (is likeley)
    std::vector<double> emcal_cell_area;
    std::vector<std::set<size_t> > emcal_cell_incident;

    voronoi_area_incident(emcal_cell_area, emcal_cell_incident,
			  *reinterpret_cast<std::vector<point_2d_t> *>
			  (_emcal_cell_position));

    double sum_area_inside = 0;
    size_t count_inside = 0;

    for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
      if (inside_edge(cell_id, 1)) {
	sum_area_inside += emcal_cell_area[cell_id];
	count_inside++;
      }
    }

        const double mean_area_inside =
	  sum_area_inside / count_inside;
	
        for (int cell_id = 0; cell_id < EMCAL_NCELL; cell_id++) {
            _branch_cell_voronoi_area[cell_id] =
	      inside_edge(cell_id, 1) ?
	      emcal_cell_area[cell_id] : mean_area_inside;
        }
  }

}

void AliAnalysisTaskNTGJ::fillMuonBranches()
{

}

void AliAnalysisTaskNTGJ::fillEgNtrial()
{
	// lines 2512-2515
	// Now that the event is accepted, copy over the total counted
	// ntrials so far, and reset the ntrial counter
	_branch_eg_ntrial = _skim_sum_eg_ntrial;
	_skim_sum_eg_ntrial = 0;
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
