#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include <iostream>
#include <THashList.h>
#include <array>

using std::array;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(mEvMult, evmult, int);
DECLARE_SOA_COLUMN(mCollID, collid, int);
DECLARE_SOA_COLUMN(mVertexZ, vertexZ, float);
DECLARE_SOA_COLUMN(mCentrality, evcentrality, float);
DECLARE_SOA_COLUMN(mAbContent, avgbincontent, float);
DECLARE_SOA_COLUMN(mMsquare, msquare, float);
DECLARE_SOA_COLUMN(mF2, fq2, float);
DECLARE_SOA_COLUMN(mF3, fq3, float);
DECLARE_SOA_COLUMN(mF4, fq4, float);
DECLARE_SOA_COLUMN(mF5, fq5, float);
} // namespace full

DECLARE_SOA_TABLE(nfmFULL, "AOD", "nfmFULL", full::mEvMult, full::mCollID, full::mVertexZ, full::mCentrality, full::mAbContent, full::mMsquare, full::mF2, full::mF3, full::mF4, full::mF5);
} // namespace o2::aod

// Write information to the tree
struct AnalysisExec {
  Produces<o2::aod::nfmFULL> eventsel;
  OutputObj<THashList> fOutputList{"eventStandard"};

  Configurable<float> mcentralEta{"centralEta", 0.8, "eta limit for central tracks"};
  Configurable<float> mvertexZ{"vertexZ", 10, "z vertex limit"};
  Configurable<std::vector<float>> mpTbins{"ptCuts", {0.4f, 1.0f, 0.4f, 0.6f, 0.8f, 1.0f, 0.4f, 2.0f, 0.4f, 5.0f}, "pT cuts"};
  Configurable<int> mnumPt{"numPt", 5, "number of pT bins"};

  Filter etaTracks = nabs(aod::track::eta) < mcentralEta;
  // Filter centEvents = aod::cent::FV0A > 0 && aod::cent::FV0A < 100;

  // Histograms
  HistogramRegistry histos{
    "histos",
    {
      {"collID", "collisionID", {HistType::kTH1I, {{1000, -10000, 10000}}}},
      {"centFT0M", "centFT0C_0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"centFV0A", "centFV0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"centFT0A", "centFT0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"centFT0C", "centFT0C", {HistType::kTH1F, {{100, 0, 100}}}},
      {"vertexX", "vertexX", {HistType::kTH1F, {{100, -10, 10}}}},
      {"vertexY", "vertexY", {HistType::kTH1F, {{100, -10, 10}}}},
      {"vertexZ", "vertexZ", {HistType::kTH1F, {{100, -10, 10}}}},
      {"eta", "#eta", {HistType::kTH1F, {{1000, -2, 2}}}},
      {"pT", "#pt", {HistType::kTH1F, {{1000, -0.01, 50}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    true};

  array<int, 40> mbinsscale;
  std::vector<std::shared_ptr<TH2>> histArrR;
  std::vector<std::shared_ptr<TH2>> histArrUR;

  void init(o2::framework::InitContext&)
  {
    for (auto iM = 0; iM < 40; ++iM) {
      mbinsscale[iM] = 2 * (iM + 2);
    }
    for (auto iPt = 0; iPt < mnumPt; ++iPt) {
      histos.add(Form("bin%i/eta", iPt + 1), "#eta distribution;#eta;", HistType::kTH1F, {{1000, -2, 2}});
      histos.add(Form("bin%i/pT", iPt + 1), "#p_{T} distribution;#p_{T};", HistType::kTH1F, {{1000, -0.01, 50}});
      histos.add(Form("bin%i/phi", iPt + 1), "#phi distribution;#phi;", HistType::kTH1F, {{1000, 0, 2 * TMath::Pi()}});
      for (auto iM = 0; iM < 40; ++iM) {
        auto mHistsR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/Reset/hetaphi%i", iPt + 1, iM), Form("#eta#phi_%i for bin%i;#eta;#phi", iM, iPt + 1), HistType::kTH2F, {{mbinsscale[iM], -0.8, 0.8}, {mbinsscale[iM], 0, 2 * TMath::Pi()}}));
        histArrR.push_back(mHistsR);
        auto mHistsUR = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/UReset/hetaphi%i", iPt + 1, iM), Form("#eta#phi_%i for bin%i;#eta;#phi", iM, iPt + 1), HistType::kTH2F, {{mbinsscale[iM], -0.8, 0.8}, {mbinsscale[iM], 0, 2 * TMath::Pi()}}));
        histArrUR.push_back(mHistsUR);
      }
    }
  }

  template <class T>
  void checkpT(const T& track)
  {
    for (auto iPt = 0; iPt < mnumPt; ++iPt) {
      if (track.pt() > mpTbins.value[2 * iPt] && track.pt() < mpTbins.value[2 * iPt + 1]) {
        histos.fill(HIST(Form("bin%i/eta", iPt + 1)), track.eta());
        histos.fill(HIST(Form("bin%i/pT", iPt + 1)), track.pt());
        histos.fill(HIST(Form("bin%i/phi", iPt + 1)), track.phi());
        for (auto iM = 0; iM < 40; ++iM) {
          histArrR[iPt * 40 + iM]->Fill(track.eta(), track.phi());
          histArrUR[iPt * 40 + iM]->Fill(track.eta(), track.phi());
        }
      }
    }
  }

  using CollisionFMs = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  void process(CollisionFMs::iterator const& coll, soa::Filtered<aod::Tracks> const& tracks)
  {
    histos.fill(HIST("vertexX"), coll.posX());
    histos.fill(HIST("vertexY"), coll.posY());
    histos.fill(HIST("vertexZ"), coll.posZ());
    histos.fill(HIST("centFT0M"), coll.centFT0M());
    histos.fill(HIST("centFV0A"), coll.centFV0A());
    histos.fill(HIST("centFT0A"), coll.centFT0A());
    histos.fill(HIST("centFT0C"), coll.centFT0C());

    for (auto const& h : histArrR) {
      h->Reset();
    }
    for (auto const& track : tracks) {
      histos.fill(HIST("collID"), track.collisionId());
      histos.fill(HIST("eta"), track.eta());
      histos.fill(HIST("pT"), track.pt());
      checkpT(track);
    }

    // Calculate the normalized factorial moments
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisExec>(cfgc),
  };
}