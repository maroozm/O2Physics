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

  Configurable<float> mCentralEta{"centralEta", 0.8, "eta limit for central tracks"};
  Configurable<float> mVertexZ{"vertexZ", 10, "z vertex limit"};
  Configurable<std::vector<float>> mPtbins{"ptCuts", {0.4f, 1.0f, 0.4f, 0.6f, 0.8f, 1.0f, 0.4f, 2.0f, 0.4f, 5.0f}, "pT cuts"};

  Filter etaTracks = nabs(aod::track::eta) < mCentralEta;
  // Filter centEvents = aod::cent::FV0A > 0 && aod::cent::FV0A < 100;

  array<int, 40> mbinsscale;
  std::vector<std::shared_ptr<TH2>> histArr;

  // Histograms
  HistogramRegistry histos{
    "histos",
    {
      {"collID", "collisionID", {HistType::kTH1I, {{1000, -10000, 10000}}}},
      {"centFT0C0A", "centFT0C_0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"centFV0A", "centFV0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"centFT0A", "centFT0A", {HistType::kTH1F, {{100, 0, 100}}}},
      {"vertexX", "vertexX", {HistType::kTH1F, {{100, -10, 10}}}},
      {"vertexY", "vertexY", {HistType::kTH1F, {{100, -10, 10}}}},
      {"vertexZ", "vertexZ", {HistType::kTH1F, {{100, -10, 10}}}},
      {"eta", "#eta", {HistType::kTH1F, {{100, -2, 2}}}},
      {"pT", "#pt", {HistType::kTH1F, {{100, -0.01, 20}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    true // enable booking
  };

  void init(o2::framework::InitContext&)
  {
    for (auto iM = 0; iM < 40; ++iM) {
      mbinsscale[iM] = 2 * (iM + 2);
    }
    for (auto iPt = 0; iPt < 5; ++iPt) {
      for (auto iM = 0; iM < 40; ++iM) {
        auto mScalers = std::get<std::shared_ptr<TH2>>(histos.add(Form("bin%i/hetaphi%i", iPt, iM), Form("#eta#phi_%iforbin%i", iM, iPt), HistType::kTH2F, {{mbinsscale[iM], -0.8, 0.8}, {mbinsscale[iM], 0, 2 * TMath::Pi()}}));
        histArr.push_back(mScalers);
      }
    }
  }

  template <class T>
  void checkpT(const T& track)
  {
    for (auto iPt = 0; iPt < 5; ++iPt) {
      if (track.pt() > mPtbins.value[2 * iPt] && track.pt() < mPtbins.value[2 * iPt + 1]) {
        for (auto iM = 0; iM < 40; ++iM) {
          histArr[iPt * 40 + iM]->Fill(track.eta(), track.phi());
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
    histos.fill(HIST("centFT0C0A"), coll.centFT0M());
    histos.fill(HIST("centFV0A"), coll.centFV0A());
    histos.fill(HIST("centFT0A"), coll.centFT0A());

    for (auto const& h : histArr) {
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