import FWCore.ParameterSet.Config as cms

l1tkmuMaker = cms.EDProducer(
    "L1TkMuMaker",
    aliasPrefix = cms.untracked.string("l1tkmus"),
    l1tkmuInputTag = cms.InputTag("L1TkMuonsMerge")
    )
