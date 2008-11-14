import FWCore.ParameterSet.Config as cms

process = cms.Process("DigiToRaw1")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource", fileNames =  cms.untracked.vstring('file:digis.root'))

process.load("Geometry.TrackerSimData.trackerSimGeometryXML_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")


#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("CalibTracker.Configuration.SiPixel_FakeConditions_cff")

process.load("CalibTracker.Configuration.SiPixelCabling.SiPixelCabling_SQLite_cff")
process.siPixelCabling.connect = 'sqlite_file:cabling.db'
process.siPixelCabling.toGet = cms.VPSet(cms.PSet(
    record = cms.string('SiPixelFedCablingMapRcd'),
    tag = cms.string('SiPixelFedCablingMap_v14')
))

process.load("EventFilter.SiPixelRawToDigi.SiPixelDigiToRaw_cfi")
#process.siPixelRawData.InputLabel = 'siPixelDigis'

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siPixelRawData'),
    destinations = cms.untracked.vstring('d2r'),
    d2r = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG'))
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName =  cms.untracked.string('file:rawdata.root'),
    outputCommands = cms.untracked.vstring("drop *","keep *_siPixelRawData_*_*")
)


process.p = cms.Path(process.siPixelRawData)
process.ep = cms.EndPath(process.out)



