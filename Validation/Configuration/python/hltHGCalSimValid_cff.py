import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import lcAssocByEnergyScoreProducer as _lcAssocByEnergyScoreProducer
from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import scAssocByEnergyScoreProducer as _scAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation as _layerClusterSimClusterAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation as _layerClusterCaloParticleAssociationProducer

from SimCalorimetry.HGCalAssociatorProducers.SimClusterToCaloParticleAssociation_cfi import SimClusterToCaloParticleAssociation
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import  allTrackstersToSimTrackstersAssociationsByLCs as _allTrackstersToSimTrackstersAssociationsByLCs
from SimCalorimetry.HGCalAssociatorProducers.hitToSimClusterCaloParticleAssociator_cfi import hitToSimClusterCaloParticleAssociator as _hitToSimClusterCaloParticleAssociator

from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabels as _hltTiclIterLabels

from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cff import recHitMapProducer as _recHitMapProducer

from Validation.Configuration.hltBarrelSimValid_cff import hltBarrelRecHitMapProducer as _hltBarrelRecHitMapProducer
from Validation.Configuration.hltBarrelSimValid_cff import barrel_hits
hgcal_hits = ["hltHGCalRecHit:HGCEERecHits", "hltHGCalRecHit:HGCHEFRecHits", "hltHGCalRecHit:HGCHEBRecHits"]
hltRecHitMapProducer = _hltBarrelRecHitMapProducer.clone()

hltHGCalRecHitMapProducer = _hltBarrelRecHitMapProducer.clone(
    hits = hgcal_hits,
    hgcalOnly = True,
)
from Configuration.Eras.Modifier_phase2_common_cff import phase2_common
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel
(phase2_common & ~ticl_barrel).toReplaceWith(hltRecHitMapProducer, hltHGCalRecHitMapProducer)

(phase2_common & ticl_barrel).toModify(hltRecHitMapProducer,
                                       hits = [*hgcal_hits, *barrel_hits],
                                       )

hltLcAssocByEnergyScoreProducer = _lcAssocByEnergyScoreProducer.clone(
    hits = 'hltHGCalRecHitMapProducer:RefProdVectorHGCRecHitCollection',
    hitMapTag = 'hltHGCalRecHitMapProducer:hgcalRecHitMap'
)

hltScAssocByEnergyScoreProducer = _scAssocByEnergyScoreProducer.clone(
    hits = 'hltHGCalRecHitMapProducer:RefProdVectorHGCRecHitCollection',
    hitMapTag = 'hltHGCalRecHitMapProducer:hgcalRecHitMap',
)

hltLayerClusterCaloParticleAssociationProducer = _layerClusterCaloParticleAssociationProducer.clone(
    associator = 'hltLcAssocByEnergyScoreProducer',
    label_lc = 'hltMergeLayerClusters'
)

hltLayerClusterSimClusterAssociationProducer = _layerClusterSimClusterAssociationProducer.clone(
    associator = 'hltScAssocByEnergyScoreProducer',
    label_lcl = 'hltMergeLayerClusters'
)

from SimCalorimetry.HGCalAssociatorProducers.AllLayerClusterToTracksterAssociatorsProducer_cfi import AllLayerClusterToTracksterAssociatorsProducer as _AllLayerClusterToTracksterAssociatorsProducer

hltAllLayerClusterToTracksterAssociations = _AllLayerClusterToTracksterAssociatorsProducer.clone(
    layer_clusters = 'hltMergeLayerClusters',
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels],
        cms.InputTag("hltTiclSimTracksters"),
        cms.InputTag("hltTiclSimTracksters", "fromCPs"),
    )
)

hltAllTrackstersToSimTrackstersAssociationsByLCs = _allTrackstersToSimTrackstersAssociationsByLCs.clone(
    allLCtoTSAccoc = 'hltAllLayerClusterToTracksterAssociations',
    layerClusters = 'hltMergeLayerClusters',
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels]
    ),
    simTracksterCollections = cms.VInputTag(
      cms.InputTag('hltTiclSimTracksters'),
      cms.InputTag('hltTiclSimTracksters','fromCPs')
    ),
)

from SimCalorimetry.HGCalAssociatorProducers.AllTracksterToSimTracksterAssociatorsByHitsProducer_cfi import AllTracksterToSimTracksterAssociatorsByHitsProducer as _AllTracksterToSimTracksterAssociatorsByHitsProducer

hltHitToSimClusterCaloParticleAssociator = _hitToSimClusterCaloParticleAssociator.clone(
    hitMap = 'hltHGCalRecHitMapProducer:hgcalRecHitMap',
    hits = 'hltHGCalRecHitMapProducer:RefProdVectorHGCRecHitCollection'
)

from SimCalorimetry.HGCalAssociatorProducers.AllHitToTracksterAssociatorsProducer_cfi import AllHitToTracksterAssociatorsProducer as _AllHitToTracksterAssociatorsProducer

hltAllHitToTracksterAssociations =  _AllHitToTracksterAssociatorsProducer.clone(
    hitMapTag = 'hltHGCalRecHitMapProducer:hgcalRecHitMap',
    hits = 'hltHGCalRecHitMapProducer:RefProdVectorHGCRecHitCollection',
    layerClusters = 'hltMergeLayerClusters',
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels],
        cms.InputTag("hltTiclSimTracksters"),
        cms.InputTag("hltTiclSimTracksters", "fromCPs"),
    )
)

hltAllTrackstersToSimTrackstersAssociationsByHits = _AllTracksterToSimTracksterAssociatorsByHitsProducer.clone(
    allHitToTSAccoc = 'hltAllHitToTracksterAssociations',
    hitToCaloParticleMap = 'hltHitToSimClusterCaloParticleAssociator:hitToCaloParticleMap',
    hitToSimClusterMap = 'hltHitToSimClusterCaloParticleAssociator:hitToSimClusterMap',
    hits = 'hltHGCalRecHitMapProducer:RefProdVectorHGCRecHitCollection',
    tracksterCollections = cms.VInputTag(
        *[cms.InputTag(label) for label in _hltTiclIterLabels]
    ),
    simTracksterCollections = cms.VInputTag(
      'hltTiclSimTracksters',
      'hltTiclSimTracksters:fromCPs'
    ),
)

hltHgcalAssociatorsTask = cms.Task(hltHGCalRecHitMapProducer,
                                   hltLcAssocByEnergyScoreProducer,
                                   hltScAssocByEnergyScoreProducer,
                                   SimClusterToCaloParticleAssociation,
                                   hltLayerClusterCaloParticleAssociationProducer,
                                   hltLayerClusterSimClusterAssociationProducer,
                                   hltAllLayerClusterToTracksterAssociations,
                                   hltAllTrackstersToSimTrackstersAssociationsByLCs,
                                   hltAllHitToTracksterAssociations,
                                   hltHitToSimClusterCaloParticleAssociator,
                                   hltAllTrackstersToSimTrackstersAssociationsByHits
                                   )

from SimCalorimetry.HGCalAssociatorProducers.hltLCToCPAssociation_cfi import (hltHGCalLCToCPAssociatorByEnergyScoreProducer,
                                                                              hltHGCalLayerClusterCaloParticleAssociation)
from SimCalorimetry.HGCalAssociatorProducers.hltLCToSCAssociation_cfi import (hltHGCalLCToSCAssociatorByEnergyScoreProducer,
                                                                              hltHGCalLayerClusterSimClusterAssociation)

hltHgcalPrevalidation = cms.Sequence(
    hltHGCalLCToCPAssociatorByEnergyScoreProducer *
    hltHGCalLCToSCAssociatorByEnergyScoreProducer *
    hltHGCalLayerClusterCaloParticleAssociation *
    hltHGCalLayerClusterSimClusterAssociation
)
