
from __future__ import division

import copy
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp
import ROOT as r
import rootpy.plotting.root2matplotlib as rplt

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from scipy.stats import norm

class PhotoNuclearAnalysis(object) : 

    def __init__(self) : 

        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 20})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1

        self.initialize()

    def initialize(self) :

        self.events = []
        self.pn_gamma_energy = []
        self.pn_particle_mult = []
        self.total_energy_ecal = []
        self.total_energy_hcal = []
        self.neutral_hadron_lead = []
        self.target_pn_energy = []
        self.target_recoil_energy = []

    def process(self, event) :
        
        sim_particles = event.getCollection("SimParticles")

        recoil_e = None
        pn_gamma = None

        for particle_n in xrange(0, sim_particles.GetEntriesFast()) : 
            particle = sim_particles.At(particle_n)
            
            print "[ PhotoNuclearAnalysis ]: Particle PDG ID: %s" % particle.getPdgID()
            print "[ PhotoNuclearAnalysis ]: Parent count: %s" % particle.getParentCount()
            print "[ PhotoNuclearAnalysis ]: Daughter count: %s" % particle.getDaughterCount()
            print "[ PhotoNuclearAnalysis ]: Generator Status: %s" % particle.getGenStatus()
           
            if self.is_pn_gamma(particle) :
                pn_gamma = particle
                break
        
        if pn_gamma : 

            #print "[ PhotoNuclearAnalysis ]: Particle PDG ID: %s" % particle.getPdgID()
            #print "[ PhotoNuclearAnalysis ]: Parent count: %s" % particle.getParentCount()
            #print "[ PhotoNuclearAnalysis ]: Daughter count: %s" % particle.getDaughterCount()
            #print "[ PhotoNuclearAnalysis ]: Generator Status: %s" % particle.getGenStatus()
            self.pn_gamma_energy.append(pn_gamma.getEnergy())
            self.pn_particle_mult.append(pn_gamma.getDaughterCount())

        ecal_hits = event.getCollection("EcalSimHits")
        total_ecal_energy = 0
        for ecal_hit_n in xrange(0, ecal_hits.GetEntriesFast()) : 
            ecal_hit = ecal_hits.At(ecal_hit_n)
            total_ecal_energy += ecal_hit.getEdep()

        self.total_energy_ecal.append(total_ecal_energy)

        hcal_hits = event.getCollection("HcalSimHits")
        total_hcal_energy = 0
        for hcal_hit_n in xrange(0, hcal_hits.GetEntriesFast()) : 
            hcal_hit = hcal_hits.At(hcal_hit_n)
            total_hcal_energy += hcal_hit.getEdep()

        self.total_energy_hcal.append(total_hcal_energy)

        target_hits = event.getCollection("TargetSimHits")
        total_pn_target_energy = 0
        for index in xrange(0, target_hits.GetEntriesFast()) : 
            target_hit = target_hits.At(index)
            if target_hit.getSimParticle().getParentCount() != 0 : 
                if target_hit.getSimParticle().getParent(0) == pn_gamma :
                    total_pn_target_energy += target_hit.getEdep()

        self.target_pn_energy.append(total_pn_target_energy)

    def finalize(self) : 
        gev = 1/1000

        with PdfPages("photo_nuclear_analysis.pdf") as pdf : 

            #
            # Photonuclear gamma energy
            #

            # Create the figure
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            bins = np.linspace(0, 5, 51)
            ax0.hist(np.array(self.pn_gamma_energy)*gev, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel('$E(\gamma)$ (GeV)')
            ax0.set_yscale('symlog')
            pdf.savefig(bbox_inches='tight')
            plt.close()
           
            #
            # PN Multiplicity
            #
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            bins = np.linspace(0, 100, 101)
            ax0.hist(self.pn_particle_mult, bins, histtype='step', linewidth=2)
            ax0.set_xlabel('PN Multiplicity')
            ax0.set_yscale('symlog')
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # PN Multiplicity vs Photonuclear gamma energy
            #

            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            
            bins_x = np.linspace(0, 5, 51)
            bins_y = np.linspace(0, 100, 101)
            im = ax0.hist2d(np.array(self.pn_gamma_energy)*gev,
                            self.pn_particle_mult, bins=[bins_x, bins_y],  norm=LogNorm())
            fig.colorbar(im[3], ax=ax0) 
            ax0.set_xlabel('$E(\gamma)$ (GeV)')
            ax0.set_ylabel('PN Multiplicity')
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            bins = np.linspace(0, 200, 101)
            ax0.hist(self.total_energy_ecal, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel('Total energy deposited in Ecal (MeV)')
            ax0.set_yscale('symlog')
            pdf.savefig(bbox_inches='tight')
            plt.close()
        
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            bins = np.linspace(0, 200, 101)
            ax0.hist(self.total_energy_hcal, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel('Total energy deposited in Hcal (MeV)')
            ax0.set_yscale('symlog')
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            im = ax0.hist2d(self.total_energy_ecal,
                            self.total_energy_hcal, bins=bins,  norm=LogNorm())
            ax0.set_xlabel('Total energy deposited in Ecal (MeV)')
            ax0.set_ylabel('Total energy deposited in Hcal (MeV)')
            pdf.savefig(bbox_inches='tight')
            plt.close()
          
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            bins = np.linspace(0, 200, 101)
            ax0.hist(self.target_pn_energy, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel('PN energy in target (MeV)')
            ax0.set_yscale('symlog')
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # PN Multiplicity
            #
            #fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
            #bins = np.linspace(0, 100, 101)
            #ax0.hist(self.pn_particle_mult, bins, histtype='step', linewidth=2)
            #ax0.hist(np.array(self.pn_particle_mult)[
            #        (np.array(self.total_energy_ecal) == 0) &
            #        (np.array(self.total_energy_hcal) == 0)]
            #        , bins, histtype='step', linewidth=2)
            #ax0.set_xlabel('PN Multiplicity')
            #ax0.set_yscale('symlog')
            #pdf.savefig(bbox_inches='tight')
            #plt.close()


    def is_recoil(self, sim_particle) : 
        return ((sim_particle.getPdgID() == 11) &
               (sim_particle.getParentCount() == 0))

    def is_pn_gamma(self, sim_particle) : 
        return ((sim_particle.getPdgID() == 22) &
               (sim_particle.getParentCount() == 0))

    def created_within_target(self, sim_particle) : 
        
        if abs(sim_particle.getVertex()[2]) < 0.175 : return True 
        return False
