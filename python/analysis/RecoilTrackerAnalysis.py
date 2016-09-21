
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

class RecoilTrackerAnalysis(object) : 

    def __init__(self) : 

        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 12})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1

        self.results = []

        self.recoil_p        = []
        self.recoil_truth_p  = []
        self.recoil_pt       = []
        self.recoil_truth_pt = []
        self.n_recoil_hits   = []

    def binomial_error(self, total, passing) :
        return (1/total)*math.sqrt(passing*(1-passing/total))


    def process(self, root_file) : 

        print '[ RecoilTrackerAnalysis ]: Processing file: ' + str(root_file)
        self.results = rnp.root2array(root_file) 
        self.n_tracks = self.results["n_tracks"][self.results['recoil_is_found'] == 1]
        self.recoil_p  = self.results["recoil_p"][self.results['recoil_is_found'] == 1]
        self.recoil_px  = self.results["recoil_px"][self.results['recoil_is_found'] == 1]
        self.recoil_py  = self.results["recoil_py"][self.results['recoil_is_found'] == 1]
        self.recoil_pz  = self.results["recoil_pz"][self.results['recoil_is_found'] == 1]
        self.recoil_pt = self.results["recoil_pt"][self.results['recoil_is_found'] == 1]
        self.recoil_truth_p  = self.results["recoil_truth_p_last"][self.results['recoil_is_found'] == 1]
        self.recoil_truth_pt = self.results["recoil_truth_pt_first"][self.results['recoil_is_found'] == 1]
        self.recoil_truth_px = self.results["recoil_truth_px_first"][self.results['recoil_is_found'] == 1]
        self.recoil_truth_py = self.results["recoil_truth_py_first"][self.results['recoil_is_found'] == 1]
        self.n_recoil_hits = self.results["n_recoil_hits"][self.results['recoil_is_found'] == 1]
        self.n_3d_hits = self.results["n_3d_hits"][self.results['recoil_is_found'] == 1]
        self.recoil_chi2 = self.results['recoil_chi2'][self.results['recoil_is_found'] == 1]
        self.n_recoil_mishits = self.results["n_recoil_mishits"][self.results['recoil_is_found'] == 1]
        self.recoil_cluster_size = self.results['recoil_cluster_size'][self.results['recoil_is_found'] == 1]
        self.recoil_cluster_layer = self.results['recoil_cluster_layer'][self.results['recoil_is_found'] == 1]
        self.hardest_brem_pos_z = self.results['hardest_brem_pos_z'][self.results['recoil_is_found'] == 1]
        self.hardest_brem_energy =  self.results['hardest_brem_energy'][self.results['recoil_is_found'] == 1]
        

    def finalize(self) : 
 
        with PdfPages("recoil_tracker_analysis.pdf") as pdf : 

            #
            # Active target energy
            #
            target_energy = self.results['target_energy']
            #[self.results['target_energy'] >= 0]

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(0, 30, 301)
            ax.hist(target_energy, bins=bins, histtype='step', linewidth=2)
            ax.set_xlabel("Target Energy (MeV)")
            ax.set_yscale("symlog")

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            #
            # Trigger scintillator energy
            #
            trigger_pad_energy = self.results['trigger_pad_energy']
            #[self.results['trigger_pad_energy'] >= 0]

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(0, 30, 301)
            ax.hist(target_energy, bins=bins, histtype='step', linewidth=2)
            ax.set_xlabel("Trigger Pad Energy (MeV)")
            ax.set_yscale("symlog")

            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            #
            #
            energy_cut = 0
            rejection_efficiency = []
            energy_cuts = []
            for index in xrange(0, 80) : 
                rejection_efficiency.append((len(target_energy[target_energy <= energy_cut])/len(target_energy))*100)
                energy_cuts.append(energy_cut)
                energy_cut += 0.05
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            ax.errorbar(energy_cuts, rejection_efficiency, marker='o', markersize=7, linestyle='--')
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Track multiplicity
            #
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(0, 5, 6)
            ax.hist(self.results["n_tracks"], bins=bins, histtype='step', linewidth=2)
            ax.set_xlabel("Track multiplicity")
            ax.set_yscale("symlog")
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Number of hits associated with recoil track
            #
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))

            # Set the number of bins
            bins = np.linspace(0, 8, 9)
            ax.hist(self.n_recoil_hits, bins=bins, histtype='step', linewidth=2)
            ax.set_xlabel("Hits per recoil track")
            ax.set_yscale("symlog")
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Recoil track p and pt
            #


            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))

            bins = np.linspace(0., 300., 301)
            ax0.hist(self.recoil_p, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")

            bins = np.linspace(0, 2.0, 201)
            ax1.hist(self.recoil_pt, bins=bins, histtype='step', linewidth=2)
            ax1.set_xlabel("Recoil $p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")

            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Recoil track p and pt
            #


            fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, figsize=(7.5, 22.5))

            bins = np.linspace(0., 7., 701)
            ax0.hist(self.recoil_px, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel("Recoil $p_{x}(e^{-})$ GeV")
            ax0.set_yscale("symlog")

            bins = np.linspace(-.2, .1, 501)
            ax1.hist(self.recoil_py, bins=bins, histtype='step', linewidth=2)
            ax1.set_xlabel("Recoil $p_y(e^-)$ GeV")
            ax1.set_yscale("symlog")

            bins = np.linspace(-.1, .1, 501)
            ax2.hist(self.recoil_pz, bins=bins, histtype='step', linewidth=2)
            ax2.set_xlabel("Recoil $p_z(e^-)$ GeV")
            ax2.set_yscale("symlog")

            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            bins = np.linspace(0., 7., 701)
            ax0.hist(self.recoil_p, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")

            bins = np.linspace(0, .1, 501)
            ax1.hist(self.recoil_pt, bins=bins, histtype='step', linewidth=2)
            ax1.set_xlabel("Recoil $p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")

            pdf.savefig(bbox_inches='tight')
            plt.close()


            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            bins = np.linspace(0., 7., 701)
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_mishits == 0]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='0 mishits')
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_mishits >= 1]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='>=1 mishits')
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")
            ax0.legend(loc=2)

            print "Total tracks: " + str(len(self.recoil_p))
            print "Total tracks with mishits: " + str(len(recoil_p_hit_cut))

            bins = np.linspace(0, .1, 501)
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_mishits == 0]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='0 mishits')
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_mishits >= 1]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='>= 1 mishits')
            ax1.set_xlabel("$p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")
            ax1.legend(loc=1)

            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            bins = np.linspace(0., 7., 701)
            recoil_p_hit_cut = self.recoil_p[(self.n_recoil_mishits == 0) &
                                            (self.n_tracks == 1)]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='0 mishits, # tracks = 1')
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_mishits >= 1
                                            & (self.n_tracks == 1)]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='>=1 mishits, # tracks = 1')
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")
            ax0.legend(loc=2)

            bins = np.linspace(0, .1, 501)
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_mishits == 0
                                              & (self.n_tracks == 1)]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='0 mishits, # tracks = 1')
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_mishits >= 1
                                              & (self.n_tracks == 1)]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='>= 1 mishits, # tracks = 1')
            ax1.set_xlabel("$p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")
            ax1.legend(loc=1)

            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            bins_x = np.linspace(0, 40, 41)
            bins_y = np.linspace(0, 10, 11)
            ax0.hist2d(self.n_3d_hits, self.n_recoil_mishits, bins=[bins_x, bins_y], norm=LogNorm())
            ax0.set_xlabel("Total stereo hits per event")
            ax0.set_ylabel("Total recoil mishits")
            
            bins_x = np.linspace(0, 100, 101)
            bins_y = np.linspace(0, 10, 11)
            ax1.hist2d(self.recoil_chi2, self.n_recoil_mishits, bins=[bins_x, bins_y], norm=LogNorm())
            ax1.set_xlabel("Total stereo hits per event")
            ax1.set_ylabel("Total recoil mishits")

            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            bins = np.linspace(0., 300., 301)
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_hits == 6]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_hits == 5]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_hits == 4]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")
            ax0.legend(loc=2)

            bins = np.linspace(0, 2.0, 201)
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_hits == 6]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_hits == 5]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_hits == 4]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax1.set_xlabel("$p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")
            ax1.legend(loc=1)

            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            bins = np.linspace(0., 7., 701)
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_hits == 6]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_hits == 5]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            recoil_p_hit_cut = self.recoil_p[self.n_recoil_hits == 4]
            ax0.hist(recoil_p_hit_cut, bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")
            ax0.legend(loc=2)

            bins = np.linspace(0, .1, 501)
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_hits == 6]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_hits == 5]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            recoil_pt_hit_cut = self.recoil_pt[self.n_recoil_hits == 4]
            ax1.hist(recoil_pt_hit_cut, bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax1.set_xlabel("$p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")
            ax1.legend(loc=1)

            pdf.savefig(bbox_inches='tight')
            plt.close()


            #
            #
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            recoil_p_l9_single_10_single = []
            recoil_p_l9_double_10_single = []
            recoil_p_l9_single_10_double = []
            recoil_p_l9_double_10_double = []

            l9_cluster_size = 0
            l10_cluster_size = 0
            for index_i in xrange(0, len(self.recoil_cluster_size)) : 
                for index_j in xrange(0, len(self.recoil_cluster_size[index_i])) :
                    if self.recoil_cluster_layer[index_i][index_j] == 9 : 
                        l9_cluster_size = self.recoil_cluster_size[index_i][index_j]
                    elif self.recoil_cluster_layer[index_i][index_j] == 10 : 
                        l10_cluster_size = self.recoil_cluster_size[index_i][index_j]
                
                if l9_cluster_size == 1 and l10_cluster_size == 1 : 
                    recoil_p_l9_single_10_single.append(self.recoil_p[index_i])
                if l9_cluster_size == 1 and l10_cluster_size == 2 : 
                    recoil_p_l9_single_10_double.append(self.recoil_p[index_i])
                if l9_cluster_size == 2 and l10_cluster_size == 1 : 
                    recoil_p_l9_double_10_single.append(self.recoil_p[index_i])
                if l9_cluster_size == 2 and l10_cluster_size == 2 : 
                    recoil_p_l9_double_10_double.append(self.recoil_p[index_i])

            bins = np.linspace(0., 7., 701)
            ax.hist(recoil_p_l9_single_10_single, bins=bins, histtype='step', linewidth=2, label='single, single')
            ax.hist(recoil_p_l9_double_10_single, bins=bins, histtype='step', linewidth=2, label='double, single')
            ax.hist(recoil_p_l9_single_10_double, bins=bins, histtype='step', linewidth=2, label='single, double')
            ax.hist(recoil_p_l9_double_10_double, bins=bins, histtype='step', linewidth=2, label='double, double')
            ax.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax.set_yscale("symlog")
            ax.legend(loc=2)
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Recoil track p and pt
            #


            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))

            bins = np.linspace(0., 6., 601)
            ax0.hist(self.recoil_truth_p, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel("Recoil $p(e^{-})$ GeV")
            ax0.set_yscale("symlog")

            bins = np.linspace(0, .1, 501)
            ax1.hist(self.recoil_truth_pt, bins=bins, histtype='step', linewidth=2)
            ax1.set_xlabel("Recoil $p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")

            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Recoil track parameters
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(0, 100, 101)
            recoil_chi2 = self.results['recoil_chi2'][self.results['recoil_is_found'] == 1]
            ax.hist(recoil_chi2[self.n_recoil_hits == 6], 
                    bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            ax.hist(recoil_chi2[self.n_recoil_hits == 5], 
                    bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            ax.hist(recoil_chi2[self.n_recoil_hits == 4], 
                    bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax.set_xlabel("Recoil track $\chi^2$")
            ax.set_yscale("symlog")
            ax.legend(loc=1)
            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(-15, 15, 121)
            recoil_d0 = self.results['recoil_d0'][self.results['recoil_is_found'] == 1]
            ax.hist(recoil_d0[self.n_recoil_hits == 6],
                    bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            ax.hist(recoil_d0[self.n_recoil_hits == 5],
                    bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            ax.hist(recoil_d0[self.n_recoil_hits == 4],
                    bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax.set_xlabel("Recoil track D0")
            ax.set_yscale("symlog")
            ax.legend(loc=2)
            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(-15, 15, 121)
            recoil_z0 = self.results['recoil_z0'][self.results['recoil_is_found'] == 1]
            ax.hist(recoil_z0[self.n_recoil_hits == 6],
                    bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            ax.hist(recoil_z0[self.n_recoil_hits == 5],
                    bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            ax.hist(recoil_z0[self.n_recoil_hits == 4],
                    bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax.set_xlabel("Recoil track Z0")
            ax.set_yscale("symlog")
            ax.legend(loc=2)
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(0, 20, 121)
            recoil_phi0 = self.results['recoil_phi0'][self.results['recoil_is_found'] == 1]
            ax.hist(np.sin(recoil_phi0[self.n_recoil_hits == 6])*180/np.pi,
                    bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            ax.hist(np.sin(recoil_phi0[self.n_recoil_hits == 5])*180/np.pi,
                    bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            ax.hist(np.sin(recoil_phi0[self.n_recoil_hits == 4])*180/np.pi,
                    bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax.set_xlabel("Recoil track $\sin{\phi_0}$ (Degrees)")
            ax.set_yscale("symlog")
            ax.legend(loc=1)
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(-0.02, 0.02, 201)
            recoil_omega = self.results['recoil_omega'][self.results['recoil_is_found'] == 1]
            ax.hist(recoil_omega[self.n_recoil_hits == 6],
                    bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            ax.hist(recoil_omega[self.n_recoil_hits == 5],
                    bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            ax.hist(recoil_omega[self.n_recoil_hits == 4],
                    bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax.set_xlabel("Recoil track $\Omega$")
            ax.set_yscale("symlog")
            ax.legend(loc=2)
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            bins = np.linspace(-0.5, 0.5, 101)
            recoil_tan_lambda = self.results['recoil_tan_lambda'][self.results['recoil_is_found'] == 1]
            ax.hist(recoil_tan_lambda[self.n_recoil_hits == 6],
                    bins=bins, histtype='step', linewidth=2, label='6 hit tracks')
            ax.hist(recoil_tan_lambda[self.n_recoil_hits == 5],
                    bins=bins, histtype='step', linewidth=2, label='5 hit tracks')
            ax.hist(recoil_tan_lambda[self.n_recoil_hits == 4],
                    bins=bins, histtype='step', linewidth=2, label='4 hit tracks')
            ax.set_xlabel("Recoil track $\tan{\lambda}$")
            ax.set_yscale("symlog")
            ax.legend(loc=1)
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Recoil track p and pt resolution
            #
            p_res = np.subtract(self.recoil_p, self.recoil_truth_p)

            fig, axes = plt.subplots(nrows=16, ncols=1, figsize=(7.5, 120))
            ax = axes.flatten() 
            
            bins = np.linspace(-2, 2, 100)
            bin_edge = 0
            p_res_points = []
            p_bin_centers = []
            for index in xrange(0, 16) : 
                p_res_region = p_res[(self.recoil_truth_p < (bin_edge + 0.5)) 
                                     & (self.recoil_truth_p >= bin_edge)]
                if len(p_res_region) != 0 :
                    ax[index].hist(p_res_region, bins=bins, histtype='step', linewidth=2, normed=True)
                
                    mu, std = norm.fit(p_res_region[np.abs(p_res_region) < .4])
                    p_res_points.append((std/(bin_edge + 0.5/2))*100)
                    p = norm.pdf(bins, mu, std)
                    ax[index].plot(bins, p, 'k', linewidth=2)

                    ax[index].set_xlabel("Recoil track $p(e^-)$ - truth $p(e^-)$")
                    p_bin_centers.append(bin_edge + 0.5/2)
                bin_edge += 0.5
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Recoil track p and pt resolution
            #
            pt_res = np.subtract(self.recoil_pt, self.recoil_truth_pt)

            fig, axes = plt.subplots(nrows=16, ncols=1, figsize=(7.5, 120))
            ax = axes.flatten() 
            
            bins = np.linspace(-.05, .05, 101)
            bin_edge = 0
            pt_res_points = []
            pt_bin_centers = []
            for index in xrange(0, 16) : 
                pt_res_region = pt_res[(self.recoil_truth_p < (bin_edge + 0.5)) 
                                     & (self.recoil_truth_p >= bin_edge)]
                if len(pt_res_region) != 0 :
                    ax[index].hist(pt_res_region, bins=bins, histtype='step', linewidth=2, normed=True)
                
                    mean = pt_res_region.mean()
                    std_dev = pt_res_region.std()

                    mu, std = norm.fit(pt_res_region[(pt_res_region > (mean - 2*std_dev))
                                                     (pt_res_region < (mean + 2*std_dev))])
                    pt_res_points.append(std)
                    p = norm.pdf(bins, mu, std)
                    ax[index].plot(bins, p, 'k', linewidth=2)

                    ax[index].set_xlabel("Recoil track $p_t(e^-)$ - truth $p_t(e^-)$")
                    pt_bin_centers.append(bin_edge + 0.5/2)
                bin_edge += 0.5
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

 
            #
            # Recoil track p and pt resolution
            #
            '''
            px_res = np.subtract(self.recoil_px, self.recoil_truth_px)

            fig, axes = plt.subplots(nrows=16, ncols=1, figsize=(7.5, 120))
            ax = axes.flatten() 
            
            bins = np.linspace(-1, 1, 101)
            bin_edge = 0
            px_res_points = []
            px_bin_centers = []
            for index in xrange(0, 16) : 
                px_res_region = px_res[(self.recoil_truth_p < (bin_edge + 0.5)) 
                                     & (self.recoil_truth_p >= bin_edge)]
                if len(px_res_region) != 0 :
                    ax[index].hist(px_res_region, bins=bins, histtype='step', linewidth=2, normed=True)
                
                    mu, std = norm.fit(px_res_region[np.abs(px_res_region) < .15])
                    px_res_points.append(std)
                    p = norm.pdf(bins, mu, std)
                    ax[index].plot(bins, p, 'k', linewidth=2)

                    ax[index].set_xlabel("Recoil track $p_x(e^-)$ - truth $p_x(e^-)$")
                    px_bin_centers.append(bin_edge + 0.5/2)
                bin_edge += 0.5
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, axes = plt.subplots(nrows=16, ncols=1, figsize=(7.5, 120))
            ax = axes.flatten() 
            
            bins = np.linspace(-.05, .05, 101)
            bin_edge = 0
            py_res_points = []
            py_bin_centers = []
            for index in xrange(0, 16) : 
                py_res_region = py_res[(self.recoil_truth_p < (bin_edge + 0.5)) 
                                     & (self.recoil_truth_p >= bin_edge)]
                if len(py_res_region) != 0 :
                    ax[index].hist(py_res_region, bins=bins, histtype='step', linewidth=2, normed=True)
                
                    mu, std = norm.fit(py_res_region[np.abs(px_res_region) < .15])
                    py_res_points.append(std)
                    p = norm.pdf(bins, mu, std)
                    ax[index].plot(bins, p, 'k', linewidth=2)

                    ax[index].set_xlabel("Recoil track $p_y(e^-)$ - truth $p_y(e^-)$")
                    py_bin_centers.append(bin_edge + 0.5/2)
                bin_edge += 0.5

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 15))
            
            ax0.plot(px_bin_centers, px_res_points, marker='o', markersize=10, linestyle='--')
            ax0.set_xlabel("Recoil track truth $px(e^-)$  GeV")
            ax0.set_ylabel("Recoil track $\sigma_{px(e^-)}$")
            
            ax1.plot(py_bin_centers, py_res_points, marker='o', markersize=10, linestyle='--')
            ax1.set_xlabel("Recoil track truth $py(e^-)$  GeV")
            ax1.set_ylabel("Recoil track $\sigma_{py(e^-)}$")
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
            '''

            fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2, figsize=(15, 15))
            
            bins = np.linspace(0., 7., 701)
            im = ax0.hist2d(self.recoil_p, self.recoil_truth_p, bins=bins, norm=LogNorm())
            fig.colorbar(im[3], ax=ax0, extend='min') 
            ax0.set_xlabel("Recoil track $p(e^-)$ (GeV)")
            ax0.set_ylabel("Recoil track truth $p(e^-)$ (GeV)")

            ax1.plot(p_bin_centers, p_res_points, marker='o', markersize=10, linestyle='--')
            ax1.set_xlabel("Recoil track truth $p(e^-)$  GeV")
            ax1.set_ylabel("Recoil track $\sigma_{p(e^-)}$")

            pt_res = np.subtract(self.recoil_pt, self.recoil_truth_pt)

            bins = np.linspace(0, .1, 251)
            im = ax2.hist2d(self.recoil_pt, self.recoil_truth_pt, bins=bins, norm=LogNorm())
            fig.colorbar(im[3], ax=ax2, extend='min')  
            ax2.set_xlabel("Recoil track $p_t(e^-)$")
            ax2.set_ylabel("Recoil track truth $p_t(e^-)$")

            ax3.plot(pt_bin_centers, pt_res_points, marker='o', markersize=10, linestyle='--')
            ax3.set_xlabel("Recoil track truth $p_t(e^-)$  GeV")
            ax3.set_ylabel("Recoil track $\sigma_{p_t(e^-)}$")

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            #
            # Tracking efficiency
            #
            
            bins = np.linspace(0, 4, 41)
            
            recoil_truth_p_find = self.results['recoil_truth_p_last'][self.results['recoil_is_findable'] == 1]
            print "Findable tracks: " + str(len(recoil_truth_p_find))
            p_truth_find_hist, bin_edges = np.histogram(recoil_truth_p_find, bins=bins)
            print p_truth_find_hist

            recoil_truth_p_found = self.results['recoil_truth_p_last'][
                (self.results['recoil_is_findable'] == 1) & (self.results['recoil_is_found'] == 1)]
            print "Found tracks: " + str(len(recoil_truth_p_found))
            p_truth_found_hist, bin_edges = np.histogram(recoil_truth_p_found, bins=bins)


            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.

            if len(p_truth_find_hist) != len(p_truth_found_hist) : print "Arrays are not the same size"
            
            recoil_trk_eff = []
            y_err = []
            for index in xrange(0, len(p_truth_find_hist)) : 
                recoil_trk_eff.append(float(p_truth_found_hist[index]/p_truth_find_hist[index])*100)
                y_err.append(self.binomial_error(p_truth_find_hist[index], p_truth_found_hist[index])*100)
            #recoil_trk_eff = np.divide(p_truth_found_hist, p_truth_find_hist)
            #print recoil_trk_eff

         
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            ax.errorbar(bin_centers, recoil_trk_eff, yerr=y_err, marker='o', markersize=7, linestyle='--')
            ax.set_xlabel('Recoil track $p(e^{-})$ GeV')
            ax.set_ylabel('Tracking Efficiency (%)')

            pdf.savefig(bbox_inches='tight')
            plt.close()
            #
            #
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bins_x = np.linspace(0, 5, 251)
            bins_y = np.linspace(0, 20, 401)
           
            recoil_p_cut = self.recoil_p[self.hardest_brem_pos_z > -.175]
            hardest_brem_pos_z_cut = self.hardest_brem_pos_z[self.hardest_brem_pos_z > -.175]
            im = ax.hist2d(recoil_p_cut, hardest_brem_pos_z_cut, bins=[bins_x, bins_y], norm=LogNorm())
            fig.colorbar(im[3], ax=ax, extend='min')  
            ax.set_xlabel("Recoil track $p(e^-)$")
            ax.set_ylabel("Hardest brem z (cm)")

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bins = np.linspace(0, 5, 251)
            
            hardest_brem_energy_cut = self.hardest_brem_energy[self.hardest_brem_pos_z > -.175]
            im = ax.hist2d(recoil_p_cut, hardest_brem_energy_cut, bins=bins, norm=LogNorm())
            fig.colorbar(im[3], ax=ax, extend='min')  
            ax.set_xlabel("Recoil track $p(e^-)$")
            ax.set_ylabel("Hardest brem energy (GeV)")

            pdf.savefig(bbox_inches='tight')
            plt.close()


            #
            #
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bins_x = np.linspace(0, 5, 251)
            bins_y = np.linspace(0, 20, 401)
           
            recoil_p_cut = self.recoil_p[(self.hardest_brem_pos_z > -.175)
                                        & (self.hardest_brem_pos_z < .175)]
            hardest_brem_pos_z_cut = self.hardest_brem_pos_z[(self.hardest_brem_pos_z > -.175)
                                        & (self.hardest_brem_pos_z < .175)]
            im = ax.hist2d(recoil_p_cut, hardest_brem_pos_z_cut, bins=[bins_x, bins_y], norm=LogNorm())
            fig.colorbar(im[3], ax=ax, extend='min')  
            ax.set_xlabel("Recoil track $p(e^-)$")
            ax.set_ylabel("Hardest brem z (cm)")

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bins = np.linspace(0, 5, 251)
            
            hardest_brem_energy_cut = self.hardest_brem_energy[(self.hardest_brem_pos_z > -.175)
                                        & (self.hardest_brem_pos_z < .175)]
            im = ax.hist2d(recoil_p_cut, hardest_brem_energy_cut, bins=bins, norm=LogNorm())
            fig.colorbar(im[3], ax=ax, extend='min')  
            ax.set_xlabel("Recoil track $p(e^-)$")
            ax.set_ylabel("Hardest brem energy (GeV)")

            pdf.savefig(bbox_inches='tight')
            plt.close()
