
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

class SignalAnalysis(object) : 
    
    def __init__(self) : 

        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 12})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1

        self.results = []
   
    def binomial_error(self, total, passing) :
        return (1/total)*math.sqrt(passing*(1-passing/total))

    def process(self, root_file) : 
        # Load the records containing the A' MC mass resolution measurements
        path = "/home/omoreno/work/ldmx/analysis/signal_mc/eventNumbers.txt"
        self.event_numbers = np.genfromtxt(path,
                              dtype=[('event', 'f8'), 
                                     ('triggered', 'f8'),
                                     ('ecal_energy', 'f8')],
                              delimiter=" "
        )

        self.results = rnp.root2array(root_file)
        self.recoil_is_found = self.results['recoil_is_found']
        self.hardest_brem_energy = self.results['hardest_brem_energy']
        
        self.recoil_p = self.results["recoil_p"][
            (self.recoil_is_found == 1)
            & (self.hardest_brem_energy == -10000)
        ]
        
        self.recoil_pt = self.results["recoil_pt"][
            (self.recoil_is_found == 1)
            & (self.hardest_brem_energy == -10000)
        ]

        self.recoil_truth_p = self.results["recoil_truth_p"]
        self.recoil_truth_py = self.results["recoil_truth_px"]
        self.recoil_truth_pz = self.results["recoil_truth_pz"]
        self.recoil_truth_pt = self.results["recoil_truth_pt"]
        self.recoil_is_findable = self.results['recoil_is_findable']
        self.theta = np.abs(np.arctan(np.divide(self.recoil_truth_py, self.recoil_truth_pz))*(180/math.pi))
        self.ap_mass = self.results['ap_mass']
        #self.recoil_ecal_sp_x = self.results['recoil_ecal_sp_x']
        self.recoil_vertex_x = self.results['recoil_vertex_x']
        self.x_target = self.results['x_target']
        self.recoil_vertex_y = self.results['recoil_vertex_y']
        self.y_target = self.results['y_target']


    def finalize(self) :

        # Get a unique set of A' masses
        ap_masses = np.unique(self.ap_mass)

        tracker_acceptance = []
        tracker_acceptance_err = []

        tracker_acceptance_energy_cut = []
        tracker_acceptance_energy_cut_err = []

        combined_acceptance = []
        combined_acceptance_err = []
        
        combined_acceptance_transverse = []
        combined_acceptance_transverse_err = []
        
        masses = []
        for mass in ap_masses : 
            print "[ Signal Analysis ]: Processing mass: " + str(mass)

            ap_acceptance = self.recoil_is_findable[self.ap_mass == mass]
            
            tracker_acceptance.append(len(ap_acceptance[ap_acceptance == 1])/len(ap_acceptance))
            tracker_acceptance_err.append(self.binomial_error(len(ap_acceptance),
                                                              len(ap_acceptance[ap_acceptance == 1])))

            ap_recoil_truth_p = self.recoil_truth_p[self.ap_mass == mass]
            ap_acceptance_energy_cut = ap_acceptance[ap_recoil_truth_p < 1.2]

            #ap_recoil_vertex_x = self.recoil_vertex_x[self.ap_mass == mass]
            #ap_recoil_ecal_sp_x = self.recoil_ecal_sp_x[self.ap_mass == mass]
            #dx = np.subtract(ap_recoil_vertex_x, ap_recoil_ecal_sp_x)
            #dx_energy = dx[ap_recoil_truth_p < 1.2]
            #print dx_energy

            tracker_acceptance_energy_cut.append(
                len(ap_acceptance_energy_cut[ap_acceptance_energy_cut == 1])/len(ap_acceptance))
            tracker_acceptance_energy_cut_err.append(
                self.binomial_error(len(ap_acceptance),
                                    len(ap_acceptance_energy_cut[ap_acceptance_energy_cut == 1])))

            '''
            if (int(mass*1000) == 200) or (int(mass*1000) == 1) : 
                print "Skipping mass " + str(mass)
                continue

            masses.append(mass)

            ap_events = self.results['event'][mass == self.ap_mass]
            
            ap_ecal_selection_all = self.event_numbers['triggered'][
                (self.event_numbers['event'] >= ap_events[0])
                & (self.event_numbers['event'] <= ap_events[len(ap_events) - 1])]
            
            ap_ecal_selection_event = self.event_numbers['event'][
                (self.event_numbers['event'] >= ap_events[0])
                & (self.event_numbers['event'] <= ap_events[len(ap_events) - 1])]
 
            ap_ecal_selection = []
            for index in xrange(0, len(ap_ecal_selection_event)) : 
                if ap_ecal_selection_event[index] in ap_events : 
                    ap_ecal_selection.append(ap_ecal_selection_all[index])

            combined = 0
            combined_trans = 0

            for index in xrange(0, len(ap_acceptance_energy_cut)) :
                
                ecal_selection = ap_ecal_selection[index] 

                if ecal_selection.size == 0 : continue 
                
                if (ecal_selection == 1) &  (ap_acceptance_energy_cut[index] == 1): 
                    combined += 1
                    if dx_energy[index] > 40 : 
                        combined_trans += 1

            combined_acceptance.append(combined/len(ap_acceptance))
            combined_acceptance_err.append(self.binomial_error(len(ap_acceptance), combined))

            combined_acceptance_transverse.append(combined_trans/len(ap_acceptance))
            combined_acceptance_transverse_err.append(self.binomial_error(len(ap_acceptance), combined_trans))
            '''

        tracker_acceptance = np.array(tracker_acceptance)*100
        tracker_acceptance_err = np.array(tracker_acceptance_err)*100

        tracker_acceptance_energy_cut = np.array(tracker_acceptance_energy_cut)*100
        tracker_acceptance_energy_cut_err = np.array(tracker_acceptance_energy_cut_err)*100
        
        '''
        combined_acceptance = np.array(combined_acceptance)*100
        combined_acceptance_err = np.array(combined_acceptance_err)*100
        '''
        with PdfPages("signal_analysis.pdf") as pdf : 
            
            plot_masses = [10, 100, 200, 400, 1000, 1500]
            
            #
            # Kinematics
            #
            fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, figsize=(7.5, 22.5))

            bins_p  = np.linspace(0, 4, 120)
            bins_pt = np.linspace(0, 1.5, 90)
            bins_theta = np.linspace(0, 90, 180)

            for mass in ap_masses :

                if int(round(mass*1000, 0)) not in plot_masses : continue 

                ax0.hist(self.recoil_truth_p[self.ap_mass == mass], bins=bins_p, histtype='step', 
                         linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
                ax1.hist(self.recoil_truth_pt[self.ap_mass == mass], bins=bins_pt, alpha=0.7,
                         histtype='step', linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
                ax2.hist(self.theta[self.ap_mass == mass], 
                         bins=bins_theta, alpha=0.7, histtype='step', linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
   
            ax0.set_xlabel("$p(e^{-})$ (GeV)")
            ax0.set_yscale("log")
            ax0.legend(fontsize=10, ncol=2)

            ax1.set_xlabel("$p_{t}(e^{-})$ (GeV)")
            ax1.set_yscale("log")
            ax1.legend(fontsize=10)

            ax2.set_xlabel("$\theta(e^{-})$ (degrees)")
            ax2.set_yscale("log")
            ax2.legend(fontsize=10, ncol=2)

            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Kinematics within acceptance
            #
            fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, figsize=(7.5, 22.5))

            bins_p  = np.linspace(0, 4, 120)
            bins_pt = np.linspace(0, 1.5, 90)
            bins_theta = np.linspace(0, 90, 180)

            for mass in ap_masses :

                if int(round(mass*1000, 0)) not in plot_masses : continue 
                
                ax0.hist(self.recoil_truth_p[
                                (self.ap_mass == mass) & (self.recoil_is_findable == 1)],
                         bins=bins_p, histtype='step', 
                         linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
                ax1.hist(self.recoil_truth_pt[
                                (self.ap_mass == mass) & (self.recoil_is_findable == 1)],
                         bins=bins_pt, alpha=0.7,
                         histtype='step', linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
                ax2.hist(self.theta[(self.ap_mass == mass) & (self.recoil_is_findable == 1)],
                         bins=bins_theta, alpha=0.7, histtype='step', linewidth=2,
                         normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
    
            ax0.set_title("Within Tracker Acceptance")
            ax0.set_xlabel("$p(e^{-})$ (GeV)")
            ax0.set_yscale("log")
            ax0.legend(fontsize=10, ncol=2)

            ax1.set_xlabel("$p_{t}(e^{-})$ (GeV)")
            ax1.set_yscale("log")
            ax1.legend(fontsize=10)

            ax2.set_xlabel("$\theta(e^{-})$ (degrees)")
            ax2.set_yscale("log")
            ax2.legend(fontsize=10, ncol=2)

            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Kinematics within acceptance
            #
            fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, figsize=(7.5, 22.5))

            bins_p  = np.linspace(0, 4, 120)
            bins_pt = np.linspace(0, 1.5, 90)
            bins_theta = np.linspace(0, 90, 180)

            for mass in ap_masses :

                if int(round(mass*1000, 0)) not in plot_masses : continue 
                
                ax0.hist(self.recoil_truth_p[
                                (self.ap_mass == mass) & (self.recoil_is_findable == 1)
                                & (self.recoil_truth_p < 1.2)],
                         bins=bins_p, histtype='step', normed=True, 
                         linewidth=2, label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
                ax1.hist(self.recoil_truth_pt[
                                (self.ap_mass == mass) & (self.recoil_is_findable == 1)
                                & (self.recoil_truth_p < 1.2)],
                         bins=bins_pt, alpha=0.7,
                         histtype='step', linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
                ax2.hist(self.theta[
                                (self.ap_mass == mass) & (self.recoil_is_findable == 1)
                                & (self.recoil_truth_p < 1.2)],
                         bins=bins_theta, alpha=0.7, histtype='step', linewidth=2, normed=True,
                         label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
    
            ax0.set_title("Within Tracker Acceptance and recoil $e^-$ < 1.2 GeV")
            ax0.set_xlabel("$p(e^{-})$ (GeV)")
            ax0.set_yscale("log")
            ax0.legend(fontsize=10, ncol=2)

            ax1.set_xlabel("$p_{t}(e^{-})$ (GeV)")
            ax1.set_yscale("log")
            ax1.legend(fontsize=10)

            ax2.set_xlabel("$\theta(e^{-})$ (degrees)")
            ax2.set_yscale("log")
            ax2.legend(fontsize=10, ncol=2)

            pdf.savefig(bbox_inches='tight')
            plt.close()

            '''
            #
            # Transverse seperation
            #
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bins = np.linspace(-500, 500, 251)
            diff_map = {}
            for mass in ap_masses :
                
                if int(round(mass*1000, 0)) not in plot_masses : continue 

                ap_recoil_ecal_sp_x = self.recoil_ecal_sp_x[
                   (self.ap_mass == mass) & (self.recoil_is_findable == 1)]
                ap_recoil_vertex_x = self.recoil_vertex_x[
                   (self.ap_mass == mass) & (self.recoil_is_findable == 1)]
                
                diff = np.subtract(ap_recoil_vertex_x, ap_recoil_ecal_sp_x)
                diff_map[mass] = diff

                ax0.hist(diff, bins=bins, histtype='step', normed=True, 
                         linewidth=2, label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
            ax0.set_xlabel("Transverse Seperation (mm)")
            ax0.set_yscale("log")
            ax0.legend(fontsize=10, loc=2)

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            for mass in ap_masses :
                
                if int(round(mass*1000, 0)) not in plot_masses : continue 
                
                dx_cut = 0
                efficiency = []
                cuts = []
                for index in xrange(0, 80) : 
                    diff = diff_map[mass]
                    efficiency.append((len(diff[np.abs(diff) >= dx_cut])/len(diff))*100)
                    cuts.append(dx_cut)
                    dx_cut += 2

                ax0.errorbar(cuts, efficiency, marker='o', markersize=7, linestyle='--',
                            label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
            ax0.set_xlabel("Transverse Seperation Cut (mm)")
            ax0.set_ylabel("Efficiency (%)")
            ax0.legend(fontsize=10, loc=1)


            pdf.savefig(bbox_inches='tight')
            plt.close()

            #
            # Transverse seperation
            #
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bins = np.linspace(-500, 500, 251)
            diff_map = {}
            for mass in ap_masses :
                
                if int(round(mass*1000, 0)) not in plot_masses : continue 

                ap_recoil_ecal_sp_x = self.recoil_ecal_sp_x[
                   (self.ap_mass == mass) & (self.recoil_is_findable == 1)
                   & (self.recoil_truth_p  < 1.2)
                ]
                ap_recoil_vertex_x = self.recoil_vertex_x[
                   (self.ap_mass == mass) & (self.recoil_is_findable == 1)
                   & (self.recoil_truth_p  < 1.2)
                ]
                
                diff = np.subtract(ap_recoil_vertex_x, ap_recoil_ecal_sp_x)
                diff_map[mass] = diff

                ax0.hist(diff, bins=bins, histtype='step', normed=True, 
                         linewidth=2, label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
            ax0.set_title("Within tracker acceptance and recoil $e^-$ < 1.2 GeV")
            ax0.set_xlabel("Transverse Seperation (mm)")
            ax0.set_yscale("log")
            ax0.legend(fontsize=10, loc=2)

            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            for mass in ap_masses :
                
                if int(round(mass*1000, 0)) not in plot_masses : continue 
                
                dx_cut = 0
                efficiency = []
                cuts = []
                for index in xrange(0, 80) : 
                    diff = diff_map[mass]
                    efficiency.append((len(diff[np.abs(diff) >= dx_cut])/len(diff))*100)
                    cuts.append(dx_cut)
                    dx_cut += 2

                ax0.errorbar(cuts, efficiency, marker='o', markersize=7, linestyle='--',
                            label="$A'$ mass = " + str(int(round(mass*1000, 0))) + " MeV")
            ax0.set_title("Within tracker acceptance and recoil $e^-$ < 1.2 GeV")
            ax0.set_xlabel("Transverse Seperation Cut (mm)")
            ax0.set_ylabel("Efficiency (%)")
            ax0.legend(fontsize=10, loc=1)


            pdf.savefig(bbox_inches='tight')
            plt.close()
            '''


            #
            # Acceptance
            #

            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))

            ax0.errorbar(np.array(ap_masses)*1000, tracker_acceptance, yerr=tracker_acceptance_err,
                         fmt='', marker='o', markersize=6,
                         linestyle='', label='Tracker Only')
            ax0.set_xlabel("$A'$ Mass (MeV)")
            ax0.set_ylabel("Acceptance (%)")
            ax0.set_xscale("symlog")
            ax0.legend()
            
            pdf.savefig(bbox_inches='tight')
            plt.close()


            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))

            ax0.errorbar(np.array(ap_masses)*1000, tracker_acceptance_energy_cut, 
                         yerr=tracker_acceptance_energy_cut_err,
                         fmt='', marker='o', markersize=6,
                         linestyle='', label='Tracker Only & p < 1.2 GeV')
            ax0.set_xlabel("$A'$ Mass (MeV)")
            ax0.set_ylabel("Acceptance (%)")
            ax0.set_xscale("symlog")
            ax0.legend()
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
           
            bins = np.linspace(0, 4, 41)
            
            recoil_truth_p_find = self.results['recoil_truth_p'][self.results['recoil_is_findable'] == 1]
            print "Findable tracks: " + str(len(recoil_truth_p_find))
            p_truth_find_hist, bin_edges = np.histogram(recoil_truth_p_find, bins=bins)
            print p_truth_find_hist

            recoil_truth_p_found = self.results['recoil_truth_p'][
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

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))

            ax0.errorbar(np.array(ap_masses)*1000, tracker_acceptance, yerr=tracker_acceptance_err,
                         fmt='', marker='o', markersize=6,
                         linestyle='', label='Tracker Only')
            ax0.set_xlabel("$A'$ Mass (MeV)")
            ax0.set_ylabel("Acceptance (%)")
            ax0.set_xscale("symlog")

            ax1.errorbar(bin_centers, recoil_trk_eff, yerr=y_err, marker='o', markersize=7, linestyle='--')
            ax1.set_xlabel('Recoil track $p(e^{-})$ GeV')
            ax1.set_ylabel('Tracking Efficiency (%)')

            pdf.savefig(bbox_inches='tight')
            plt.close()

            '''
            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))

            ax0.errorbar(np.array(masses)*1000, combined_acceptance, 
                         yerr=combined_acceptance_err,
                         fmt='', marker='o', markersize=6,
                         linestyle='', label='Tracker Only & p < 1.2 GeV & Ecal')
            ax0.set_xlabel("$A'$ Mass (MeV)")
            ax0.set_ylabel("Acceptance (%)")
            ax0.set_xscale("symlog")
            ax0.legend()
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))

            ax0.errorbar(np.array(masses)*1000, combined_acceptance_transverse, 
                         yerr=combined_acceptance_transverse_err,
                         fmt='', marker='o', markersize=6,
                         linestyle='', label='Tracker Only,p < 1.2 GeV, Ecal, dx > 4 cm')
            ax0.set_xlabel("$A'$ Mass (MeV)")
            ax0.set_ylabel("Acceptance (%)")
            ax0.set_xscale("symlog")
            ax0.legend()
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
            '''



            ''' 
            p_res = np.subtract(self.recoil_p, self.recoil_truth_p)
            bins_x = np.linspace(0, 4, 81)
            bins_y = np.linspace(-2, 2, 160)
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            ax.hist2d(self.recoil_truth_p, p_res, bins=[bins_x, bins_y], norm=LogNorm())
            ax.set_xlabel("Recoil p($e^-$) GeV")
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
            
            #
            # Recoil track p and pt resolution
            #

            fig, axes = plt.subplots(nrows=16, ncols=1, figsize=(7.5, 120))
            ax = axes.flatten() 
            
            bins = np.linspace(-1, 1, 80)
            bin_edge = 0
            #p_res_std = [1.20997e-02, 4.29245e-02, 7.76853e-02, 1.05983e-01, 
            #             1.43243e-01, 1.75590e-01, 2.18814e-01, 2.38578e-01]
            p_res_std = [6.18790e-02, 1.24762e-01, 1.88175e-01, 2.30077e-01]
            #p_bin_centers = [0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75 ]
            p_bin_centers = [1.0, 2.0, 3.0, 3.5]
            p_res_std = np.divide(p_res_std, p_bin_centers)*100
            '''
            '''
            for index in xrange(0, 16) : 
                region = p_res[
                    (self.recoil_truth_p < (bin_edge + 0.5))
                    & (self.recoil_truth_p >= bin_edge)]
                
                if len(region) != 0 :
                    ax[index].hist(region, bins=bins, histtype='step',
                                   linewidth=2, normed=True)
                    ax[index].set_xlabel("Recoil track $p(e^-)$ - truth $p(e^-)$")
            
                    mean = region.mean()
                    std_dev = region.std()

                    mu, std = norm.fit(region[
                        (region > (mean - 1*std_dev))
                        & (region < (mean + 1*std_dev))
                    ])
                    #[(region > (mean - 2*std_dev))
                    #                                 (region < (mean + 2*std_dev))])

                    p = norm.pdf(bins, mu, std)
                    ax[index].plot(bins, p, 'k', linewidth=2)
                    
                    p_res_points.append((std/(bin_edge + 0.5/2))*100)
                    p_bin_centers.append(bin_edge + 0.5/2)
                
                bin_edge += 0.5
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))
            
            ax0.plot(p_bin_centers, p_res_std, marker='o', markersize=10, linestyle='--')
            ax0.set_xlabel("Recoil track truth $p(e^-)$  GeV")
            ax0.set_ylabel("Recoil track $\sigma_{p(e^-)}$")

            pdf.savefig(bbox_inches='tight')
            plt.close()

            '''
