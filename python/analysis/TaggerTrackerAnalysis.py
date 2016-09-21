
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

class TaggerTrackerAnalysis(object) : 

    def __init__(self) : 

        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 20})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1

        self.results = []

        self.p = []
        self.pt = []

    def binomial_error(self, total, passing) :
        return (1/total)*math.sqrt(passing*(1-passing/total))


    def process(self, root_file) : 

        print '[ TaggerTrackerAnalysis ]: Processing file: ' + str(root_file)
        self.results = rnp.root2array(root_file) 
        self.p   = self.results['p'][self.results['n_tracks'] != 0]
        self.px  = self.results['px'][self.results['n_tracks'] != 0]
        self.py  = self.results['py'][self.results['n_tracks'] != 0]
        self.pz  = self.results['pz'][self.results['n_tracks'] != 0]
        self.pt  = np.sqrt(np.add(np.power(self.py, 2), np.power(self.pz, 2))) 

    def finalize(self) : 
 
        with PdfPages("tagger_tracker_analysis.pdf") as pdf : 

            # Create the figure
            fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(15, 7.5))

            # Set the number of bins
            bins = np.linspace(0., 4., 251)

            ax0.hist(self.p, bins=bins, histtype='step', linewidth=2)
            ax0.set_xlabel("$p(e^{-})$ GeV")
            ax0.set_yscale("symlog")

            # Set the number of bins
            bins = np.linspace(0.06, .2, 251)
            ax1.hist(self.pt, bins=bins, histtype='step', linewidth=2)
            ax1.set_xlabel("$p_t(e^-)$ GeV")
            ax1.set_yscale("symlog")
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
