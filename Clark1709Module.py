########################################################################
# Boolean network simulator v2.0
# Erik Clark
########################################################################

import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
import matplotlib.lines as mlines
from matplotlib import animation
from IPython.display import display, clear_output

########################################################################

class Gene(object):
    """
    Holds control logic
    """

    def __init__(self, l):
        self.label = l
        self.control_logic = None


######################################################################

class Network(object):
    """
    Holds regulatory system (genes, state space, stable states)
    """

    def __init__(self):
        self.genes = OrderedDict()
        self.state_spaces = {}
        self.stable_states = []
        self.state_space_dims = None
        self.up_to_date = False

    def add_gene(self, g):
        assert g.label not in self.genes
        self.genes[g.label] = g
        self.up_to_date = False

    def remove_gene(self, label):
        assert label in self.genes
        del self.genes[label]
        self.up_to_date = False

    def update(self):
        self.calculate_state_spaces()
        self.calculate_stable_states()
        self.up_to_date = True

    def calculate_state_spaces(self):
        self.state_spaces.clear()
        # calculate dimensions of state space
        dims = np.array(
            [2 for gene in self.genes.values()], dtype=int)
        self.state_space_dims = dims
        num_genes = len(dims)
        state_dict = OrderedDict((gene,0) for gene in self.genes.keys())
        # calculate state space for each gene
        for gene in self.genes.values():
            state_space = StateSpace(dims, self, gene.label)
            # calculate output for each coordinate and fill array
            for input_state, x in np.ndenumerate(state_space.array):
                # create dictionary of gene label : input state
                for i in range(num_genes):
                    sdk = list(state_dict.keys())
                    state_dict[sdk[i]] = input_state[i]
                state_space.array[input_state] = gene.control_logic(state_dict)
            state_space.filled = True
            self.state_spaces[gene.label] = state_space

    def calculate_stable_states(self):
        self.stable_states = []
        state_space = StateSpace(self.state_space_dims, self)
        for input_state, x in np.ndenumerate(state_space.array):
            stable = True
            # iterate through all genes and check if states match
            for i in range(len(self.state_space_dims)):
                ssi = list(state_space.inputs)
                gene_state_space = self.state_spaces[ssi[i]]
                assert gene_state_space.filled
                if gene_state_space.array[input_state] != input_state[i]:
                    stable = False
            if stable == True:
                self.stable_states += [input_state]

    def print_stable_states(self):
        num_stable_states = len(self.stable_states)
        if num_stable_states == 0:
            print("There are no stable states")
            return
        for i in range(len(self.stable_states)):
            print("Stable state %i" %(i+1))
            stable_state = self.stable_states[i]
            for j in range(len(stable_state)):
                sgk = self.genes.keys()
                print(sgk[j] + '\t%i' % stable_state[j])

    def transcription(self, gene, cell_state):
        return self.state_spaces[gene].array[cell_state]
        

#######################################################################

class StateSpace(object):

    def __init__(self, dims, net, l=""):
        self.array = np.zeros(dims, dtype=int)
        self.network = net
        self.label = l
        self.inputs = self.network.genes.keys()
        self.filled = False
        

#######################################################################

class Cell(object):
    """
    Holds information about cell state (RNA, protein, delays)
    """

    def __init__(self, t, i):
        self.tissue = t
        self.position = i
        self.genome = self.tissue.genome
        self.protein_state = np.zeros(len(self.genome), dtype=int)
        self.rna_state = np.zeros(len(self.genome), dtype=int)
        self.synthesis_delays = np.zeros(len(self.genome), dtype=int)
        self.decay_delays = np.zeros(len(self.genome), dtype=int)
        self.rna_ages = np.zeros(len(self.genome), dtype=int)
        self.protein_ages = np.zeros(len(self.genome), dtype=int)

    def update_state(self, transcription_states):
        for i in range(len(self.genome)):
            transcription = transcription_states[i]
            protein = self.protein_state[i]
            protein_age = self.protein_ages[i]
            decay_delay = self.decay_delays[i]
            synth_delay = self.synthesis_delays[i]
            rna_age = self.rna_ages[i]
            if not transcription:
                self.rna_ages[i] = 0
                if protein:
                    if protein_age >= decay_delay-1:
                        self.protein_state[i] = 0
                        self.protein_ages[i] = 0
                    else:
                        self.protein_ages[i] += 1
            else:
                if not protein:
                    if rna_age >= synth_delay-1:
                        self.protein_state[i] = 1
                    else:
                        self.rna_ages[i] += 1
                else:
                    self.rna_ages[i] = synth_delay
                    self.protein_ages[i] = 0
            self.rna_state[i] = transcription


######################################################################

class Tissue(object):
    """
    Container for a 1D array of cells.
    """

    def __init__(self, l, n, net):
        self.label = l
        self.num_cells = n
        self.network = net
        self.genome = list(self.network.genes.keys())
        self.cells = []
        for i in range(n):
            self.cells.append(Cell(self, i))
        self.simulations = OrderedDict()
        self.timepoint = 0
        self.synthesis_delay = 4
        self.decay_delay = 4
        self.colors = {}


    def set_protein_state(self, gene, states):
        # states should be i length array/list
        assert len(states) == self.num_cells
        assert gene in self.genome
        gene_index = self.genome.index(gene)
        for i in range(self.num_cells):
            self.cells[i].protein_state[gene_index] = states[i]

    def set_rna_state(self, gene, states):
        # states should be i length array/list
        assert len(states) == self.num_cells
        assert gene in  self.genome
        gene_index = self.genome.index(gene)
        for i in range(self.num_cells):
            self.cells[i].rna_state[gene_index] = states[i]

    def set_rna_ages(self, gene, states):
        # states should be i length array/list
        assert len(states) == self.num_cells
        assert gene in  self.genome
        gene_index = self.genome.index(gene)
        for i in range(self.num_cells):
            self.cells[i].rna_ages[gene_index] = states[i]

    def set_protein_ages(self, gene, states):
        # states should be i length array/list
        assert len(states) == self.num_cells
        assert gene in  self.genome
        gene_index = self.genome.index(gene)
        for i in range(self.num_cells):
            self.cells[i].protein_ages[gene_index] = states[i]

    def set_uniform_delays(self, sd, dd):
        self.synthesis_delay = sd
        self.decay_delay = dd
        for i in range(self.num_cells):
            for j in range(len(self.genome)):
                self.cells[i].synthesis_delays[j] = self.synthesis_delay
                self.cells[i].decay_delays[j] = self.decay_delay

    def set_delays(self, gene, sd, dd):
        gene_index = self.genome.index(gene)
        for i in range(self.num_cells):
            self.cells[i].synthesis_delays[gene_index] = sd
            self.cells[i].decay_delays[gene_index] = dd

    def simulate(self, l, d):
        assert l not in self.simulations
        self.timepoint = 0
        sim = Simulation(self, l, d)
        sim.record_state()
        while self.timepoint < d:
            self.advance_timepoint()
            sim.record_state()
        self.simulations[l] = sim

    def advance_timepoint(self):
        self.timepoint += 1
        num_cells = self.num_cells
        num_genes = len(self.genome)
        for i in range(num_cells):
            cell = self.cells[i]
            curr_protein = cell.protein_state
            transcription_state = np.zeros(num_genes, dtype=int)
            for j in range(num_genes):
                gene = self.genome[j]
                transcription_state[j] = self.network.transcription(gene, tuple(curr_protein))
            cell.update_state(transcription_state)

    def animate(self, l, framedelay=100):
        diagram = Visualisation(self, l)
        diagram.draw_diagram(framedelay)

    def set_colors(self, colors):
        self.colors = colors

        
####################################################################

class Simulation(object):

    """
    Stores data for a given simulation.
    """

    def __init__(self, tiss, l, d):
        self.tissue = tiss
        self.num_cells = self.tissue.num_cells
        self.num_genes = len(self.tissue.genome)
        self.gene_labels = self.tissue.genome
        self.label = l
        self.duration = d
        self.rna_history = np.zeros(
            (self.duration+1,self.num_cells, self.num_genes),dtype=int)
        self.protein_history = np.zeros(
            (self.duration+1, self.num_cells, self.num_genes),dtype=int)
        self.protein_age_history = np.zeros(
            (self.duration+1, self.num_cells, self.num_genes),dtype=int)
        self.rna_age_history = np.zeros(
            (self.duration+1, self.num_cells, self.num_genes),dtype=int)
        self.synthesis_delay_history = np.zeros(
            (self.duration+1, self.num_cells, self.num_genes),dtype=int)
        self.decay_delay_history = np.zeros(
            (self.duration+1, self.num_cells, self.num_genes),dtype=int)
        self.timepoint = 0

    def record_state(self):
        t = self.tissue.timepoint
        for i in range(self.num_cells):
            cell = self.tissue.cells[i]
            self.rna_history[t, i, :] = cell.rna_state
            self.protein_history[t, i, :] = cell.protein_state
            self.rna_age_history[t, i, :] = cell.rna_ages
            self.protein_age_history[t, i, :] = cell.protein_ages
            self.synthesis_delay_history[t, i, :] = cell.synthesis_delays
            self.decay_delay_history[t, i, :] = cell.decay_delays

            
####################################################################

class Visualisation(object):
    """
    Draws and animates simulation data
    """
    
    def __init__(self, t, l):
        self.tissue = t
        self.label = l
        self.font = "sans-serif"
        self.num_cells = self.tissue.num_cells
        self.num_genes = len(self.tissue.genome)
        self.genes = self.tissue.genome
        self.rna_patches = None
        self.protein_patches = None
        self.colors = OrderedDict()
        for gene in self.tissue.genome:
            self.colors[gene] = 'b'
        self.colors.update(self.tissue.colors)
        self.timestamp = None

    def draw_diagram(self, framedelay):
        cols = self.num_cells
        rows = self.num_genes
        # set up figure canvas
        fig = plt.figure(figsize=(cols + 1.25, rows + 1.5))
        ax = plt.axes([0,0,1,1])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-1.25, cols + .5)
        ax.set_ylim(-.5, rows + 1.5)
        fig.suptitle(self.label, family=self.font, size=24)
        # set up frame
        border = mpatches.Rectangle((0,0),
                                    cols, rows + .5, fill=False, ec='black', lw=2)
        ax.add_patch(border)
        lines = [mlines.Line2D([0,cols],[rows, rows], lw=2, color='black')]
        for i in range(1,cols):
            line = mlines.Line2D([i,i], [0, rows+.5], lw=2, color='black')
            lines.append(line)
        for line in lines:
            ax.add_line(line)
        del line
        # Set up label text
        for i in range(cols):
            plt.text(i+.25, rows+.15, "C "+str(i+1), 
                     family=self.font, size=20)
        for i in range(rows):
            gene = self.genes[i]
            plt.text(-.9, rows-i-.6, gene, 
                      family=self.font, size=20)
        self.timestamp = plt.text(
            self.num_cells-2.2, self.num_genes+.7, "Timepoint: ", family = self.font, size=20)           
        # Set up RNA / protein patches
        rna_patches = []
        protein_patches = []
        for i in range(cols):
            cell_rna = []
            cell_protein = []
            for j in range(rows):
                face_color = self.colors[self.tissue.genome[j]]
                rna_patch = mpatches.Rectangle((i+.05, rows-j-.5), .9, .4, 
                                               fill=False,
                                               fc=face_color,
                                               ec=None, lw=0)
                protein_patch = mpatches.Rectangle((i+.05, rows-j-.95), .9, .4, 
                                                   fill=False,
                                                   fc=face_color,
                                                   ec=None, lw=0)
                cell_rna.append(rna_patch)
                cell_protein.append(protein_patch)
                ax.add_patch(rna_patch)
                ax.add_patch(protein_patch)
            rna_patches.append(cell_rna)
            protein_patches.append(cell_protein)
        self.rna_patches = rna_patches
        self.protein_patches = protein_patches
        # Animate
        #anim = animation.FuncAnimation(fig, self.animation_function, 
        #                               init_func = self.init_function,
        #                               frames=self.tissue.timepoint+1, 
        #                               interval=framedelay, blit=True, repeat=False)

        #plt.show()
        self.init_function()
        for timepoint in range(self.tissue.timepoint + 1):    
            self.animation_function(timepoint)
            display(fig)
            clear_output(wait = True)
            input('Press key...')

    def init_function(self):
        for j in range(self.num_cells):
            for k in range(self.num_genes):
                self.rna_patches[j][k].set_fill(False)
                self.protein_patches[j][k].set_fill(False)
        self.timestamp.set_text("Timepoint: ")
        l = self.rna_patches + self.protein_patches
        return [item for sublist in l for item in sublist] + [self.timestamp]
    
    def animation_function(self, i):
        sim = self.tissue.simulations[self.label]
        for j in range(self.num_cells):
            for k in range(self.num_genes):
                # Set rna patches
                state_r = int(sim.rna_history[i][j][k])
                fill_state_r = bool(state_r)
                self.rna_patches[j][k].set_fill(fill_state_r)
                self.rna_patches[j][k].set_alpha(0.4)
                # Set protein patches
                state_p = int(sim.protein_history[i][j][k])
                fill_state_p = bool(state_p)
                self.protein_patches[j][k].set_fill(fill_state_p)
        l = self.rna_patches + self.protein_patches
        self.timestamp.set_text("Timepoint: "+str(i))
        return [item for sublist in l for item in sublist] + [self.timestamp]
