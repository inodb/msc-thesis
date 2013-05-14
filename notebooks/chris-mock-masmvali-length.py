# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <markdowncell>

# ##Determine length distribution of chris-mock metassemble assemblies
#
#
# Uses masmvali for contig length calculation.

# <codecell>

#!echo "`cd /glob/inod/github/masm-vali/masmvali/ && git log -1`"

# <codecell>

from numpy import *
from matplotlib.pyplot import *

import sys
sys.path.append('/glob/inod/github/masm-vali/masmvali/')
import assembly
reload(assembly)
import nucmer
reload(nucmer)
import refgenome
reload(refgenome)

# <markdowncell>

# Name to color, shapes and filename mappings for both even and unbalanced library. Take assembly with best L50 from each assembly strategy (manually determined).

# <codecell>

names = ['velvetnoscaf','velvetscaf','velvetnoscafminimus2','velvetnoscafnewbler','velvetnoscafbambus2','velvetnoscafminimus2bambus2','velvetnoscafnewblerbambus2',
         'metavelvetnoscaf','metavelvetscaf','metavelvetnoscafminimus2','metavelvetnoscafnewbler','metavelvetnoscafbambus2','metavelvetnoscafminimus2bambus2','metavelvetnoscafnewblerbambus2',
         'raynoscaf','rayscaf','raynoscafminimus2','raynoscafnewbler']

colord = dict(zip(names,
              ['k','k','k','k','k','k','k',
               'g','g','g','g','g','g','g',
               'r','r','r','r']
              )
          )

shapesd = dict(zip(names,
               ['-','-','--','-','-','--','-',
                '-','-','--','-','-','--','-',
                '-','-','--','-']
               )
           )

# Even
metassembledir = '/proj/b2010008/nobackup/projects/chris-mock/Project_ID793_dAmore/Sample_50ng_even/metassemble/'

filesd = {}
filesd["Sample_50ng_even"] = dict(zip(names,
              [metassembledir + 'assemblies/velvet/noscaf/noscaf_25/ma-contigs.fa',
               metassembledir + 'assemblies/velvet/scaf/scaf_27/ma-scaffolds.fa',
               metassembledir + 'assemblies/velvet/noscaf/minimus2/ma-merge.fa',
               metassembledir + 'assemblies/velvet/noscaf/newbler/ma-merge.fa',
               metassembledir + 'assemblies/velvet/noscaf/noscaf_31/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/velvet/noscaf/minimus2/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/velvet/noscaf/newbler/bambus2/bambus2.scaffold.linear.fasta',

               metassembledir + 'assemblies/metavelvet/noscaf/noscaf_25/ma-contigs.fa',
               metassembledir + 'assemblies/metavelvet/scaf/scaf_27/ma-scaffolds.fa',
               metassembledir + 'assemblies/metavelvet/noscaf/minimus2/ma-merge.fa',
               metassembledir + 'assemblies/metavelvet/noscaf/newbler/ma-merge.fa',
               metassembledir + 'assemblies/metavelvet/noscaf/noscaf_31/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/metavelvet/noscaf/minimus2/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/metavelvet/noscaf/newbler/bambus2/bambus2.scaffold.linear.fasta',

               metassembledir + 'assemblies/ray/noscaf/noscaf_21/ma-contigs.fa',
               metassembledir + 'assemblies/ray/scaf/scaf_21/ma-scaffolds.fa',
               metassembledir + 'assemblies/ray/noscaf/minimus2/ma-merge.fa',
               metassembledir + 'assemblies/ray/noscaf/newbler/ma-merge.fa']
              )
          )
coordsd = {}
coordsd["Sample_50ng_even"] = dict([(k, "/".join(v.split('/')[:-1]) + '/val/nucmer.coords') for k, v in filesd["Sample_50ng_even"].items()])

# Unbalanced
metassembledir = '/proj/b2010008/nobackup/projects/chris-mock/Project_ID793_dAmore/Sample_50ng_unbalanced/metassemble/'

filesd["Sample_50ng_unbalanced"] = dict(zip(names,
              [metassembledir + 'assemblies/velvet/noscaf/noscaf_27/ma-contigs.fa',
               metassembledir + 'assemblies/velvet/scaf/scaf_27/ma-scaffolds.fa',
               metassembledir + 'assemblies/velvet/noscaf/minimus2/ma-merge.fa',
               metassembledir + 'assemblies/velvet/noscaf/newbler/ma-merge.fa',
               metassembledir + 'assemblies/velvet/noscaf/noscaf_31/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/velvet/noscaf/minimus2/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/velvet/noscaf/newbler/bambus2/bambus2.scaffold.linear.fasta',

               metassembledir + 'assemblies/metavelvet/noscaf/noscaf_29/ma-contigs.fa',
               metassembledir + 'assemblies/metavelvet/scaf/scaf_29/ma-scaffolds.fa',
               metassembledir + 'assemblies/metavelvet/noscaf/minimus2/ma-merge.fa',
               metassembledir + 'assemblies/metavelvet/noscaf/newbler/ma-merge.fa',
               metassembledir + 'assemblies/metavelvet/noscaf/noscaf_31/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/metavelvet/noscaf/minimus2/bambus2/bambus2.scaffold.linear.fasta',
               metassembledir + 'assemblies/metavelvet/noscaf/newbler/bambus2/bambus2.scaffold.linear.fasta',

               metassembledir + 'assemblies/ray/noscaf/noscaf_21/ma-contigs.fa',
               metassembledir + 'assemblies/ray/scaf/scaf_21/ma-scaffolds.fa',
               metassembledir + 'assemblies/ray/noscaf/minimus2/ma-merge.fa',
               metassembledir + 'assemblies/ray/noscaf/newbler/ma-merge.fa']
              )
          )

coordsd["Sample_50ng_unbalanced"] = dict([(k, "/".join(v.split('/')[:-1]) + '/val/nucmer.coords') for k, v in filesd["Sample_50ng_unbalanced"].items()])

# <markdowncell>

# Name sets of names, so we can compare each set separately. Otherwise you don't see anything in the plot

# <codecell>

asm_sets = {}
asm_sets['contig'] = [n for n in names if n.endswith('noscaf')]
asm_sets['scaf'] = [n for n in names if n.endswith('scaf') and not n.endswith('noscaf')]
asm_sets['contig-merge'] = [n for n in names if n.endswith('newbler') or n.endswith('minimus2')]
asm_sets['contig-scaf'] = [n for n in names if n.endswith('noscafbambus2')]
asm_sets['contig-merge-scaf'] = [n for n in names if n.endswith('newblerbambus2') or n.endswith('minimus2bambus2')]
asm_sets

# <codecell>

METAGENOME_LENGTH=195024525

# <markdowncell>

# Plot length distribution function, uses masmvali

# <codecell>

def plot_length(assemblies, names, colors, linestyles, cut_off=100):
    ax = subplot(1,1,1)
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        clacs = np.cumsum(cla[cla > cut_off])

        # Plot contig length vs sum of bases
        plot(cla[cla > cut_off], clacs, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=2.0)

        # Add L50
        l50 = cla[cla > cut_off][clacs >= clacs[len(clacs)-1]/2][0]
        plot(l50, clacs[clacs >= clacs[len(clacs)-1]/2][0], color=colors[i], marker='o', markersize=10)
    ax.set_xscale('log')
    ax.set_ylim([0,METAGENOME_LENGTH])
    gca().invert_xaxis()
    # Labels
    xlabel('Contig/Scaffold length')
    ylabel('Sum of bases')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper left') # or lower right
    #ax.set_yscale('log')
    ax.set_xlim([10**6,10**2])

# <markdowncell>

# Compare all sets

# <codecell>

# Change figure size
rcParams['figure.figsize'] = 8, 6
rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 100
outputdir = '/glob/inod/github/msc-thesis/figures/notebooks/chris-mock-length-distribution'

# <codecell>

for k in asm_sets:
    print "Comparing: " + k
    plot_length([filesd["Sample_50ng_even"][name] for name in asm_sets[k]], asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdir + '/' + 'even-' + k + '.png')
    savefig(outputdir + '/' + 'even-' + k + '.pdf')
    show()
    close()

# <markdowncell>

# ###Same for unbalanced

# <codecell>

for k in asm_sets:
    print "Comparing: " + k
    plot_length([filesd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdir + '/' + 'uneven-' + k + '.png')
    savefig(outputdir + '/' + 'uneven-' + k + '.pdf')
    show()
    close()

# <markdowncell>

# #Assemblathon 2 style length distribution plot

# <markdowncell>

# Way to get value from an array where conditions are based on another array

# <codecell>

# all existing values
r = arange(1,101,1)
s = array([2,5,7,8])
for i in xrange(len(s)):
    s[i] = r[where(r == s[i])][0]
print s

# unexisting values
r = arange(1,101,1)
s = array([2,5,7,8,500])
for i in xrange(len(s)):
    v = r[where(r == s[i])]
    s[i] = v[0] if v else 0
print s

# <markdowncell>

# Function to plot assemblathon 2 style length distribtion

# <codecell>

def plot_length_asmbla2(assemblies, names, colors, linestyles, mg_ref_len=195000000):
    ax = subplot(1,1,1)
    lmg_prc = arange(0,101,1)
    lmg_len = arange(0,101,1)

    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        clcsa = cumsum(cla)
        for j in xrange(len(lmg_prc)):
            v = cla[where(clcsa >= (lmg_prc[j] * mg_ref_len / 100))]
            # if 0 is used instead of 1, the plot doesn't plot the point
            lmg_len[j] = v[0] if v.any() else 1
        plot(lmg_prc, lmg_len, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=2.0)
    # Labels
    xlabel('LMG(X) %')
    ylabel('Contig/Scaffold LMG(X)')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    ax.set_yscale('log')

# <markdowncell>

# ### 50 ng even

# <codecell>

#for k in asm_sets:
#    print "Comparing: " + k
#    plot_length_asmbla2([filesd["Sample_50ng_even"][name] for name in asm_sets[k]], asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
#    show()

# <markdowncell>

# ###50ng unbalanced

# <codecell>

#for k in asm_sets:
#    print "Comparing: " + k
#    plot_length_asmbla2([filesd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
#    show()

# <markdowncell>

# Now plot against non-overlapping coverage

# <codecell>

import sqlite3

# <codecell>

coordsd["Sample_50ng_unbalanced"]["velvetnoscaf"]

# <codecell>

refphylfile = "/glob/inod/github/masm-vali/masmvali/test/data/chris-mock/reference/phylogeny-references.tsv"
refstatsfiled = {}

refstatsfiled["Sample_50ng_even"] = "/proj/b2010008/nobackup/projects/chris-mock/Project_ID793_dAmore/Sample_50ng_even/metassemble/reference-stats/ref.stats"
refstatsfiled["Sample_50ng_unbalanced"] = "/proj/b2010008/nobackup/projects/chris-mock/Project_ID793_dAmore/Sample_50ng_unbalanced/metassemble/reference-stats/ref.stats"

refsets = {}
refsets["Sample_50ng_even"] = refgenome.ReferenceSet(refphylfile, refstatsfiled["Sample_50ng_even"])
refsets["Sample_50ng_unbalanced"] = refgenome.ReferenceSet(refphylfile, refstatsfiled["Sample_50ng_unbalanced"])

# <markdowncell>

# Plot contig length vs metagenome coverage. This function does one point for every contig, but goes extremely slow, so I decided to do only 7 points. See below.

# <codecell>

def plot_length_vs_metagenome_cov(assemblies, coordsd, refset, names, colors, linestyles):
    ax = subplot(1,1,1)
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        c_cova = np.copy(cla)
        c = nucmer.Coords(coordsd[i])
        for j in xrange(len(c_cova)):
            c.calc_genome_contig_cov_in_bases(refset, cut_off=int(cla[j]), count_contig_once=True)
            c_cova[j] = float(c.ref_sum_cov) / c.ref_sum_bases
        plot(cla, c_cova, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=2.0)
    gca().invert_xaxis()
    # Labels
    xlabel('Contig/Scaffold length')
    ylabel('Non-overlapping reference metagenome coverage')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    #ax.set_yscale('log')

# <codecell>

print refsets["Sample_50ng_unbalanced"]

# <markdowncell>

# Plot length vs non overlapping metagenome coverage of each contig. Take best alignment for each contig. 7 points per assembly.

# <codecell>

POINT_SHAPES = dict(zip(names,
              ['ko','kv','ko','kd','k^','k*','kD',
               'go','gv','go','gd','g^','g*','gD',
               'ro','rv','ro','rd']
              )
          )

def plot_length_vs_metagenome_cov2(assemblies, coordsd, refset, names, colors, linestyles, cut_off=100):

    #lengthsa = np.array([item for sublist in map(lambda x: range(10**x,10**(x+1),10**x),range(2,7)) for item in sublist])
    ax = subplot(1,1,1)
    for i, asm in enumerate(assemblies):

        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)

        lengthsa = sort([item for sublist in map(lambda x: (range(10**x, 10**(x+1), (10**(x+1)-10**x) /4)), range(2,6)) for item in sublist ] + [10**6] + [contigs.l50])

        # Determine coverage
        c_cova = np.copy(lengthsa)
        c = nucmer.Coords(coordsd[i])
        for j in xrange(len(c_cova)):
            c.calc_genome_contig_cov_in_bases(refset, cut_off=int(lengthsa[j]), count_contig_once=True)
            c_cova[j] = (c.ref_sum_cov * 100) / c.ref_sum_bases
        plot(lengthsa, c_cova, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=2.0)
        plot(contigs.l50, c_cova[lengthsa == contigs.l50], 'o', color=colors[i], markersize=10)

    ax.set_xscale('log')
    gca().invert_xaxis()
    # Labels
    ylabel('Metagenome coverage')
    xlabel('Contig/Scaffold length')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, numpoints=1, loc='upper left')

    ax.set_ylim([0,100])
    ax.set_xlim([10**6,0])

# <codecell>

# Change figure size
rcParams['figure.figsize'] = 8, 6
rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 100
outputdir = '/glob/inod/github/msc-thesis/figures/notebooks/chris-mock-length-distribution-cov'

for k in asm_sets:
    print "Comparing: " + k
    plot_length_vs_metagenome_cov2([filesd["Sample_50ng_even"][name] for name in asm_sets[k]], [coordsd["Sample_50ng_even"][name] for name in asm_sets[k]], refsets["Sample_50ng_even"].refs, asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdir + '/' + 'even-' + k + '.png')
    savefig(outputdir + '/' + 'even-' + k + '.pdf')
    show()
    close()

    plot_length_vs_metagenome_cov2([filesd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], [coordsd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], refsets["Sample_50ng_unbalanced"].refs, asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdir + '/' + 'uneven-' + k + '.png')
    savefig(outputdir + '/' + 'uneven-' + k + '.pdf')
    show()
    close()

# <markdowncell>

# #FRC style length distribution
# ### 50 ng unbalanced

# <codecell>

def plot_length_frc_style(assemblies, names, colors, linestyles, cut_off=100):
    ax = subplot(1,1,1)
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        nr_of_contigs = len(cla[cla > cut_off])
        # x contig lengths, 100 dots
        xarray = np.arange(nr_of_contigs, step=nr_of_contigs/100)
        # sum of bases
        yarray = np.cumsum(cla[cla > cut_off])[xarray]
        plot(xarray, yarray, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=2.0)
    # Labels
    xlabel('Contig n')
    ylabel('Sum of bases')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='lower right')
    ax.set_ylim([0,METAGENOME_LENGTH])
    ax.set_xlim(0,350000)

# <codecell>

def plot_length_frc_style_metagenome_cov(assemblies, coordsd, refset, names, colors, linestyles, cut_off=100, number_of_points_per_line=20):
    ax = subplot(1,1,1)
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        nr_of_contigs = len(cla[cla > cut_off])
        # x contig lengths
        xarray = np.arange(nr_of_contigs, step=nr_of_contigs/number_of_points_per_line)

        # metagenome cov
        yarray = np.copy(xarray)
        c = nucmer.Coords(coordsd[i])
        for j in xrange(len(xarray)):
            clength = cla[cla > cut_off][xarray[j]]
            c.calc_genome_contig_cov_in_bases(refset, cut_off=int(clength), count_contig_once=True)
            yarray[j] = (c.ref_sum_cov * 100) / c.ref_sum_bases
        plot(xarray, yarray, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=2.0)
    # Labels
    xlabel('Contig n')
    ylabel('Metagenome coverage')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='lower right')
    ax.set_ylim([0,100])
    ax.set_xlim(0,350000)

# <codecell>

# Change figure size
rcParams['figure.figsize'] = 8, 6
rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 100
outputdirlen = '/glob/inod/github/msc-thesis/figures/notebooks/chris-mock-length-distribution-roc'
outputdircov = '/glob/inod/github/msc-thesis/figures/notebooks/chris-mock-length-distribution-cov-roc'

# <markdowncell>

# ###even

# <codecell>

for k in asm_sets:
    print "Comparing: " + k
    plot_length_frc_style([filesd["Sample_50ng_even"][name] for name in asm_sets[k]], asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdirlen + '/' + 'even-' + k + '.png')
    savefig(outputdirlen + '/' + 'even-' + k + '.pdf')
    show()
    close()
    plot_length_frc_style_metagenome_cov([filesd["Sample_50ng_even"][name] for name in asm_sets[k]], [coordsd["Sample_50ng_even"][name] for name in asm_sets[k]], refsets["Sample_50ng_even"].refs, asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdircov + '/' + 'even-' + k + '.png')
    savefig(outputdircov + '/' + 'even-' + k + '.pdf')
    show()
    close()

# <markdowncell>

# ###uneven

# <codecell>

for k in asm_sets:
    print "Comparing: " + k
    plot_length_frc_style([filesd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdirlen + '/' + 'uneven-' + k + '.png')
    savefig(outputdirlen + '/' + 'uneven-' + k + '.pdf')
    show()
    close()
    plot_length_frc_style_metagenome_cov([filesd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], [coordsd["Sample_50ng_unbalanced"][name] for name in asm_sets[k]], refsets["Sample_50ng_unbalanced"].refs, asm_sets[k], [colord[name] for name in asm_sets[k]], [shapesd[name] for name in asm_sets[k]])
    savefig(outputdircov + '/' + 'uneven-' + k + '.png')
    savefig(outputdircov + '/' + 'uneven-' + k + '.pdf')
    show()
    close()

