from __future__ import division, print_function

import os
import numpy as np
import pandas as pd
import Tkinter as tk
import matplotlib 
matplotlib.use("TkAgg")

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import hugs_pipe as hp
import lsst.daf.persistence

class GUI(object):

    def __init__(self, master, cat_fn):

        frame = tk.Frame(master)
        self.fig = Figure(figsize=(12,8))
        self.ax = self.fig.add_subplot(111)

        # get catalog and remove duplicate entries
        self.cat = pd.read_csv(cat_fn)
        hp.cattools.remove_duplicates(self.cat)

        # sort by tract and patch
        self.cat.sort_values(['tract', 'patch'], inplace=True)
        self.cat.reset_index(drop=True, inplace=True)

        # initialize attributes
        self.current_idx = 0
        self.rgb_kw = dict(Q=8, dataRange=0.25)
        self.subplots_kw = dict(figsize=(8,8))
        self.tract_prev = None
        self.patch_prev = None
        self.patch_cat = None
        self.is_new_idx = True
        self.is_new_patch = True
        self.patch_visible = False
        self._coord = None
        self._viewer = hp.Viewer()

        # create entry
        tk.Label(frame, text='Source Index:').pack()
        self.idx_evar= tk.IntVar()
        self.idx_evar.set(self.current_idx)
        idx_entry = tk.Entry(frame, textvariable=self.idx_evar)
        idx_entry.pack()
        idx_entry.bind('<Return>', self.set_idx)

        # create buttons
        kws = dict(height=2, width=10)
        self.prev_button = tk.Button(
            frame, text='Previous', command=self.prev_idx, **kws)
        self.prev_button.pack(side='left')

        self.next_button = tk.Button(
            frame, text='Next', command=self.next_idx, **kws)
        self.next_button.pack(side='left', padx=10)

        self.source_button = tk.Button(
            frame, text='Toggle Sources', command=self.toggle_source, **kws)
        self.source_button.pack(side='top', padx=20)

        # build canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.display_images()
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        frame.pack()

    @property
    def viewer(self):
        return self._viewer

    @property
    def coord(self):
        if self.is_new_idx:
            cols = ['x_hsc', 'y_hsc']
            x, y = self.cat.ix[self.current_idx, cols]
            self._coord = (x, y)
            self.is_new_idx = False
        return self._coord

    def update_data_id(self):
        cols = ['tract', 'patch']
        self.tract, self.patch = self.cat.ix[self.current_idx, cols]
        if (self.tract!=self.tract_prev) or (self.patch!=self.patch_prev):
            self.data_id = dict(tract=self.tract, 
                                patch=self.patch, 
                                filter='HSC-I')
            self._viewer.set_data_id(self.data_id)
            self.tract_prev = self.tract
            self.patch_prev = self.patch
            self.is_new_patch = True
            cut = self.cat['tract']==self.tract
            cut &= self.cat['patch']==self.patch
            self.patch_cat = self.cat[cut]
        else:
            self.is_new_patch = False

    def next_idx(self):
        self.prev_idx = self.current_idx
        self.current_idx += 1
        self.is_new_idx = True
        self.idx_evar.set(self.current_idx)
        self.display_images()

    def prev_idx(self):
        if self.current_idx >0:
            self.prev_idx = self.current_idx
            self.current_idx -= 1
            self.is_new_idx = True
            self.idx_evar.set(self.current_idx)
            self.display_images()

    def set_idx(self, event=None):
        self.current_idx = self.idx_evar.get()
        self.is_new_idx = True
        self.display_images()

    def display_images(self):
        self.update_data_id()

        if self.is_new_patch:
            self.viewer.ds9_display_patch(frame=1)
        self.viewer.ds9_display_cutout(self.coord, frame=2)
        self.viewer.ds9_draw_ell(self.patch_cat, frame=2)

        self.ax.cla()
        self.ax.set(xticks=[], yticks=[])
        self.viewer.mpl_display_cutout(self.coord, rgb_kw=self.rgb_kw,
                                       subplots=(self.fig, self.ax))
        self.viewer.mpl_draw_ell(self.patch_cat, coord=self.coord)
        self.canvas.draw()

    def toggle_source(self):
        for p in self.viewer.current_axis.ax.patches:
            p.set_visible(self.patch_visible)
        for p in self.viewer.current_axis.points:
            p.set_visible(self.patch_visible)
        self.patch_visible = not self.patch_visible
        self.canvas.draw()

if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser('view candidates')
    parser.add_argument('cat_fn', type=str, help='cat file name')
    args = parser.parse_args()

    root = tk.Tk()
    root.title('hugs-pipe viz inspect')
    gui = GUI(root, args.cat_fn)
    root.mainloop()
