from __future__ import division, print_function

import os
import numpy as np
import pandas as pd
import Tkinter as tk
import matplotlib 
from functools import partial
matplotlib.use("TkAgg")

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import hugs_pipe as hp
import lsst.daf.persistence

viz_dir = os.path.join(hp.io, 'viz-inspect-results')

class GUI(object):

    def __init__(self, master, cat_fn, out_fn, group_id, apply_cuts, use_ds9):

        self.master = master
        top_fr = tk.Frame(master)
        mid_fr = tk.Frame(master)
        bot_fr = tk.Frame(master)

        # status info
        self.status = tk.Label(
            master, text='', bd=1, relief='sunken', anchor='w')
        self.status.pack(side='bottom', fill='x', anchor='w')

        top_fr.pack(side='top', expand=0)
        bot_fr.pack(side='bottom', expand=0, pady=10)
        mid_fr.pack(side='bottom', expand=0)

        self.fig = Figure(figsize=(8, 7))
        self.ax = self.fig.add_subplot(111)

        # get catalog and remove duplicate entries
        self.cat = pd.read_csv(cat_fn)
        if apply_cuts:
            hp.cattools.cutter(self.cat, group_id=group_id, inplace=True)

        # sort by tract and patch
        self.cat.sort_values(['tract', 'patch'], inplace=True)
        self.cat.reset_index(drop=True, inplace=True)
        self.cat['is_candy'] = -1

        # initialize attributes
        self.current_idx = 0
        self.rgb_kw = dict(Q=8, dataRange=0.25)
        self.tract_prev = None
        self.patch_prev = None
        self.patch_cat = None
        self.is_new_idx = True
        self.is_new_patch = True
        self.ell_visible = False
        self.out_fn = out_fn
        self._coord = None
        self._viewer = hp.Viewer()

        # create bottom buttons
        tk.Label(bot_fr, text='Source Index').grid(
            row=0, column=0, sticky='w')
        self.idx_evar= tk.IntVar()
        self.idx_evar.set(self.current_idx)
        idx_entry = tk.Entry(bot_fr, textvariable=self.idx_evar)
        idx_entry.grid(row=0, column=1, sticky='w', padx=8)
        idx_entry.bind('<Return>', self.set_idx)

        prev_button = tk.Button(
            bot_fr, text='Previous', command=self.prev_idx)
        prev_button.grid(row=0, column=4, padx=12, sticky='w', columnspan=5)

        next_button = tk.Button(
            bot_fr, text='Next', command=self.next_idx)
        next_button.grid(row=0, column=9, sticky='w', columnspan=5)

        # create middle buttons
        padx = 15
        up_flag = partial(self.is_candy, 1)
        fn = os.path.join(viz_dir, 'buttons/up.gif')
        sig_img = tk.PhotoImage(file=fn)
        signal_button = tk.Button(
            mid_fr, image=sig_img, width='100', 
            height='100', command=up_flag) 
        signal_button.image = sig_img
        signal_button.grid(row=0, column=0, sticky='w', padx=padx)

        down_flag = partial(self.is_candy, 0)
        fn = os.path.join(viz_dir, 'buttons/down.gif')
        noise_img = tk.PhotoImage(file=fn)
        noise_button = tk.Button(
            mid_fr, image=noise_img, width='100', 
            height='100', command=down_flag)
        noise_button.image = noise_img
        noise_button.grid(row=0, column=1, sticky='w', padx=padx)

        question_flag = partial(self.is_candy, 2)
        fn = os.path.join(viz_dir, 'buttons//question-mark.gif')
        question_img = tk.PhotoImage(file=fn)
        question_button = tk.Button(
            mid_fr, image=question_img, width='100', 
            height='100', command=question_flag)
        question_button.image = question_img 
        question_button.grid(row=0, column=2, sticky='w', padx=padx)

        # create top buttons
        self.use_ds9 = tk.IntVar()
        ds9_button = tk.Checkbutton(
            top_fr, text='Show ds9', variable=self.use_ds9)
        self.use_ds9.set(use_ds9)
        ds9_button.pack(side='left', padx=padx)

        save_button = tk.Button(
            top_fr, text='Save Progess', command=self.save_progress)
        save_button.pack(side='left', padx=padx)

        source_button = tk.Button(
            top_fr, text='Toggle Sources', command=self.toggle_source)
        source_button.pack(side='left', padx=padx)

        quit_button = tk.Button(top_fr, text='Quit', command=self.quit)
        quit_button.pack(side='left', padx=padx)

        # useful key bindings
        master.bind('y', up_flag)
        master.bind('n', down_flag)
        master.bind('q', question_flag)
        master.bind('t', self.toggle_source)
        master.bind('s', self.save_progress)
        master.bind('<Left>', self.prev_idx)
        master.bind('<Right>', self.next_idx)

        # build canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.display_images()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        self.canvas.show()

    @property
    def viewer(self):
        return self._viewer

    @property
    def coord(self):
        if self.is_new_idx:
            c = ['x_hsc', 'y_hsc', 'ra', 'dec', 
                 'r_circ_seg', 'mu_2_i', 'is_candy']
            x, y, ra, dec, r_circ, mu, flag = self.cat.ix[self.current_idx, c]
            self._coord = (x, y)
            self.is_new_idx = False
            self.update_candy((ra, dec), r_circ, mu, flag)
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

    def next_idx(self, event=None):
        self.current_idx += 1
        self.is_new_idx = True
        self.idx_evar.set(self.current_idx)
        self.display_images()

    def prev_idx(self, event=None):
        if self.current_idx >0:
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

        if self.use_ds9.get():
            if self.is_new_patch:
                self.viewer.ds9_display_patch(frame=1)
            self.viewer.ds9_display_cutout(self.coord, frame=2)
            self.viewer.ds9_draw_ell(self.patch_cat, frame=2)

        self.ax.cla()
        self.ax.set(xticks=[], yticks=[])
        self.viewer.mpl_display_cutout(self.coord, rgb_kw=self.rgb_kw,
                                       subplots=(self.fig, self.ax))
        self.viewer.mpl_draw_ell(self.patch_cat, coord=self.coord)
        for p in self.viewer.current_axis.ax.patches:
            p.set_visible(self.ell_visible)
        for p in self.viewer.current_axis.points:
            p.set_visible(self.ell_visible)
        self.canvas.draw()

    def toggle_source(self, event=None):
        self.ell_visible = not self.ell_visible
        for p in self.viewer.current_axis.ax.patches:
            p.set_visible(self.ell_visible)
        for p in self.viewer.current_axis.points:
            p.set_visible(self.ell_visible)
        self.canvas.draw()

    def update_candy(self, coord, r_circ, mu, flag):
        txt = 'tract: {}   -    patch: {}   -   coord: {}, {}   -   '
        txt += 'r: {}   -   mu: {}   -   flag: {}'
        txt = txt.format(
            self.tract, self.patch, round(coord[0], 4), round(coord[1], 4), 
            round(r_circ, 2), round(mu, 2), flag)
        self.status.configure(text=txt)

    def is_candy(self, flag, event=None):
        self.cat.loc[self.current_idx, 'is_candy'] = flag
        print(self.cat.loc[self.current_idx, 'is_candy'])
        self.next_idx()

    def save_progress(self, event=None):
        self.cat.to_csv(self.out_fn, index=False)

    def quit(self):
        self.save_progress()
        self.master.destroy()


if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser('view candidates')

    parser.add_argument('group_id', type=str, help='group id')
    parser.add_argument('-o', '--out_fn', default=None)
    parser.add_argument('-c', '--cat_fn', default=None)
    parser.add_argument(
        '--no-cuts', dest='apply_cuts', action='store_false', default=True)
    parser.add_argument(
        '--no-ds9', dest='use_ds9', action='store_false', default=True)

    args = parser.parse_args()

    if args.out_fn is None:
        outdir = viz_dir
        args.out_fn = 'viz-inspect-group-'+args.group_id+'.csv'
        args.out_fn = os.path.join(outdir, args.out_fn)
    
    if args.cat_fn is None:
        catdir = os.path.join(hp.io, 'run-results')
        catdir = os.path.join(catdir, 'group_'+args.group_id)
        args.cat_fn = os.path.join(catdir, 'hugs-pipe-cat.csv')

    root = tk.Toplevel()
    root.title('hugs-pipe viz inspect')
    gui = GUI(root, args.cat_fn, args.out_fn, args.group_id, 
              args.apply_cuts, args.use_ds9)
    root.mainloop()
