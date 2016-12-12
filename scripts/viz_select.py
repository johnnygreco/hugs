from __future__ import division, print_function

import os, sys
import pandas as pd
import Tkinter as tk
import tkMessageBox
from functools import partial

import matplotlib 
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import hugs_pipe as hp
import lsst.daf.persistence

viz_dir = os.path.join(hp.io, 'viz-inspect-results')

class GUI(object):

    def __init__(self, root, master, cat_fn, out_fn, group_id, apply_cuts):

        # initialize attributes
        self.root = root
        self.master = master
        self.save_delta_t = 5 # minutes
        self.current_idx = 0
        self.group_id = group_id
        self.rgb_kw = dict(Q=8, dataRange=0.3)
        self.tract = None
        self.patch = None
        self.tract_prev = None
        self.patch_prev = None
        self.patch_cat = None
        self.is_new_idx = True
        self.is_new_patch = True
        self.ell_visible = False
        self.out_fn = out_fn
        self._coord = None
        self._viewer = hp.Viewer()
        self.master.withdraw()

        # if output catalog exists, check if we want to reload progress
        if os.path.isfile(out_fn):
            msg = 'Output catalog exists. Want to start where you left off?'
            answer = tkMessageBox.askyesno('HSC-HUGs Message', msg)
            if answer:
                self.cat = pd.read_csv(out_fn)
                flagged = self.cat['candy']==-1
                if flagged.sum()==0:
                    msg = 'All sources have been classified.'
                    tkMessageBox.showinfo('HSC-HUGs Message', msg)
                    sys.exit('Exiting without changing anything...')
                self.current_idx = self.cat[self.cat['candy']==-1].index[0]
            else:
                verify = tkMessageBox.askyesno(
                    'Verify', 'Progress will be lost, okay?', default='no')
                if verify:
                    self._load_cat(cat_fn, apply_cuts)
                else:
                    self.root.destroy()
                    sys.exit('Exiting without changing anything...')
        else:
            self._load_cat(cat_fn, apply_cuts)

        top_fr = tk.Frame(self.master)
        mid_fr = tk.Frame(self.master)
        bot_fr = tk.Frame(self.master)
        
        # status info
        self.status = tk.Text(
            self.master, height=1, relief='sunken', bg=self.master.cget('bg'))
        self.status.pack(side='bottom', fill='x', anchor='w')

        top_fr.pack(side='top', expand=0)
        bot_fr.pack(side='bottom', expand=0, pady=10)
        mid_fr.pack(side='bottom', expand=0)

        self.fig = Figure(figsize=(9, 7))
        self.ax = self.fig.add_subplot(111)

        # create bottom buttons
        tk.Label(bot_fr, text='Source Index').grid(
            row=0, column=0, sticky='w')
        self.idx_evar= tk.IntVar()
        self.idx_evar.set(self.current_idx)
        idx_entry = tk.Entry(bot_fr, textvariable=self.idx_evar, takefocus=0)
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
        up_flag = partial(self.set_flag, 'candy')
        fn = os.path.join(viz_dir, 'buttons/up.gif')
        sig_img = tk.PhotoImage(file=fn)
        signal_button = tk.Button(
            mid_fr, image=sig_img, width='100', 
            height='100', command=up_flag) 
        signal_button.image = sig_img
        signal_button.grid(row=0, column=0, sticky='w', padx=padx)

        down_flag = partial(self.set_flag, 'junk')
        fn = os.path.join(viz_dir, 'buttons/down.gif')
        noise_img = tk.PhotoImage(file=fn)
        noise_button = tk.Button(
            mid_fr, image=noise_img, width='100', 
            height='100', command=down_flag)
        noise_button.image = noise_img
        noise_button.grid(row=0, column=1, sticky='w', padx=padx)

        question_flag = partial(self.set_flag, 'ambiguous')
        fn = os.path.join(viz_dir, 'buttons/question-mark.gif')
        question_img = tk.PhotoImage(file=fn)
        question_button = tk.Button(
            mid_fr, image=question_img, width='100', 
            height='100', command=question_flag)
        question_button.image = question_img 
        question_button.grid(row=0, column=2, sticky='w', padx=padx)

        # create top buttons
        ds9_button = tk.Button(
            top_fr, text='Display with ds9', command=self.show_ds9)
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
        self.master.bind('1', up_flag)
        self.master.bind('2', down_flag)
        self.master.bind('3', question_flag)
        self.master.bind('t', self.toggle_source)
        self.master.bind('s', self.save_progress)
        self.master.bind('<Down>', self.prev_idx)
        self.master.bind('<Up>', self.next_idx)
        self.master.bind_all('<1>', lambda event: event.widget.focus_set())

        # build canvas
        self.fig.set_tight_layout(True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.display_image()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        self.canvas.show()
        self.master.deiconify()

        # save progress every save_delta_t minutes
        self.root.after(self.save_delta_t*60*1000, self.save_progress)

    @property
    def viewer(self):
        return self._viewer

    @property
    def coord(self):
        if self.is_new_idx:
            x, y = self.cat.ix[self.current_idx, ['x_hsc', 'y_hsc']]
            self._coord = (x, y)
            self.is_new_idx = False
            self.update_info()
        return self._coord

    def _load_cat(self, cat_fn, apply_cuts):
        self.cat = pd.read_csv(cat_fn)
        if apply_cuts:
            hp.cattools.cutter(self.cat, group_id=self.group_id, inplace=True)
        # sort by tract and patch
        self.cat.sort_values(['tract', 'patch'], inplace=True)
        self.cat.reset_index(drop=True, inplace=True)
        self.cat['candy'] = -1
        self.cat['junk'] = -1
        self.cat['ambiguous'] = -1

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
            cut = self.cat['tract']==self.tract
            cut &= self.cat['patch']==self.patch
            self.patch_cat = self.cat[cut]
            self.is_new_patch = True

    def next_idx(self, event=None):
        self.current_idx += 1
        if self.current_idx > (len(self.cat)-1):
            msg = 'Congrats, visual inspection for group '
            msg += self.group_id+' is complete!'
            w = tk.Tk()
            w.withdraw()
            tkMessageBox.showinfo('HSC-HUGs Message', msg, parent=w)
            w.destroy()
            self.quit()
        else:
            self.is_new_idx = True
            self.idx_evar.set(self.current_idx)
            self.display_image()

    def prev_idx(self, event=None):
        if self.current_idx >0:
            self.current_idx -= 1
            self.is_new_idx = True
            self.idx_evar.set(self.current_idx)
            self.display_image()

    def set_idx(self, event=None):
        self.current_idx = self.idx_evar.get()
        self.is_new_idx = True
        self.display_image()

    def show_ds9(self):
        if self.is_new_patch:
            self.viewer.ds9_display_patch(frame=1)
            self.is_new_patch = False
        self.viewer.ds9_display_cutout(self.coord, frame=2)
        self.viewer.ds9_draw_ell(self.patch_cat, frame=2)

    def display_image(self):
        self.update_data_id()
        self.ax.cla()
        self.ax.set(xticks=[], yticks=[])
        self.viewer.mpl_display_cutout(self.coord, rgb_kw=self.rgb_kw,
                                       subplots=(self.fig, self.ax), 
                                       size=50)
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

    def update_info(self):
        txt = 'tract: {}  -   patch: {}  -  coord: {:.5f} {:.5f}  -  '
        txt += 'r: {:.2f}  -  mu: {:.2f}  -  flag: {}'
        cols = ['ra', 'dec', 'a_2_sig', 'mu_2_i']
        cols += ['candy', 'junk', 'ambiguous']
        info = self.cat.ix[self.current_idx, cols]
        ra, dec, size, mu = info[['ra', 'dec', 'a_2_sig', 'mu_2_i']]
        flags = info[['candy', 'junk', 'ambiguous']]
        flag = flags[flags==1]
        if len(flag)==1:
            flag = flag.index[0]
        else:
            flag = 'n/a'
        txt = txt.format(self.tract, self.patch, ra, dec, size, mu, flag)
        self.status.config(state='normal')
        self.status.delete(1.0, 'end')
        self.status.insert('insert', txt)
        self.status.config(state='disabled')
        #self.status.configure(text=txt)

    def set_flag(self, flag, event=None):
        self.cat.loc[self.current_idx, flag] = 1
        others = [f for f in ['candy', 'junk', 'ambiguous'] if f!=flag]
        self.cat.loc[self.current_idx, others] = 0
        self.next_idx()

    def save_progress(self, event=None):
        self.cat.to_csv(self.out_fn, index=False)
        self.root.after(self.save_delta_t*60*1000, self.save_progress)

    def quit(self):
        self.save_progress()
        self.root.destroy()


if __name__=='__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser('view candidates')
    parser.add_argument('group_id', type=str, help='group id')
    parser.add_argument('-o', '--out_fn', default=None)
    parser.add_argument('-c', '--cat_fn', default=None)
    parser.add_argument(
        '--no-cuts', dest='apply_cuts', action='store_false', default=True)
    args = parser.parse_args()

    if args.out_fn is None:
        outdir = viz_dir
        args.out_fn = 'viz-inspect-group-'+args.group_id+'.csv'
        args.out_fn = os.path.join(outdir, args.out_fn)
    
    if args.cat_fn is None:
        catdir = os.path.join(hp.io, 'run-results')
        catdir = os.path.join(catdir, 'group-'+args.group_id)
        args.cat_fn = os.path.join(catdir, 'hugs-pipe-cat.csv')

    root = tk.Tk()
    root.withdraw()
    top = tk.Toplevel(root)
    top.protocol("WM_DELETE_WINDOW", root.destroy)
    top.title('hugs-pipe viz inspect group '+args.group_id)
    gui = GUI(root, top, args.cat_fn, args.out_fn, args.group_id, 
              args.apply_cuts)
    root.mainloop()
