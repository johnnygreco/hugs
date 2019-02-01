import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


__all__ = ['plot_sep_sources']


def plot_sep_sources(image, catalog, ec='lime', scale_ell=6, subplots=None,
                     mark_centers=False, subplot_kw=dict(figsize=(10, 10)),
                     ell_type='ab', per_lo=1.0, per_hi=99.0, mask=None,
                     mask_kws=dict(cmap='Blues_r', alpha=0.5)):

    fig, ax = subplots if subplots is not None else plt.subplots(**subplot_kw)

    if image is not None:
        if len(image.shape)==2:
            vmin, vmax = np.percentile(image, [per_lo, per_hi])
            ax.imshow(image, vmin=vmin, vmax=vmax,
                      origin='lower', cmap='gray_r')
        else:
            ax.imshow(image, origin='lower')
        ax.set_xlim(0, image.shape[1])
        ax.set_ylim(0, image.shape[0])

    ax.set(xticks=[], yticks=[])

    for src in catalog:
        if ell_type == 'ab':
            a = src['a']
            b = src['b']
        elif 'kronrad' in ell_type:
            kronrad = src[ell_type]
            a = kronrad * src['a']
            b = kronrad * src['b']
        e = Ellipse((src['x'], src['y']),
                     width=scale_ell*a,
                     height=scale_ell*b,
                     angle=src['theta']*180/np.pi, fc='none', ec=ec)
        ax.add_patch(e)

        if mark_centers:
            ax.plot(x_c, y_c, 'r+')

    if mask is not None:
        mask = mask.astype(float)
        mask[mask==0.0] = np.nan
        ax.imshow(mask, **mask_kws)

    return fig, ax


