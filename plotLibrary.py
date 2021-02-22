'''
@Date: 2019-05-24 12:06:13
LastEditors: cany
LastEditTime: 2021-02-19 17:22:45
'''

import matplotlib.pyplot as plt
import pandas as pd
import logomaker, os


def main():

    df = pd.read_csv(csvfile)
    ww_df = df.drop(['pos'], axis=1)
    ww_logo = logomaker.Logo(ww_df,
                             font_name='Arial',
                             color_scheme='NajafabadiEtAl2017',
                             vpad=.05,
                             width=.8)

    # style using Axes methods

    ww_logo.style_spines(spines=['left', 'right'], visible=False)
    ww_logo.ax.set_xticks(range(len(ww_df)))

    ### add xlabel
    ww_logo.ax.set_xticklabels(df['pos'].tolist())
    ww_logo.ax.set_yticks([])
    ww_logo.ax.set_title(
        '{}'.format(os.path.basename(csvfile).rsplit('.', 1)[0]),
        {'fontsize': 18})
    plt.savefig('{}.png'.format(os.path.basename(csvfile).rsplit('.', 1)[0]),
                dpi=600)
    ww_logo.fig.tight_layout()


if __name__ == "__main__":
    import sys
    csvfile = sys.argv[1]
    main()
