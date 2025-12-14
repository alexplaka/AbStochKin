"""  Graphing for AbStochKin simulations. """

#  Copyright (c) 2024-2026, Alex Plakantonakis.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
from pathlib import Path
from typing import Literal, Self

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import plotly.graph_objects as go

from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


class Graph:
    """
    Graphing class for displaying the results of AbStochKin simulations.

    Notes
    -----
    To successfully use the LaTeX engine for rendering text on Linux,
    run the following command in a terminal: `sudo apt install cm-super`.
    """

    # First, set some global matplotlib settings
    plt.rcParams['figure.autolayout'] = True  # tight layout
    plt.rcParams["figure.facecolor"] = 'lightgray'
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.titlesize'] = 9
    plt.rcParams['axes.labelsize'] = 7
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    plt.rcParams["legend.fontsize"] = 7
    plt.rcParams["legend.framealpha"] = 0.65

    def __init__(self,
                 /,
                 nrows=1,
                 ncols=1,
                 figsize=(5, 5),
                 dpi=300,
                 *,
                 backend: Literal['matplotlib', 'plotly'] = 'matplotlib',
                 **kwargs):
        self.backend = backend

        if self.backend == 'matplotlib':
            self.fig, self.ax = plt.subplots(nrows=nrows, ncols=ncols,
                                             figsize=figsize, dpi=dpi,
                                             **kwargs)
        elif self.backend == 'plotly':
            self.fig = go.Figure()
            self.fig.update_layout(
                width=figsize[0] * dpi * 0.6,  # Convert figsize to pixels
                height=figsize[1] * dpi * 0.6,  # Convert figsize to pixels
            )
            self.fig2 = None  # optional second figure
        else:
            error_msg = f"Unknown backend: {self.backend}. " \
                        f"Please choose from 'matplotlib' (default), 'plotly'."
            logger.error(error_msg)
            raise ValueError(error_msg)

    def setup_spines_ticks(self, ax_loc):
        """
        Set up the spines and ticks in a `matplotlib` graph.
        Make only the left and bottom spines/axes visible on the graph
        and place major ticks on them. Also set the minor ticks.
        """
        if self.backend == 'matplotlib':
            axs = self.ax if len(ax_loc) == 0 else self.ax[ax_loc]
            axs.spines[['left']].set_position('zero')
            axs.spines[['top', 'right']].set_visible(False)
            axs.xaxis.set_ticks_position('bottom')
            axs.yaxis.set_ticks_position('left')
            axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())

        else:  # self.backend == 'plotly'
            self.fig.update_layout(
                xaxis=dict(zeroline=True, showline=True, linewidth=2,
                           zerolinecolor="black", showticklabels=True,
                           tickmode='auto', ticks='outside'),
                yaxis=dict(zeroline=True, showline=True, linewidth=2,
                           zerolinecolor="black", showticklabels=True,
                           tickmode='auto', ticks='outside')
            )

    def plot_ODEs(self,
                  de_data,
                  *,
                  num_pts: int = 1000,
                  species: list[str] | tuple[str] | set[str] | None = None,
                  show_plot: bool = True,
                  ax_loc: tuple = ()) -> Self:
        """
        Plot the deterministic trajectories of all species obtained
        by obtaining the solution to a system of ODEs.

        Parameters
        ----------
        de_data : DEcalcs object
                 Data structure containing all the data related to
                 solving the system of ODEs.

        num_pts : int, default: 1000, optional
                 Number of points used to calculate DE curves at.
                 Used to approximate a smooth/continuous curve.

        species : sequence of strings, default: None, optional
                 An iterable sequence of strings specifying the species
                 names to plot. If no species are specified (the default),
                 then all species trajectories are plotted.

        show_plot : bool, default: True, optional
                 If True, show the plot.

        ax_loc : tuple, optional
                If the figure is made up of subplots, specify the location
                of the axis to draw the data at.
                Ex: for two subplots, the possible values of `ax_loc`
                are (0, ) and (1, ). That's because the `self.ax` object is
                a 1-D array. For figures with multiple rows and columns of
                subplots, a 2-D tuple is needed.
        """
        species = list(de_data.odes.keys()) if not species else species
        # t, y = ode_sol.t, ode_sol.y.T  # values at precomputed time pts
        t = np.linspace(de_data.odes_sol.t[0], de_data.odes_sol.t[-1],
                        num_pts)  # time points for obtaining...
        y = de_data.odes_sol.sol(t).T  # an approximately continuous solution

        self.setup_spines_ticks(ax_loc)

        if self.backend == 'matplotlib':
            axs = self.ax if len(ax_loc) == 0 else self.ax[ax_loc]
            axs.set(xlim=(0, de_data.odes_sol.t[-1]))

            for i, sp in enumerate(list(de_data.odes.keys())):
                if sp in species:
                    axs.plot(t, y[:, i], label=f"${sp}_{{DE}}$",
                             linestyle='--', linewidth=0.75, alpha=0.75)

            # axs.set(title="Deterministic trajectories")
            axs.set(xlabel=f"$t$ ({de_data.time_unit})", ylabel="$N$")
            axs.legend(loc='upper right')
            self.fig.tight_layout()

            if show_plot:
                plt.show()

        else:  # self.backend == 'plotly'
            for i, sp in enumerate(list(de_data.odes.keys())):
                if sp in species:
                    self.fig.add_trace(
                        go.Scatter(
                            x=t.tolist(),
                            y=y[:, i].tolist(),
                            mode='lines',
                            name=f"${sp}_{{DE}}$",
                            line=dict(dash='dash', width=1)
                        )
                    )

            self.fig.update_layout(
                xaxis=dict(title=f"$t \\, ({de_data.time_unit})$",
                           range=[-0.01 * de_data.odes_sol.t[-1], de_data.odes_sol.t[-1]]),
                yaxis=dict(title="$N$"),
            )

            if show_plot:
                self.fig.show()

        return self

    def plot_trajectories(self,
                          time,
                          data,
                          *,
                          species: list[str] | tuple[str] | set[str] | None = None,
                          show_plot: bool = True,
                          ax_loc: tuple = ()
                          ) -> Self:
        """ Graph simulation time trajectories. """
        self.setup_spines_ticks(ax_loc)
        species = list(data.keys()) if not species else species

        if self.backend == 'matplotlib':
            axs = self.ax if len(ax_loc) == 0 else self.ax[ax_loc]
            axs.set(xlim=(0, time[-1]))

            for sp, sp_data in data.items():
                if sp in species:
                    trajs = sp_data['N'].T
                    for traj in trajs:
                        axs.plot(time, traj, linewidth=0.25)

            axs.set(title="ABK trajectories")
            axs.set(xlabel="$t$ (sec)", ylabel="$N$")
            # axs.legend(loc='best')

            if show_plot:
                plt.show()

        else:  # self.backend == 'plotly'
            for sp, sp_data in data.items():
                if sp in species:
                    trajs = sp_data['N'].T
                    for i, traj in enumerate(trajs):
                        self.fig.add_trace(
                            go.Scatter(
                                x=time.tolist(),
                                y=traj.tolist(),
                                mode='lines',
                                name=f"${sp} \\; Run \\, {i}$",
                                line=dict(width=0.5)
                            )
                        )
            self.fig.update_layout(
                xaxis=dict(title="$t \\; (sec)$",
                           range=[-0.01 * time[-1], time[-1]]),
                yaxis=dict(title="$N$")
            )

            if show_plot:
                self.fig.show()

        return self

    def plot_avg_std(self,
                     time,
                     data,
                     *,
                     species: list[str] | tuple[str] | set[str] | None = None,
                     show_plot: bool = True,
                     ax_loc: tuple = ()) -> Self:
        """
        Graph simulation average trajectories and
        1-standard-deviation envelopes.
        """
        self.setup_spines_ticks(ax_loc)
        species = list(data.keys()) if not species else species

        if self.backend == 'matplotlib':
            axs = self.ax if len(ax_loc) == 0 else self.ax[ax_loc]
            axs.set(xlim=(0, time[-1]))

            for sp, sp_data in data.items():
                if sp in species:
                    axs.plot(time, sp_data['N_avg'],
                             linewidth=1.5, label=f"$<{sp}>$", alpha=0.5)
                    axs.fill_between(time,
                                     sp_data['N_avg'] - sp_data['N_std'],
                                     sp_data['N_avg'] + sp_data['N_std'],
                                     alpha=0.5, linewidth=0)

            axs.set(xlabel="$t$ (sec)", ylabel="$N$")
            axs.legend(loc='upper right')
            self.fig.tight_layout()

            if show_plot:
                plt.show()

        else:  # self.backend == 'plotly'
            for sp, sp_data in data.items():
                if sp in species:
                    self.fig.add_trace(
                        go.Scatter(
                            x=time.tolist(),
                            y=sp_data['N_avg'].tolist(),
                            mode='lines',
                            name=f"$<{sp}>$",
                            line=dict(width=2)
                        )
                    )

                    self.fig.add_trace(
                        go.Scatter(
                            x=time.tolist(),
                            y=(sp_data['N_avg'] + sp_data['N_std']).tolist(),
                            mode='lines',
                            name=f"$<{sp}> + \\sigma$",
                            line=dict(width=0),
                            showlegend=False
                        )
                    )

                    self.fig.add_trace(
                        go.Scatter(
                            x=time.tolist(),
                            y=(sp_data['N_avg'] - sp_data['N_std']).tolist(),
                            mode='lines',
                            line=dict(width=0),
                            name=f"$<{sp}> - \\sigma$",
                            fill='tonexty',
                            showlegend=False
                        )
                    )

            self.fig.update_layout(
                xaxis=dict(title="$t \\; (sec)$",
                           range=[-0.01 * time[-1], time[-1]]),
                yaxis=dict(title="$N$")
            )

            if show_plot:
                self.fig.show()

        return self

    def plot_eta(self,
                 time,
                 data,
                 *,
                 species: list[str] | tuple[str] | set[str] | None = None,
                 show_plot: bool = True,
                 ax_loc: tuple = ()) -> Self:
        """ Graph the coefficient of variation. """
        self.setup_spines_ticks(ax_loc)
        species = list(data.keys()) if not species else species

        if self.backend == "matplotlib":
            axs = self.ax if len(ax_loc) == 0 else self.ax[ax_loc]
            axs.set(xlim=(0, time[-1]))

            for sp, sp_data in data.items():
                if sp in species:
                    axs.plot(time, sp_data['eta'], linewidth=1.5, label=f"${sp}$")
                    axs.plot(time, sp_data['eta_p'], linewidth=1, linestyle='--',
                             label=f"${sp}_{{Poisson}}$", color=(0.5, 0.5, 0.5))

            axs.set(title="Coefficient of Variation, $\\eta$")
            axs.set(xlabel="$t$ (sec)", ylabel="$\\eta$")
            axs.legend(loc='upper right')

            if show_plot:
                plt.show()

        else:  # self.backend == 'plotly'
            for sp, sp_data in data.items():
                if sp in species:
                    self.fig.add_trace(
                        go.Scatter(
                            x=time.tolist(),
                            y=sp_data['eta'].tolist(),
                            mode='lines',
                            name=f"${sp}$",
                            line=dict(width=2)
                        ))

                    self.fig.add_trace(
                        go.Scatter(
                            x=time.tolist(),
                            y=sp_data['eta_p'].tolist(),
                            mode='lines',
                            name=f"${sp}_{{Poisson}}$",
                            line=dict(width=2, dash="dash")
                        ))

            self.fig.update_layout(
                title="Coefficient of Variation",
                xaxis=dict(title="$t \\; (sec)$",
                           range=[-0.01 * time[-1], time[-1]]),
                yaxis=dict(title="$\\eta$")
            )

            if show_plot:
                self.fig.show()

        return self

    def plot_het_metrics(self,
                         time,
                         proc_str: tuple[str, str],
                         proc_data: dict,
                         *,
                         het_attr='k',
                         show_plot: bool = True,
                         ax_loc: tuple = ()
                         ) -> Self:
        """
        Graph species- and process-specific metrics of population heterogeneity.
        """
        self.setup_spines_ticks(ax_loc)
        title = f"${proc_str[0].split(';')[0].replace(' ,', '\\hspace{10pt} ,').replace('->', '\\rightarrow')}$"

        if self.backend == 'matplotlib':
            axs = self.ax if len(ax_loc) == 0 else self.ax[ax_loc]
            axs.set(xlim=(0, time[-1]))
            axs.set(xlim=(0, time[-1]),
                    ylim=(0, 1.5 * np.max(proc_data[f"<{het_attr}_avg>"] + proc_data[
                        f"<{het_attr}_std>"])))

            axs.plot(time, proc_data[f'<{het_attr}_avg>'],
                     linewidth=1.5, label=f"$<{het_attr}>$", alpha=0.5)
            axs.fill_between(time,
                             proc_data[f"<{het_attr}_avg>"] - proc_data[f"<{het_attr}_std>"],
                             proc_data[f"<{het_attr}_avg>"] + proc_data[f"<{het_attr}_std>"],
                             alpha=0.5, linewidth=0)
            axs.tick_params(axis='y', labelcolor='blue')

            axs.set(title=title + (f"$, {proc_str[1]}$" if proc_str[1] != "" else ""),
                    xlabel="$t$ (sec)")

            if het_attr == 'Km':
                axs.set_ylabel("$K_m$", color='blue')
            elif het_attr == 'K50':
                axs.set_ylabel("$K_{50}$", color='blue')
            else:
                axs.set_ylabel(f"${het_attr}$", color='blue')

            # Second y-axis
            axs2 = axs.twinx()
            axs2.spines[['bottom']].set_position('zero')  # x axis
            axs2.spines[['right']].set_position(('axes', 1))  # y axis
            axs2.spines[['top', 'left', 'bottom']].set_visible(False)
            # axs2.spines['right'].set_color('red')
            axs2.set(ylim=(0, 1))
            axs2.tick_params(axis='y', labelcolor='red')
            axs2.yaxis.set_ticks_position('right')
            axs2.set_yticks([i for i in np.arange(0, 1.1, 0.1)])
            axs2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs2.grid(which='major', axis='y', color='r',
                      linestyle='--', linewidth=0.25, alpha=0.25)

            axs2.plot(time, proc_data['psi_avg'],
                      linewidth=1.5, label='$<\\psi>$', color='red', alpha=0.5)
            axs2.fill_between(time,
                              proc_data['psi_avg'] - proc_data['psi_std'],
                              proc_data['psi_avg'] + proc_data['psi_std'],
                              color='red', alpha=0.5, linewidth=0)
            axs2.set_ylabel("$\\psi$", color='red')

            if show_plot:
                plt.show()

            return self

        else:  # self.backend == 'plotly'
            self.fig.add_trace(
                go.Scatter(
                    x=time.tolist(),
                    y=proc_data[f'<{het_attr}_avg>'].tolist(),
                    mode='lines',
                    name=f"$<{het_attr}>$",
                    line=dict(width=2, color='blue'),
                )
            )

            self.fig.add_trace(
                go.Scatter(
                    x=time.tolist(),
                    y=(proc_data[f"<{het_attr}_avg>"] + proc_data[f"<{het_attr}_std>"]).tolist(),
                    mode='lines',
                    name=f"$<{het_attr}> + \\sigma$",
                    line=dict(width=0),
                    showlegend=False
                )
            )

            self.fig.add_trace(
                go.Scatter(
                    x=time.tolist(),
                    y=(proc_data[f"<{het_attr}_avg>"] - proc_data[f"<{het_attr}_std>"]).tolist(),
                    mode='lines',
                    name=f"$<{het_attr}> - \\sigma$",
                    line=dict(width=0),
                    fill='tonexty',
                    showlegend=False
                )
            )

            self.fig.update_layout(
                title=title + (f"$, {proc_str[1]}$" if proc_str[1] != "" else ""),
                xaxis=dict(title="$t \\; (sec)$",
                           range=[-0.01 * time[-1], time[-1]]),
                yaxis=dict(title="$K_m$" if het_attr == "Km" else "$K_{50}$" if het_attr == "K50" else f"${het_attr}$",
                           color="blue",
                           range=[0, 1.5 * np.max(proc_data[f"<{het_attr}_avg>"] + proc_data[f"<{het_attr}_std>"])])
            )

            if show_plot:
                self.fig.show()

            # Plot psi on a separate figure
            self.fig2 = go.Figure()
            self.fig2.update_layout(
                width=self.fig.layout.width,
                height=self.fig.layout.height,
                xaxis=dict(zeroline=True, showline=True, linewidth=2,
                           zerolinecolor="black", showticklabels=True,
                           tickmode='auto', ticks='outside'),
                yaxis=dict(zeroline=True, showline=True, linewidth=2,
                           zerolinecolor="black", showticklabels=True,
                           tickmode='auto', ticks='outside')
            )

            self.fig2.add_trace(
                go.Scatter(
                    x=time.tolist(),
                    y=proc_data['psi_avg'],
                    mode='lines',
                    name="$<\\psi>$",
                    line=dict(width=2, color='red'),
                )
            )

            self.fig2.add_trace(
                go.Scatter(
                    x=time.tolist(),
                    y=proc_data['psi_avg'] + proc_data['psi_std'],
                    mode='lines',
                    name="$<\\psi> + \\sigma$",
                    line=dict(width=0),
                    showlegend=False
                )
            )

            self.fig2.add_trace(
                go.Scatter(
                    x=time.tolist(),
                    y=proc_data['psi_avg'] - proc_data['psi_std'],
                    mode='lines',
                    name="$<\\psi> + \\sigma$",
                    line=dict(width=0),
                    fill='tonexty',
                    showlegend=False
                )
            )

            self.fig2.update_layout(
                title=title + (f"$, {proc_str[1]}$" if proc_str[1] != "" else ""),
                xaxis=dict(title="$t \\; (sec)$",
                           range=[-0.01 * time[-1], time[-1]]),
                yaxis=dict(title="$\\psi$",
                           color="red",
                           range=[0, 1],
                           tick0=0,
                           dtick=0.1)
            )

            if show_plot:
                self.fig2.show()

            return self

    def savefig(self, filename: str, image_format: str = 'svg', **kwargs):
        """ Save the figure as a file. """
        graph_path = Path('.') / 'plots_output'
        graph_path.mkdir(exist_ok=True)
        graph_path_svg = graph_path / f"{filename}.{image_format}"
        if self.backend == 'matplotlib':
            self.fig.savefig(graph_path_svg, format=image_format, **kwargs)
            plt.close(self.fig)
        else:  # self.backend == 'plotly'
            self.fig.write_image(graph_path_svg, format=image_format, **kwargs)
