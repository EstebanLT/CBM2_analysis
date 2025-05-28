# -*- coding: utf-8 -*-

"""
Created on Mon Oct 29 2018

@author: Esteban López Tavera
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.spatial import distance_matrix
from datetime import datetime
import os
import matplotlib.animation as animation
from IPython.display import HTML, clear_output, display

import pandas as pd
from Bio.PDB import PDBParser

# Global variables
XYZ = ['x', 'y', 'z']

__all__ = [
        'donuts', 'read_pdb', 'get_distance_matrix', 'plot_3d', 'animate_ax'
        ]

def donuts(df, numcols=2, thres=2, maxvals=15, title=None, ytitle=None, lw=6):
    """
    Creates a donut chart for every categorical variable in a pandas dataframe.
    
    Parameters
    ----------
    df : pandas dataframe containing the variables to be represented as donut charts
    numcols : number of columns for the subplots to be generated
    thres : threshold for the minimum percentage a value must represent in order to be included in the chart individually
    maxvals : maximum number of values to be included individually without being compared against the percentage threshold
    title : main title of the figure

    Returns
    ----------
    fig : matplotlib figure containing the donut graphs
    
    """
    # Recorrer cada variable categórica
    vars_categ = list(df.select_dtypes(include='object').columns)
    numrows = math.ceil(len(vars_categ) / float(numcols))
    numcols = numcols*2
    Height = numrows*5
    Width = numcols*5
    fig = plt.figure(figsize=(Width, Height))
    for i in range(len(vars_categ)):
        col = vars_categ[i]
        # Para cada variable categórica, agrupar los valores contando las ocurrencias de cada valor.
        counts = df[col].value_counts()
        counts = counts.apply(float)
        # Para facilitar la generación de la gráfica y visualización de variables con muchos valores diferentes,
        # agrupar los valores poco frecuentes en la categoría "Otros <[threshold definido en Thres]%"
        # en las variables con más de "maxvals" valores diferentes
        if len(counts) >= maxvals:
            tot = float(counts.sum())
            counts = pd.concat([counts[counts*100/tot>=thres],pd.Series({"Other (<"+str(thres)+"%)": counts[counts*100/tot<thres].sum()})])
        # Generar las gráficas de pay.
        # Poner el título del par de gráficos actual.
        idxSubplot = (i * 2) + 1
        plt.subplot(numrows,numcols,idxSubplot)
        plt.title(col, fontsize="x-large").set_position([.8, 1])
        plt.pie(counts, startangle=40, wedgeprops = {'linewidth': lw, 'edgecolor': 'white'},
                autopct="%.2f%%", pctdistance=1.2)
        plt.legend(counts.index, loc="center left", bbox_to_anchor=(1.1, 0, 0.5, 1))
        # Dibujar un círculo para que sea una gráfica de dona (en vez de pay).
        plt.gca().add_artist(plt.Circle( (0,0), 0.65, color='white'))
    # Título principal
    if ytitle == None:
        posy = 0.88 + 1/(5*numrows)
    else:
        posy = ytitle
    fig.suptitle(title, fontsize=25, y=posy)
    plt.close()
    return fig

def read_pdb(pdb_file):
    '''Take a PDB file and convert it to a pandas dataframe'''
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("structure", pdb_file)
    data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    data.append({"modelNo": model.id,
                                 "chain": chain.id,
                                 "residue": residue.get_id()[1],
                                 "amino_acid": residue.get_resname(),
                                 "atom_name": atom.get_name(),
                                 "x": atom.get_coord()[0],
                                 "y": atom.get_coord()[1],
                                 "z": atom.get_coord()[2],
                                 })
    df = pd.DataFrame(data)
    return df

def get_distance_matrix(df, only_CA=True):
    """Calculate the distance matrix for a dataframe with x, y, z coordinates.
    Parameters
    ----------
    df : pandas dataframe
        A pandas DataFrame containing the points' coordinates. Each row represents a point,
        and the columns represent the x, y, and z coordinates respectively.
    only_CA : bool, optional
        If True, only the alpha carbon (CA) atoms will be considered. Defaults to True.
    Returns
    -------
    dists : pandas DataFrame
        A pandas DataFrame containing the distance matrix. The index and columns are the same as the input DataFrame."""
    # Make sure df only has alpha carbons
    if only_CA and "atom_name" in df:
        df = df[df.atom_name == "CA"]
    # Or keep only one atom per residue
    if only_CA and 'residue' in df:
        df = df.drop_duplicates('residue')
        df = df.set_index('residue')
    elif not only_CA and 'residue' in df:
        df = df.set_index(['residue', 'atom_name'])
    # else:
    #     df = df.reset_index(drop=True)
    # calculate the distance matrix with the function from scipy.spatial
    dists = pd.DataFrame(distance_matrix(df[['x','y','z']].values, df[['x','y','z']].values),
                         index=df.index, columns=df.index)
    # Reduce to the minimum distance between the residues
    try:
        dists = dists.groupby('residue').min().T.groupby('residue').min().T
    except Exception as e:
        print("get_distance_matrix error: {}".format(e))
    return dists
# plot 3D
def plot_3d(points_df=None, ref_s=None, plane_coef=None,
            ax=None, only_ca=True, title=None, naked=True,
            animate=False,
            savegif=False, gifname=None):
    """
    Plots a 3D scatter plot of points (C-alpha) and a reference set (if provided).
    Parameters
    ----------
    points_df : DataFrame, optional
        A pandas DataFrame containing the points' coordinates. Each row represents a point,
        and the columns represent the x, y, and z coordinates respectively. Defaults to None.

    ref_s : array-like, optional
        A reference point to indicate the origin of the coordinate system. The reference point
        will be plotted as a red dot. Defaults to None.

    plane_coef : tuple, optional
        A tuple of three coefficients representing a plane's equation (A, B, C, D) in the form Ax + By + Cz + D = 0.
        If provided, the plane will be plotted as a transparent surface. Defaults to None.

    ax : matplotlib.axes.Axes, optional
        The matplotlib axes on which to plot. If not provided, a new figure and axes will be created.
        Defaults to None.

    only_ca : bool, optional
        If True, only the alpha carbon (CA) atoms will be plotted. Defaults to True.

    title : str, optional
        The title of the plot. Defaults to None.

    animate : bool, optional
        If True, the plot will be animated. Defaults to True.

    savegif : bool, optional
        If True, the animated plot will be saved as a GIF file. Defaults to False.

    gifname : str, optional
        The name of the GIF file to be saved. Required if savegif is True. Defaults to None.
    
    naked : bool, optional
        If True, remove the axes and grid from the plot. Defaults to True.

    Returns
    -------
    ax

    Example
    -------
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    points = np.array([[0,0,0], [1,1,1], [2,2,2]])
    points_df = pd.DataFrame(points, columns=['x', 'y', 'z'])

    plot_3d(points_df)
    plt.show()
    """
    if ax is None:
        fig = plt.figure(facecolor='none')
        ax = fig.add_subplot(111, projection='3d', facecolor='none')
        
    # Plot the points (C-alpha)
    if points_df is not None:
        if only_ca:
            if 'atom_name' in points_df:
                points_df = points_df[points_df.atom_name == "CA"]
            elif "residue" in points_df:
                print('No atom_name found; dropping residue "duplicates" instead.')
                points_df = points_df.drop_duplicates("residue")
            else:
                print('No atom_name or residue index found.')
        ax.scatter(points_df['x'], points_df['y'], points_df['z'])
    if ref_s is not None:
        ax.scatter(ref_s['x'], ref_s['y'], ref_s['z'], c="r", depthshade=True)
    
    # Get z-axis limits to re-establish them after plotting the plane
    zlims = ax.get_zlim3d()
    
    # plot the plane
    if plane_coef is not None:
        plot_plane(plane_coef, ref_s, ax)
    
    # Set back the z limits to avoid plane overstretching the plot
    ax.set_zlim3d(zlims)
    
    # Customize the plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    # Make the aspect ratio 1:1:1 in data space
    ax.set_box_aspect((np.ptp(ax.get_xlim3d()), np.ptp(ax.get_ylim3d()), np.ptp(ax.get_zlim3d())))
    # Remove axes and grid
    if naked:
        ax.axis('off')
        ax.grid(False)
        plt.tight_layout()
    # Generate an animated plot
    if animate or savegif:
        animate_ax(ax, animate, savegif, gifname)
    return ax

def animate_ax(ax, show=True, savegif=False, gifname=None, dpi=100):
    """Called inside plot_3d"""
    # Define the update function for the animation
    def update(frame):
        ax.view_init(elev=15, azim=frame);  # Rotate the plot
        return ax
    
    fig = ax.get_figure()
    #fig.patch.set_alpha(0.0)
    fig.set_size_inches(4,4)

    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=np.arange(0, 360, 2), interval=100);

        # Save the animation as a GIF
    if savegif:
        if not os.path.isdir('figures'):
            os.mkdir('figures')
        if gifname is None:
            gifname = "plot_"+datetime.now().isoformat(sep="T", timespec="seconds").replace(':','.')+'.gif'
        ani.save('figures/'+gifname, writer='pillow', dpi=dpi, fps=25, 
                 #savefig_kwargs={"transparent": True, "facecolor": "none"}
                )
    
    # Show the animation
    #clear_output(wait=False)
    if show:
        plt.close()
        display(HTML(ani.to_jshtml(default_mode='loop', fps=25)))
    return
