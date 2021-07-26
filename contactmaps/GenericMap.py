# ======================================================================
#  AUTHOR:  Ben Rafferty, Purdue University
#  Copyright (c) 2010  Purdue Research Foundation
#
#  See the file "LICENSE.txt" for information on usage and
#  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# ======================================================================

import Bio.PDB
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy
import io
from .utils import remove_hetero

class GenericMap():
    """
    Generic class for 2D protein data maps.  To define a specific map, create
    a child class and override calc_matrix(model).  This class cannot be
    instantiated directly, you must use a child class.  Any additional keyword
    arguments given to the constructor will be passed to calc_matrix.
    """
    def __init__(self, model, title="", xlabel="", ylabel="", colorbar=True, \
                 colorbarlabel="", contour=False, interpolation="nearest", \
                 cmap="jet", **kwargs):
        if not isinstance(model, Bio.PDB.Model.Model):
            raise TypeError("Input argument is not of class Bio.PDB.Model")
        if len(model) < 1:
            raise ValueError("No chains found in model.")
        remove_hetero(model)
        self.matrix = self.calc_matrix(model, **kwargs)
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cmap = cmap
        self.colorbar = colorbar
        self.colorbarlabel = colorbarlabel
        self.contour = contour
        self.interpolation = interpolation

    def calc_matrix(model, **kwargs):
        """
        Override this method in any child classes.  Must perform some
        calculation on the provided model, and return a matrix representing
        the data for the map.
        """
        raise NotImplementedError( \
              "GenericMap cannot be instantiated directly.")

    def get_data(self):
        """
        Returns the matrix containing the map data.
        """
        return self.matrix

    def show(self):
        """
        Use pyplot to display the map in an interactive window.  Will block
        until the window is closed.
        """
        import matplotlib.pyplot as plt 
        self.draw(fig=plt.figure())
        plt.show()

    def print_figure(self, filename="map.png"):
        """
        Generates a plot of the map's data and saves an image at the given
        filename.
        """
        canvas = FigureCanvas(self.draw())
        canvas.print_figure(filename)

    def draw(self, fig=None):
        """
        Creates a matplotlib figure representing the map's matrix.  An existing
        figure can be provided to be drawn on, otherwise a new figure will be
        created.
        """
        if not fig: fig = Figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(self.matrix, cmap=matplotlib.cm.get_cmap(self.cmap), \
                         interpolation=self.interpolation, origin="lower")
        ax.set_title(self.title)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        if self.colorbar: 
            cbar = fig.colorbar(cax, cmap=matplotlib.cm.get_cmap(self.cmap))
            cbar.ax.set_ylabel(self.colorbarlabel)
        if self.contour:
            ax.contour(self.matrix, cmap=matplotlib.cm.get_cmap(self.cmap))
        return fig

    def tabsep_data(self):
        """
        Returns a string containing the matrix data as a tab separated
        data table.
        """
        f = io.StringIO()
        numpy.savetxt(f, self.matrix, delimiter='\t', newline='\n')
        ret = f.getvalue()
        f.close()
        return ret

    def __repr__(self):
        return repr(self.matrix)