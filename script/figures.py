""" #TODO: fix the documentation of the module
This module can be used for peaks of an NMR spectrum.

Classes:
    Peak: for peaks of the spectrum.
    PeakList: for a list of peaks of the spectrum.
    PeakOverlap: for the detection of the peak overlap.
"""

from __future__ import annotations

import logging
import os
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

logger = logging.getLogger(__name__)

# DO NOT CHANGE 'plot_settings_default'
# create a a new dictionary instead, 
# and specify the new dictionary when calling functions to create figures
plot_settings_default = {
        "figure_width": 8,
        "figure_height": 6,
        "color": 'b',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.05),
        "yaxis_label_coords": (-0.05, 1),
        "invert_xaxis": False,
        "invert_yaxis": False,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_left": 0.16,
        "subplots_adjust_right": 0.93,
        "subplots_adjust_bottom": 0.16,
        "subplots_adjust_top": 0.95,
        "scatter_color": '#00ff0080',
    }

plot_settings_cross_section_for_paper = {
        "figure_width": 6,#4.3,#6,
        "figure_height": 4.5,#3.225,#3.9,
        "color": 'b',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.05),
        "yaxis_label_coords": (-0.1, 1),#(-0.05, 1),
        "invert_xaxis": False,
        "invert_yaxis": False,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_left": 0.2,#0.16,
        "subplots_adjust_right": 0.93,
        "subplots_adjust_bottom": 0.16,
        "subplots_adjust_top": 0.90,#0.95,
    }

plot_settings_profile_for_paper = {
        "figure_width": 6,#4.35,#6,
        "figure_height": 3.9,#4.5,#3.26,#3.9,
        "color": 'b',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.12),
        "yaxis_label_coords": (-0.19, 0.95),#(-0.13, 0.95),
        "invert_xaxis": True,
        "invert_yaxis": False,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_left": 0.2,#0.16,
        "subplots_adjust_right": 0.93,
        "subplots_adjust_bottom": 0.16,
        "subplots_adjust_top": 0.90,#0.95,
    }

plot_settings_combined_plots_for_paper = {
        "figure_width": 17,
        "figure_height": 21,
        "color": 'b',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.12),
        "yaxis_label_coords": (-0.13, 0.95),
        "invert_xaxis": True,
        "invert_yaxis": False,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_wspace": 0.28,
        "subplots_adjust_hspace": 0.45,
        "subplots_adjust_left": 0.08,
        "subplots_adjust_right": 0.95,
        "subplots_adjust_bottom": 0.05,
        "subplots_adjust_top": 0.95,
        "scatter_color": '#00ff00',
        "scatter_transparency": 0.5,
    }

plot_settings_cross_section = {
        "figure_width": 8,
        "figure_height": 6,
        "color": 'b',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.05),
        "yaxis_label_coords": (-0.05, 1),
        "invert_xaxis": False,
        "invert_yaxis": False,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_left": 0.16,
        "subplots_adjust_right": 0.93,
        "subplots_adjust_bottom": 0.16,
        "subplots_adjust_top": 0.95,
    }

plot_settings_profile = {
        "figure_width": 8,
        "figure_height": 6,
        "color": 'b',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.05),
        "yaxis_label_coords": (-0.05, 1),
        "invert_xaxis": False,
        "invert_yaxis": False,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_left": 0.16,
        "subplots_adjust_right": 0.93,
        "subplots_adjust_bottom": 0.16,
        "subplots_adjust_top": 0.95,
    }


plot_settings_reference_comparison = {
        "figure_width": 6,#4.35,#6,
        "figure_height": 4.5,#3.26,#3.9,
        "color": 'g',
        "linewidth": 0.7,
        "border_linewidth": 0.5,
        "xaxis_label_coords": (1, -0.07),
        "yaxis_label_coords": (-0.13, 0.95),#(-0.13, 0.95),
        "invert_xaxis": True,
        "invert_yaxis": True,
        "tick_labelsize": 6,
        "tick_length": 2,
        "tick_pad": 0,
        "tick_width": 0.5,
        "subplots_adjust_left": 0.15,#0.16,
        "subplots_adjust_right": 0.93,
        "subplots_adjust_bottom": 0.11,
        "subplots_adjust_top": 0.95,#0.95,
        "scatter_color": 'b',
        "scatter_transparency": 1,
        "scatter_size": 1,
    }



class ScalarFormatterForceFormat(ticker.ScalarFormatter):
    def _set_format(self):
        self.format = "%1.1f"
        #self.set_powerlimits((1,1))


#TODO: think about moving this function somewhere else?
def find_closest_value(x_values, y_values, correct_shift):
    #correct_shift
    #x_values
    #if correct_shift == None: return None
    compare_value = abs(x_values[0] - correct_shift)
    for i in range(0, len(x_values)):
        distance_from_correct = abs(x_values[i] - correct_shift)
        if distance_from_correct < compare_value:
            closest_index = i
            compare_value = distance_from_correct
    plot_value = y_values[closest_index]
    return plot_value

def create_cross_section_fig_file(
    x_values,
    y_values,
    fig_file_path,
    width: float = 8,
    height: float = 6,
    ) -> None:
    """docstring
    width in cm
    height in cm
    """
    cm = 1/2.54
    fig = plt.figure(figsize=(width * cm, height * cm))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x_values, y_values, color = 'b', linewidth = 0.7)
    ax.set_yticks(ticks = [0])
    for border in ['top', 'bottom', 'left', 'right']: ax.spines[border].set_linewidth(0.5)
    
    ax.set_xticks(ticks=[len(x_values)/2], labels=['0'])
    
    #ax.set_xticks(ticks=[len(x_values)/2, peak_pos], labels=['0', "$\Omega_1$"])
    ax.set_xlabel(r"$\omega_1$", loc='right', fontsize=6)
    ax.set_ylabel("intensity", loc='top', fontsize=6)
    ax.tick_params(which='both', labelsize=6, length = 2, pad = 0, width = 0.5)
    #ax.xaxis.set_inverted(True) #if chem shifts
    ax.xaxis.set_label_coords(1, -0.05) #the same as for profiles
    ax.yaxis.set_label_coords(-0.05, 1)
    plt.subplots_adjust(left=0.16, right=0.93, bottom=0.16, top=0.95) #the same as for profiles
    fig.savefig(fig_file_path, dpi=300)
    plt.close()
    
    
def create_one_large_fig_file( #TODO: take list of 2D/cross_section/profile objects?
    list_of_x_values,
    list_of_y_values,
    fig_file_dir,
    correct_shifts,
    names,
    ) -> None:
    fig_counter = 0
    cm = 1/2.54
    for i, x_values in enumerate(list_of_x_values):
        name = names[i]
        #correct_shifts_name = correct_shifts[i][0]
        #if name != correct_shifts_name:
        #    exit(f"Names do not match: {name}, {correct_shifts_name}")
        #correct_shift = correct_shifts[i][1]
        correct_shift = correct_shifts[i]
        #print(name)
        name = name.replace("i_1", "")
        name = name.replace("i", "")
        #print(name)
        #x_values = profile[1]
        #y_values = profile[2]
        y_values = list_of_y_values[i]
        
        if (i) % 35 == 0:# or i == len(profiles_full_list) - 1:
            fig = plt.figure(figsize=(17*cm,24*cm))
            #img_path = os.path.join(dir_path, f"Figure{fig_counter}_alt.png")
            fig_file_name = f"PROFILES_{fig_counter}.png"
            #fig_file_path.replace(".png", f"{fig_counter}.png")
            fig_file_path = os.path.join(fig_file_dir, fig_file_name)
            fig_counter += 1
        
        #print(i, (i+1)%35)
        index = i%35 + 1
        #index = i 
        ax = fig.add_subplot(7, 5, index)
        ax.set_title(name, fontsize=7, pad=2)
        if correct_shift != 0:
            plot_value = find_closest_value(x_values, y_values, correct_shift)
            ax.scatter(correct_shift, plot_value, color='#00ff0080')
        ax.plot(x_values,y_values,color = 'b', linewidth = 0.7)
        ax.set_yticks(ticks=[0, 50, 100])
        #ax.set_xticks(ticks=[170, 176.5, 183])
        ax.set_xticks(ticks=[169.5, 174, 178.5, 183])
        ax.set_xlim(169.5, 183)
        ax.set_ylim(None, 100)
        ax.tick_params(labelsize=6, length = 3, pad = 2)
        ax.xaxis.set_inverted(True)
        
        #if (i) % 35 == 34 or i == len(profiles_full_list) - 1:
        if (i) % 35 == 34 or i == len(list_of_x_values) - 1:
            plt.subplots_adjust(wspace=0.28,hspace=0.35, left=0.08, right=0.95, bottom=0.05, top=0.95)
            #fig.show()
            fig.savefig(fig_file_path, dpi=300)
            plt.close()



def create_cross_section_fig_file_v2(
    x_values,
    y_values,
    fig_file_path,
    width: float = 8,
    height: float = 6,
    ) -> None:
    """docstring
    width in cm
    height in cm
    """
    cm = 1/2.54
    fig_counter = 0

    if isinstance(x_values[0], list):
        for x_values_list in x_values:
            pass
    else:
        fig = plt.figure(figsize=(width * cm, height * cm))
        ax = fig.add_subplot(1, 1, 1)

    
    #do it in one loop, establish the number of loops based on if its a list and several figures or just one figure
    ax.plot(x_values, y_values, color = 'b', linewidth = 0.7)
    ax.set_yticks(ticks = [0])
    for border in ['top', 'bottom', 'left', 'right']: ax.spines[border].set_linewidth(0.5)
    
    ax.set_xticks(ticks=[len(x_values)/2], labels=['0'])
    
    #ax.set_xticks(ticks=[len(x_values)/2, peak_pos], labels=['0', "$\Omega_1$"])
    ax.set_xlabel(r"$\omega_1$", loc='right', fontsize=6)
    ax.set_ylabel("intensity", loc='top', fontsize=6)
    ax.tick_params(which='both', labelsize=6, length = 2, pad = 0, width = 0.5)
    #ax.xaxis.set_inverted(True) #if chem shifts
    ax.xaxis.set_label_coords(1, -0.05) #the same as for profiles
    ax.yaxis.set_label_coords(-0.05, 1)
    plt.subplots_adjust(left=0.16, right=0.93, bottom=0.16, top=0.95) #the same as for profiles
    fig.savefig(fig_file_path, dpi=300)
    plt.close()


def create_2Dplot_single_file(
    x_values,
    y_values,
    fig_file_path,
    title: str = None,
    plot_type: str = None,
    #width: float = 8,
    #height: float = 6,
    #color: any = 'b',
    #plot_linewidth = 0.7,
    #border_linewidth = 0.5,
    #predefined_plot_settings = None,
    transparent = False,
    #xticks: dict[float: float|str] = None, #tick_value: label
    xticks: list[float] = None,
    xticks_labels: list[str] = None,
    yticks: list[float] = None,
    yticks_labels: list[str] = None,
    #yticks: dict[float: float|str] = None, #tick_value: label
    xlim: tuple[float, float] = None,
    ylim: tuple[float, float] = None,
    xlabel: str = "x",
    ylabel: str = "y",
    plot_settings: dict[str: any] = plot_settings_default,
    ) -> None:
    """docstring
    width in cm
    height in cm
    """

    cm = 1/2.54
    fig = plt.figure(figsize = (plot_settings.get("width") * cm, 
        plot_settings.get("height") * cm))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x_values, y_values, color = plot_settings.get("color"), 
        linewidth = plot_settings.get("linewidth"))
    
    if title != None:
        ax.set_title(title, fontsize=7, pad=2)
    
    for border in ['top', 'bottom', 'left', 'right']: 
        ax.spines[border].set_linewidth(plot_settings.get("border_linewidth"))

    if xticks != None:
        #ax.set_xticks(ticks = list(xticks.keys()),
            #labels = list(xticks.values()))
        ax.set_xticks(ticks = xticks, labels = xticks_labels)
    if yticks != None:
        #ax.set_yticks(ticks = list(yticks.keys()),
            #labels = list(yticks.values()))
        ax.set_yticks(ticks = yticks, labels = yticks_labels)
    if xlim != None: ax.set_xlim(*xlim)
    if ylim != None: ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel = xlabel, loc = 'right', fontsize = 6)
    ax.set_ylabel(ylabel = ylabel, loc = 'top', fontsize = 6)
    ax.xaxis.set_label_coords(*plot_settings.get("xaxis_label_coords")) #the same as for profiles
    ax.yaxis.set_label_coords(*plot_settings.get("yaxis_label_coords"))
    ax.xaxis.set_inverted(plot_settings.get("invert_xaxis", False))
    ax.yaxis.set_inverted(plot_settings.get("invert_yaxis", False))
    ax.tick_params(which = 'both', 
        labelsize = plot_settings.get("tick_labelsize", 6), 
        length = plot_settings.get("tick_length", 2), 
        pad = plot_settings.get("tick_pad", 0), 
        width = plot_settings.get("tick_width", 0.5))
    if True: pass
    #if predefined_plot_settings != None:
    #    ax = predefined_plot_settings(ax)    
    elif plot_type == "cross_section":
        pass
        #ax.set_xticks(ticks=[len(x_values)/2], labels=['0'])
        #ax.set_xticks(ticks=[len(x_values)/2, peak_pos], labels=['0', "$\Omega_1$"])
        #ax.set_yticks(ticks = [0])
        #ax.set_xlim(x_values[0], x_values[-1])
        #ax.set_ylim(None, 100)
        #ax.set_xlabel("$\omega_1$", loc='right', fontsize=6)
        #ax.set_ylabel("intensity", loc='top', fontsize=6)
        #ax.xaxis.set_label_coords(1, -0.05) #the same as for profiles
        #ax.yaxis.set_label_coords(-0.05, 1)
        
    elif plot_type == "profile":
        pass
        #ax.set_xticks(ticks=[168.5, 173.5, 178.5, 183.5]) #change xticks and xlim based on function input
        #ax.set_yticks(ticks=[0, 50, 100])
        #ax.set_xlim(168.5, 183.5)
        #ax.set_ylim(None, 100)
        #ax.set_xlabel("chemical shift [ppm]", loc='right', fontsize=6)
        #ax.set_ylabel("relative intensity [%]", loc='top', fontsize=6)
        #ax.xaxis.set_inverted(True) #if chem shifts
        #ax.xaxis.set_label_coords(1, -0.12)
        #ax.yaxis.set_label_coords(-0.13, 0.95)
    else:
        pass #general settings
    #ax.tick_params(which='both', labelsize=6, length = 2, pad = 0, width = 0.5)
    #ax.xaxis.set_label_coords(1, -0.05) #the same as for profiles
    #ax.yaxis.set_label_coords(-0.05, 1)
    plt.subplots_adjust(
        left = plot_settings.get("subplots_adjust_left", 0.16), 
        right = plot_settings.get("subplots_adjust_right", 0.93), 
        bottom = plot_settings.get("subplots_adjust_bottom", 0.16), 
        top = plot_settings.get("subplots_adjust_top", 0.95)
        )
    fig.savefig(fig_file_path, dpi = 300, transparent = transparent)
    plt.close()


def create_2Dplots( #TODO: take list of 2D/cross_section/profile objects?
    list_of_x_values: list[list[float]],
    list_of_y_values: list[list[float]],
    fig_file_dir: os.PathLike,
    fig_file_basename: str,
    number_of_rows: int = 1,
    number_of_columns: int = 1,
    figure_width: float = None,
    figure_height: float = None,
    #correct_points: list[float] = None,
    names: list[str] = None,
    titles: list[str] = None,
    xlabel: str = None,
    ylabel: str = None,
    xlim: tuple[float, float] = None,
    ylim: tuple[float, float] = None,
    xticks: list[float] = None,
    xticks_labels: list[str] = None,
    yticks: list[float] = None,
    yticks_labels: list[str] = None,
    transparent: bool = False,
    plot_settings: dict[str: any] = plot_settings_default,
    list_of_scatter_x_values: list[list[float]] = None,
    list_of_scatter_y_values: list[list[float]] = None,
    xlabel_global: str = None,
    ylabel_global: str = None,
    ) -> None:
    
    cm = 1/2.54
    fig_counter = 0
    number_of_plots_in_one_figure = number_of_rows * number_of_columns
    if figure_width == None:
        figure_width = plot_settings.get("figure_width", 8)
    if figure_height == None:
        figure_height = plot_settings.get("figure_height", 6)

    if len(list_of_x_values) != len(list_of_y_values):
        raise ValueError(f"The lists of x and y values have different lengths:"
            f" x: {len(list_of_x_values)}, y: {len(list_of_y_values)}")

    for plot_counter, x_values in enumerate(list_of_x_values):


        y_values = list_of_y_values[plot_counter]
        
        #correct_shifts_name = correct_shifts[i][0]
        #if name != correct_shifts_name:
        #    exit(f"Names do not match: {name}, {correct_shifts_name}")
        #correct_shift = correct_shifts[i][1]
        #correct_shift = correct_shifts[i]
        
        # start new figure when the previous one is already full
        if plot_counter % number_of_plots_in_one_figure == 0:



            fig = plt.figure(figsize=(figure_width * cm, figure_height *cm))
            fig.tight_layout()
            #fig, axs = plt.subplots(number_of_rows, number_of_columns, figsize=(figure_width * cm, figure_height *cm))
            if len(list_of_x_values) == 1:
                fig_file_name = f"{fig_file_basename}.png"
            else:
                fig_file_name = f"{fig_file_basename}_fig{fig_counter}.png"
            fig_file_path = os.path.join(fig_file_dir, fig_file_name)
            fig_counter += 1
            #axs=axs.flatten()
        
        # plotting the subplot
        #subplot_index = plot_counter % 35 + 1
        subplot_index = plot_counter % number_of_plots_in_one_figure + 1
        
        ax = fig.add_subplot(number_of_rows, number_of_columns, subplot_index)
        #ax=axs[subplot_index-1]
        ax.plot(x_values, y_values, color = plot_settings.get("color"), 
            linewidth = plot_settings.get("linewidth"))
        #plot = plt.subplot(number_of_rows, number_of_columns, subplot_index)
        #plot.set_size_inches(2, 2)
        
        if list_of_scatter_x_values != None:
            scatter_x_values = list_of_scatter_x_values[plot_counter]
            if list_of_scatter_y_values != None:
                scatter_y_values = list_of_scatter_y_values[plot_counter]
            else:
                scatter_y_values = []
                for scatter_x in scatter_x_values:
                    if scatter_x == None:
                        scatter_y = None
                    else:
                        scatter_y = find_closest_value(
                            x_values, y_values, scatter_x)
                    scatter_y_values.append(scatter_y)
            ax.scatter(scatter_x_values, scatter_y_values, 
                color = plot_settings.get("scatter_color"), 
                alpha = plot_settings.get("scatter_transparency", 1),
                s = plot_settings.get("scatter_size", None))
        
        """if correct_shift != 0:
            plot_value = find_closest_value(x_values, y_values, correct_shift)
            ax.scatter(correct_shift, plot_value, color='#00ff0080')
        """
        
        # settings for title and axis labels
        if titles != None:
            ax.set_title(titles[plot_counter], fontsize = 7, y = 1.05)#, pad = 2)
        if xlabel != None:
            ax.set_xlabel(xlabel = xlabel, loc = 'right', fontsize = 6)
            ax.xaxis.set_label_coords(*plot_settings.get("xaxis_label_coords"))
        if ylabel != None:
            ax.set_ylabel(ylabel = ylabel, loc = 'top', fontsize = 6)
            ax.yaxis.set_label_coords(*plot_settings.get("yaxis_label_coords"))
        
        # setting axis global labels
        if xlabel_global != None:
            fig.supxlabel(xlabel_global, fontsize = 8)
        if ylabel_global != None:
            fig.supylabel(ylabel_global, fontsize = 8)
        
        # setting axis limits and directions
        if xlim != None: ax.set_xlim(*xlim)
        if ylim != None: ax.set_ylim(*ylim)
        ax.xaxis.set_inverted(plot_settings.get("invert_xaxis", False))
        ax.yaxis.set_inverted(plot_settings.get("invert_yaxis", False))
        
        # setting ticks
        single_decimal_scientific_format = ScalarFormatterForceFormat()
        single_decimal_scientific_format.set_powerlimits((0,0))
        if xticks != None:
            ax.set_xticks(ticks = xticks, labels = xticks_labels)
        if yticks != None:
            ax.set_yticks(ticks = yticks, labels = yticks_labels)
        else:
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, steps = [1,2,4,10], min_n_ticks = 4))
            ax.yaxis.set_major_formatter(single_decimal_scientific_format)
        '''else:
            maximum = max(y_values)
            minimum = min(y_values)
            y_range = maximum - minimum
            y_range_str = str(y_range)
            decimal_pos = y_range_str.find(".")
            if decimal_pos == -1:
                
            precision = 
            increment = y_range / 5
            minimum + '''

        ax.tick_params(which = 'both', 
            labelsize = plot_settings.get("tick_labelsize", 6), 
            length = plot_settings.get("tick_length", 2), #3
            pad = plot_settings.get("tick_pad", 0), #2
            width = plot_settings.get("tick_width", 0.5),
            labelcolor = "black")
        
        #ax.yaxis.set_major_locator(ticker.MultipleLocator(20000))
        
        '''ax.ticklabel_format(axis='y', 
            style='scientific', 
            scilimits=(0,0),
            #useMathText=True
            )'''
        #ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.1f}"))
        #ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.1f}"))
        ax.yaxis.get_offset_text().set_fontsize(6)
        ax.yaxis.get_offset_text().set_position((0, 0.95))
        '''ax.yaxis.set_major_formatter(
            ticker.LogFormatterExponent(base=10.0, labelOnlyBase=True)
        )'''
        
        # setting border linewidth:
        for border in ['top', 'bottom', 'left', 'right']: 
            ax.spines[border].set_linewidth(plot_settings.get("border_linewidth"))    
        
        # if there are enough plots in one figure or
        # we reach the last plot
        if (plot_counter % number_of_plots_in_one_figure == 
            number_of_plots_in_one_figure - 1 
            or 
            plot_counter == len(list_of_x_values) - 1):

            plt.subplots_adjust(
                wspace = plot_settings.get("subplots_adjust_wspace", 0.28), 
                hspace = plot_settings.get("subplots_adjust_hspace", 0.35), 
                left = plot_settings.get("subplots_adjust_left", 0.08), 
                right = plot_settings.get("subplots_adjust_right", 0.95), 
                bottom = plot_settings.get("subplots_adjust_bottom", 0.05), 
                top = plot_settings.get("subplots_adjust_top", 0.8),#0.95)
                )

            fig.savefig(fig_file_path, dpi = 300)#, transparent = transparent)
            plt.close(fig)








def create_2Dplots_v2( #TODO: take list of 2D/cross_section/profile objects?
    list_of_x_values: list[list[float]],
    list_of_y_values: list[list[float]],
    fig_file_dir: os.PathLike,
    fig_file_basename: str,
    number_of_rows: int = 1,
    number_of_columns: int = 1,
    figure_width: float = None,
    figure_height: float = None,
    #correct_points: list[float] = None,
    names: list[str] = None,
    titles: list[str] = None,
    xlabel: str = None,
    ylabel: str = None,
    xlim: tuple[float, float] = None,
    ylim: tuple[float, float] = None,
    xticks: list[float] = None,
    xticks_labels: list[str] = None,
    yticks: list[float] = None,
    yticks_labels: list[str] = None,
    transparent: bool = False,
    plot_settings: dict[str: any] = plot_settings_default,
    list_of_scatter_x_values: list[list[float]] = None,
    list_of_scatter_y_values: list[list[float]] = None,
    xlabel_global: str = None,
    ylabel_global: str = None,
    plot_colors: list[str] = None,
    first_figure_plots_number: int = None,
    ) -> None:

    def adjust_height(
        number_of_rows: int,
        number_of_columns: int,
        remaining_figures: int,
        number_of_plots_in_one_figure: int,
        figure_height: float,
        top_margin: float,
        bottom_margin: float,
        plot_height: float,
        vertical_plot_spacing: float,
        flag: bool,
        ) -> None:
        orig_number_of_rows = number_of_rows
        new_number_of_rows = int(remaining_figures / number_of_plots_in_one_figure * number_of_rows)
        if remaining_figures > new_number_of_rows * number_of_columns:
            number_of_rows = new_number_of_rows + 1
        else:
            number_of_rows = new_number_of_rows
        orig_height = figure_height
        figure_height = (top_margin+bottom_margin)*figure_height + number_of_rows * plot_height + vertical_plot_spacing * (number_of_rows - 1)#number_of_rows / orig_number_of_rows * (figure_height*1.2)
        top_margin = top_margin * orig_height / figure_height
        bottom_margin = bottom_margin * orig_height / figure_height
        subplots_adjust_top = 1 - top_margin
        subplots_adjust_bottom = bottom_margin
        flag = True
        print(top_margin)
        return number_of_rows, figure_height, top_margin, bottom_margin, subplots_adjust_top, subplots_adjust_bottom, flag







    
    cm = 1/2.54
    fig_counter = 0
    number_of_plots_in_one_figure = number_of_rows * number_of_columns
    if figure_width == None:
        figure_width = plot_settings.get("figure_width", 8)
    if figure_height == None:
        figure_height = plot_settings.get("figure_height", 6)

    if len(list_of_x_values) != len(list_of_y_values):
        raise ValueError(f"The lists of x and y values have different lengths:"
            f" x: {len(list_of_x_values)}, y: {len(list_of_y_values)}")


    plot_width = 2.4 #* cm
    plot_height = 2.7 #* cm
    #horizontal_plot_spacing
    #vertical_plot_spacing
    left_margin = 0.07
    right_margin = 0.05
    top_margin = 0.05
    bottom_margin = 0.05
    horizontal_plot_spacing = (figure_width - plot_width * number_of_columns - left_margin*figure_width - right_margin*figure_width) / (number_of_columns - 1)
    vertical_plot_spacing = (figure_height - plot_height * number_of_rows - top_margin*figure_height - bottom_margin*figure_height) / (number_of_rows - 1)

    subplots_adjust_wspace = horizontal_plot_spacing / plot_width
    subplots_adjust_hspace = vertical_plot_spacing / plot_height
    subplots_adjust_left = left_margin
    subplots_adjust_right = 1 - right_margin
    subplots_adjust_top = 1 - top_margin
    subplots_adjust_bottom = bottom_margin

    orig_number_of_rows = number_of_rows
    orig_figure_height = figure_height
    orig_top_margin = top_margin
    orig_bottom_margin = bottom_margin
    orig_subplots_adjust_top = subplots_adjust_top
    orig_subplots_adjust_bottom = subplots_adjust_bottom
    orig_number_of_plots_in_one_figure = number_of_plots_in_one_figure



    already_adjusted = False #NEW
    first_fig_done = False #NEW
    first_fig_done_fully = False
    for plot_counter, x_values in enumerate(list_of_x_values):


        y_values = list_of_y_values[plot_counter]
        
        #correct_shifts_name = correct_shifts[i][0]
        #if name != correct_shifts_name:
        #    exit(f"Names do not match: {name}, {correct_shifts_name}")
        #correct_shift = correct_shifts[i][1]
        #correct_shift = correct_shifts[i]
        
        # start new figure when the previous one is already full
        #if plot_counter % number_of_plots_in_one_figure == 0:
        if first_figure_plots_number is None:
            counter = plot_counter
        elif first_fig_done_fully:
            counter = plot_counter - first_figure_plots_number
        else:
            counter = plot_counter
        if counter % number_of_plots_in_one_figure == 0:

            #NEW PART:
            remaining_figures = len(list_of_x_values) - plot_counter

            if first_figure_plots_number is not None and not first_fig_done_fully:
                number_of_rows, figure_height, top_margin, bottom_margin, subplots_adjust_top, subplots_adjust_bottom, first_fig_done = adjust_height(number_of_rows, number_of_columns, first_figure_plots_number, number_of_plots_in_one_figure, figure_height, top_margin, bottom_margin, plot_height, vertical_plot_spacing, flag=first_fig_done)
                number_of_plots_in_one_figure = first_figure_plots_number
            else:
                number_of_rows, figure_height, top_margin, bottom_margin, subplots_adjust_top, subplots_adjust_bottom = orig_number_of_rows, orig_figure_height, orig_top_margin, orig_bottom_margin, orig_subplots_adjust_top, orig_subplots_adjust_bottom
                number_of_plots_in_one_figure = orig_number_of_plots_in_one_figure
            
            if remaining_figures < number_of_plots_in_one_figure and not already_adjusted:
                number_of_rows, figure_height, top_margin, bottom_margin, subplots_adjust_top, subplots_adjust_bottom, already_adjusted = adjust_height(number_of_rows, number_of_columns, remaining_figures, number_of_plots_in_one_figure, figure_height, top_margin, bottom_margin, plot_height, vertical_plot_spacing, flag=already_adjusted)

                '''#number_of_plots_in_one_figure = remaining_figures
                orig_number_of_rows = number_of_rows
                new_number_of_rows = int(remaining_figures / number_of_plots_in_one_figure * number_of_rows)
                if remaining_figures > new_number_of_rows * number_of_columns:
                    number_of_rows = new_number_of_rows + 1
                else:
                    number_of_rows = new_number_of_rows
                orig_height = figure_height
                figure_height = (top_margin+bottom_margin)*figure_height + number_of_rows * plot_height + vertical_plot_spacing * (number_of_rows - 1)#number_of_rows / orig_number_of_rows * (figure_height*1.2)
                top_margin = top_margin * orig_height / figure_height
                bottom_margin = bottom_margin * orig_height / figure_height
                subplots_adjust_top = 1 - top_margin
                subplots_adjust_bottom = bottom_margin
                already_adjusted = True'''
            #END OF NEW PART:


            fig = plt.figure(figsize=(figure_width * cm, figure_height *cm))
            fig.tight_layout()
            #fig, axs = plt.subplots(number_of_rows, number_of_columns, figsize=(figure_width * cm, figure_height *cm))
            if len(list_of_x_values) == 1:
                fig_file_name = f"{fig_file_basename}.png"
            else:
                fig_file_name = f"{fig_file_basename}_fig{fig_counter}.png"
            fig_file_path = os.path.join(fig_file_dir, fig_file_name)
            fig_counter += 1
            #axs=axs.flatten()
        
        # plotting the subplot
        #subplot_index = plot_counter % 35 + 1
        #subplot_index = plot_counter % number_of_plots_in_one_figure + 1
        subplot_index = counter % number_of_plots_in_one_figure + 1
        
        ax = fig.add_subplot(number_of_rows, number_of_columns, subplot_index)
        #ax=axs[subplot_index-1]
        if plot_colors is None:
            plot_color = plot_settings.get("color")
        else:
            plot_color = plot_colors[plot_counter]
        ax.plot(x_values, y_values, color = plot_color,#plot_settings.get("color"), 
            linewidth = plot_settings.get("linewidth"))
        #plot = plt.subplot(number_of_rows, number_of_columns, subplot_index)
        #plot.set_size_inches(2, 2)
        
        if list_of_scatter_x_values != None:
            scatter_x_values = list_of_scatter_x_values[plot_counter]
            if list_of_scatter_y_values != None:
                scatter_y_values = list_of_scatter_y_values[plot_counter]
            else:
                scatter_y_values = []
                for scatter_x in scatter_x_values:
                    if scatter_x == None:
                        scatter_y = None
                    else:
                        scatter_y = find_closest_value(
                            x_values, y_values, scatter_x)
                    scatter_y_values.append(scatter_y)
            ax.scatter(scatter_x_values, scatter_y_values, 
                color = plot_settings.get("scatter_color"), 
                alpha = plot_settings.get("scatter_transparency", 1),
                s = plot_settings.get("scatter_size", None))
        
        """if correct_shift != 0:
            plot_value = find_closest_value(x_values, y_values, correct_shift)
            ax.scatter(correct_shift, plot_value, color='#00ff0080')
        """
        
        # settings for title and axis labels
        if titles != None:
            ax.set_title(titles[plot_counter], fontsize = 7, y = 1.05)#, pad = 2)
        if xlabel != None:
            ax.set_xlabel(xlabel = xlabel, loc = 'right', fontsize = 6)
            ax.xaxis.set_label_coords(*plot_settings.get("xaxis_label_coords"))
        if ylabel != None:
            ax.set_ylabel(ylabel = ylabel, loc = 'top', fontsize = 6)
            ax.yaxis.set_label_coords(*plot_settings.get("yaxis_label_coords"))
        
        # setting axis global labels
        if xlabel_global != None:
            fig.supxlabel(xlabel_global, fontsize = 8)#, verticalalignment = "center")
        if ylabel_global != None:
            fig.supylabel(ylabel_global, fontsize = 8)#, horizontalalignment = "center")
        
        # setting axis limits and directions
        if xlim != None: ax.set_xlim(*xlim)
        if ylim != None: ax.set_ylim(*ylim)
        ax.xaxis.set_inverted(plot_settings.get("invert_xaxis", False))
        ax.yaxis.set_inverted(plot_settings.get("invert_yaxis", False))
        
        # setting ticks
        single_decimal_scientific_format = ScalarFormatterForceFormat()
        single_decimal_scientific_format.set_powerlimits((0,0))
        if xticks != None:
            ax.set_xticks(ticks = xticks, labels = xticks_labels)
        if yticks != None:
            ax.set_yticks(ticks = yticks, labels = yticks_labels)
        else:
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, steps = [1,2,4,10], min_n_ticks = 4))
            ax.yaxis.set_major_formatter(single_decimal_scientific_format)
        '''else:
            maximum = max(y_values)
            minimum = min(y_values)
            y_range = maximum - minimum
            y_range_str = str(y_range)
            decimal_pos = y_range_str.find(".")
            if decimal_pos == -1:
                
            precision = 
            increment = y_range / 5
            minimum + '''

        ax.tick_params(which = 'both', 
            labelsize = plot_settings.get("tick_labelsize", 6), 
            length = plot_settings.get("tick_length", 2), #3
            pad = plot_settings.get("tick_pad", 0), #2
            width = plot_settings.get("tick_width", 0.5),
            labelcolor = "black")
        
        #ax.yaxis.set_major_locator(ticker.MultipleLocator(20000))
        
        '''ax.ticklabel_format(axis='y', 
            style='scientific', 
            scilimits=(0,0),
            #useMathText=True
            )'''
        #ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.1f}"))
        #ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.1f}"))
        ax.yaxis.get_offset_text().set_fontsize(6)
        ax.yaxis.get_offset_text().set_position((0, 0.95))
        '''ax.yaxis.set_major_formatter(
            ticker.LogFormatterExponent(base=10.0, labelOnlyBase=True)
        )'''
        
        # setting border linewidth:
        for border in ['top', 'bottom', 'left', 'right']: 
            ax.spines[border].set_linewidth(plot_settings.get("border_linewidth"))    
        
        # if there are enough plots in one figure or
        # we reach the last plot
        '''if (plot_counter % number_of_plots_in_one_figure == 
            number_of_plots_in_one_figure - 1 
            or 
            plot_counter == len(list_of_x_values) - 1):'''
        if (counter % number_of_plots_in_one_figure == 
            number_of_plots_in_one_figure - 1 
            or 
            plot_counter == len(list_of_x_values) - 1
            or plot_counter == first_figure_plots_number - 1):

            plt.subplots_adjust(
                wspace = subplots_adjust_wspace,#plot_settings.get("subplots_adjust_wspace", 0.28), 
                hspace = subplots_adjust_hspace,#plot_settings.get("subplots_adjust_hspace", 0.35), 
                left = subplots_adjust_left,#plot_settings.get("subplots_adjust_left", 0.08), 
                right = subplots_adjust_right,#plot_settings.get("subplots_adjust_right", 0.95), 
                bottom = subplots_adjust_bottom,#plot_settings.get("subplots_adjust_bottom", 0.05), 
                top = subplots_adjust_top,#plot_settings.get("subplots_adjust_top", 0.8),#0.95)
                )

            fig.savefig(fig_file_path, dpi = 300)#, transparent = transparent)
            plt.close(fig)
            first_fig_done_fully = True
