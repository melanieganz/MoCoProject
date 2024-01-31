# Docstrings have been generated using ChatGPT to speed-up the book-keeping.

import os
import math
import data_loader
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import probplot, pearsonr, kstest, mannwhitneyu
from matplotlib.ticker import PercentFormatter, MultipleLocator
from collections import Counter


def plot_histogram_test(m: pd.DataFrame, f: pd.DataFrame, age: int, n_bins: int = 30):
    cols = ["max_trans_x_mm", "max_trans_y_mm", "max_trans_z_mm", "max_rot_x_deg", "max_rot_y_deg", "max_rot_z_deg"]
    m = m[cols]
    f = f[cols]

    max_x_trans = max(
        int(math.ceil(m.iloc[:, 0].max() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 1].max() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 2].max() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 0].max() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 1].max() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 2].max() / 10.0)) * 10
    )

    max_x_rot = max(
        int(math.ceil(m.iloc[:, 3].max() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 4].max() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 5].max() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 3].max() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 4].max() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 5].max() / 10.0)) * 10
    )

    min_x_trans = min(
        int(math.ceil(m.iloc[:, 0].min() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 1].min() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 2].min() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 0].min() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 1].min() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 2].min() / 10.0)) * 10
    )

    min_x_rot = min(
        int(math.ceil(m.iloc[:, 3].min() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 4].min() / 10.0)) * 10,
        int(math.ceil(m.iloc[:, 5].min() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 3].min() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 4].min() / 10.0)) * 10,
        int(math.ceil(f.iloc[:, 5].min() / 10.0)) * 10
    )

    max_density_trans = max(
        np.histogram(m.iloc[:,0], bins=n_bins, density=True)[0].max(),
        np.histogram(m.iloc[:,1], bins=n_bins, density=True)[0].max(),
        np.histogram(m.iloc[:,2], bins=n_bins, density=True)[0].max(),
        np.histogram(f.iloc[:,0], bins=n_bins, density=True)[0].max(),
        np.histogram(f.iloc[:,1], bins=n_bins, density=True)[0].max(),
        np.histogram(f.iloc[:,2], bins=n_bins, density=True)[0].max(),
    )

    max_density_rot = max(
        np.histogram(m.iloc[:,3], bins=n_bins, density=True)[0].max(),
        np.histogram(m.iloc[:,4], bins=n_bins, density=True)[0].max(),
        np.histogram(m.iloc[:,5], bins=n_bins, density=True)[0].max(),
        np.histogram(f.iloc[:,3], bins=n_bins, density=True)[0].max(),
        np.histogram(f.iloc[:,4], bins=n_bins, density=True)[0].max(),
        np.histogram(f.iloc[:,5], bins=n_bins, density=True)[0].max(),
    )

    fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(15, 10))


    def plot_hist(ax, data, title, bins, color, x_label, max_x, min_x, max_density):
        if max_x % 20 == 10:
            max_x += 10
        elif min_x % 20 == 10:
            min_x -= 10

        mu = np.mean(data)
        sigma = np.std(data)

        norm_sample = np.random.normal(mu, sigma, len(data))

        # Plotting histogram
        #ax.hist(data, bins="auto", color=color, density=True)
        sns.histplot(data, bins="auto", stat='density', color=color, kde=False, ax=ax, label="actual")
        #sns.histplot(norm_sample, bins="auto", stat='density', color="cyan", kde=False, ax=ax, label="theoretical")
        #ax.hist(norm_sample, bins="auto", color="cyan", density=True)

        ax.set_xlabel(x_label)
        ax.set_ylabel("Density")
        ax.yaxis.grid(True, linestyle='--', alpha=0.5)
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(0, (max_density + 0.05))
        ax.set_xticks(np.arange(min_x, max_x, 20))
        ax.set_yticks(np.arange(0, max_density + 0.1, 0.1))
        ax.set_title(title)
        ax.legend()

    plot_hist(axs[0][0], m.iloc[:,0], f"Translation along the X-axis for {age} years old males.", n_bins, "red", "Translation (mm)", max_x_trans, min_x_trans, max_density_trans)
    plot_hist(axs[0][1], m.iloc[:,1], f"Translation along the Y-axis for {age} years old males.", n_bins, "green", "Translation (mm)", max_x_trans, min_x_trans, max_density_trans)
    plot_hist(axs[0][2], m.iloc[:,2], f"Translation along the Z-axis for {age} years old males.", n_bins, "blue", "Translation (mm)", max_x_trans, min_x_trans, max_density_trans)
    plot_hist(axs[1][0], f.iloc[:,0], f"Translation along the X-axis for {age} years old females.", n_bins, "red", "Translation (mm)", max_x_trans, min_x_trans, max_density_trans)
    plot_hist(axs[1][1], f.iloc[:,1], f"Translation along the Y-axis for {age} years old females.", n_bins, "green", "Translation (mm)", max_x_trans, min_x_trans, max_density_trans)
    plot_hist(axs[1][2], f.iloc[:,2], f"Translation along the Z-axis for {age} years old females.", n_bins, "blue", "Translation (mm)", max_x_trans, min_x_trans, max_density_trans)
    plot_hist(axs[2][0], m.iloc[:,3], f"Rotation around the X-axis for {age} years old males.", n_bins, "purple", "Rotation (deg)", max_x_rot, min_x_rot, max_density_rot)
    plot_hist(axs[2][1], m.iloc[:,4], f"Rotation around the Y-axis for {age} years old males.", n_bins, "orange", "Rotation (deg)", max_x_rot, min_x_rot, max_density_rot)
    plot_hist(axs[2][2], m.iloc[:,5], f"Rotation around the Z-axis for {age} years old males.", n_bins, "teal", "Rotation (deg)", max_x_rot, min_x_rot, max_density_rot)
    plot_hist(axs[3][0], f.iloc[:,3], f"Rotation around the X-axis for {age} years old females.", n_bins, "purple", "Rotation (deg)", max_x_rot, min_x_rot, max_density_rot)
    plot_hist(axs[3][1], f.iloc[:,4], f"Rotation around the Y-axis for {age} years old females.", n_bins, "orange", "Rotation (deg)", max_x_rot, min_x_rot, max_density_rot)
    plot_hist(axs[3][2], f.iloc[:,5], f"Rotation around the Z-axis for {age} years old females.", n_bins, "teal", "Rotation (deg)", max_x_rot, min_x_rot, max_density_rot)

    plt.tight_layout()
    plt.show()


def plot_boxplot(m: np.ndarray, f: np.ndarray, age: int, motion: str, modality: str, n_bins: int = 20) -> None:
    """
    Plots boxplots for translation and rotation data.

    Parameters:
    - m (numpy.ndarray): The male data array with shape (N, 6), where N is the number of data points.
    - f (numpy.ndarray): The female data array with shape (N, 6), where N is the number of data points.
    - age (int): The age for the data.
    - motion (str): The type of motion, either 'translation' or 'rotation'.
    - modality (str): The modality of the data, either 'fmri' or 'dti'.
    - n_bins (int, optional): The number of bins for the histogram. Default is 20.

    Returns:
    None

    This function generates a boxplot comparing translation or rotation data between males (m) and females (f).
    The left side of the plot represents female data, and the right side represents male data.

    The color-coded labels indicate different dimensions (X, Y, Z), and the y-axis is presented as a density
    to account for variations in the number of data points.

    If the modality is 'fmri', the x-axis ticks of the plot will be set at every 5 units; if 'dti', every 2 units.

    Example:
    >>> plot_boxplot(male_data, female_data, age=10, motion='translation', modality='fmri', n_bins=25)
    """
    labels = ['f-Z', 'm-Z', 'f-Y', 'm-Y', 'f-X', 'm-X']

    if motion == "translation":
        colors  = ['blue', 'blue', 'green', 'green', 'red', 'red']
        mode    = "trans"
        measure = "mm"
    elif motion == "rotation":
        colors  = ['teal', 'teal', 'orange', 'orange', 'purple', 'purple']
        mode    = "rot"
        measure = "deg"

    fig, axs = plt.subplots(figsize=(14, 10))

    bp = axs.boxplot([f[f'max_{mode}_z_{measure}'], m[f'max_{mode}_z_{measure}'],
                        f[f'max_{mode}_y_{measure}'], m[f'max_{mode}_y_{measure}'], 
                        f[f'max_{mode}_x_{measure}'], m[f'max_{mode}_x_{measure}']],
                        vert=False,  # box alignment
                        patch_artist=True,  # fill with color
                        medianprops=dict(color='black'),
                        flierprops=dict(marker='.', markerfacecolor='black', markersize=3),
                        showmeans=True,
                        meanprops=dict(marker='x', markeredgecolor='black', markerfacecolor='black'),
                        labels=labels)  # will be used to label x-ticks
    
    #axs.set_title(f"Comparisson of {motion} between {age} years old males and females.")

    if modality == 'fmri': ticks = 10    
    if modality == 'dti': ticks = 2
        
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    for ax in [axs]:
        ax.xaxis.grid(True)
        ax.set_xlabel(f'{motion.capitalize()} ({measure})')
        ax.set_ylabel("Sex-axis")

    plt.gca().xaxis.set_major_locator(MultipleLocator(ticks))

    def count_outliers(bplot):
        for line in bplot['fliers']:
            x_data = line.get_xdata()
            outliers = len(x_data)
            print(f"Number of outliers: {outliers}")

    count_outliers(bp)

    plt.tight_layout()
    plt.show()


def generate_custom_qq_plot(ax: plt.Axes, data: np.ndarray, title: str, line_color_1: str = 'blue', line_color_2: str = 'black') -> None:
    """
    Generates a customized Q-Q plot with optional line and dot color customization.

    Parameters:
    - ax (plt.Axes): The matplotlib axes where the Q-Q plot will be drawn.
    - data (numpy.ndarray): The input data array for the Q-Q plot.
    - title (str): The title of the Q-Q plot.
    - line_color_1 (str, optional): The color of the straight line on the Q-Q plot. Default is 'blue'.
    - line_color_2 (str, optional): The color of the dots on the Q-Q plot. Default is 'black'.

    Returns:
    None

    This function uses the scipy.stats.probplot to generate a Q-Q plot on the specified axes. It allows customization
    of the color of the straight line and dots on the plot.

    Example:
    >>> generate_custom_qq_plot(ax, data[:, 0], 'X-translation (mm)', line_color_1='red', line_color_2='black')
    """

    probplot(data, dist='norm', plot=ax)

    # Customize the straight line color
    line = ax.get_lines()[0]
    line.set_color(line_color_1)

    # Customize the dot color
    scatter = ax.get_lines()[1]
    scatter.set_color(line_color_2)
    
    ax.set_title(title)
    ax.set_xlabel('Theoretical Quantiles')
    ax.set_ylabel('Sample Quantiles')


def plot_qq_plot(data: np.ndarray, sex: str, age: int) -> None:
    """
    Plots a 2x3 grid of Q-Q plots for each dimension of translation and rotation data.

    Parameters:
    - data (numpy.ndarray): The input data array with shape (N, 6), where N is the number of data points.
    - sex (str): The sex of the individuals, either 'male' or 'female'.
    - age (int): The age of the individuals.

    Returns:
    None

    This function generates a 2x3 grid of Q-Q plots, with the first row showing Q-Q plots for translation data
    (X, Y, Z), and the second row showing Q-Q plots for rotation data (rotX, rotY, rotZ). Each Q-Q plot includes
    a straight line and dots representing the theoretical and sample quantiles, respectively.

    Example:
    >>> plot_qq_plot(data, sex='male', age=25)
    """

    # Create a 2x3 grid of subplots
    _, axes = plt.subplots(2, 3, figsize=(15, 8))

    # Generate customized Q-Q plots for the first line
    generate_custom_qq_plot(axes[0, 0], data['max_trans_x_mm'], f'X-translation (mm) for {age} years old {sex}.', line_color_1='red', line_color_2='black')
    generate_custom_qq_plot(axes[0, 1], data['max_trans_y_mm'], f'Y-translation (mm) for {age} years old {sex}.', line_color_1='green', line_color_2='black')
    generate_custom_qq_plot(axes[0, 2], data['max_trans_z_mm'], f'Z-translation (mm) for {age} years old {sex}.', line_color_1='blue', line_color_2='black')

    # Generate customized Q-Q plots for the second line
    generate_custom_qq_plot(axes[1, 0], data['max_rot_x_deg'], f'X-rotation (deg) for {age} years old {sex}.', line_color_1='purple', line_color_2='black')
    generate_custom_qq_plot(axes[1, 1], data['max_rot_y_deg'], f'Y-rotation (deg) for {age} years old {sex}.', line_color_1='orange', line_color_2='black')
    generate_custom_qq_plot(axes[1, 2], data['max_rot_z_deg'], f'Z-rotation (deg) for {age} years old {sex}.', line_color_1='teal', line_color_2='black')

    plt.tight_layout()
    plt.show()


def plot_correlation(data: np.ndarray) -> None:
    """
    Plot correlation maps and leave significant correlations blank.

    Parameters:
    - data (numpy.ndarray): The input data array with shape (N, 6), where N is the number of data points.
    - age (int): The age of the individuals.
    - sex (str): The sex of the individuals, either 'male' or 'female'.
    - modality (str): The modality of the data, either 'fmri' or 'dti'.

    Returns:
    None

    This function calculates the Pearson correlation coefficient and significance for pairs of variables in the
    input data. It then plots the correlation map, leaving significant correlations blank.

    Example:
    >>> plot_correlation(data, age=25, sex='male', modality='fmri')
    """
    c_val = []
    p_val = []

    cols = ["max_trans_x_mm", "max_trans_y_mm", "max_trans_z_mm", "max_rot_x_deg", "max_rot_y_deg", "max_rot_z_deg"]
    data = data[cols]

    for i in range(6):
        x = data.iloc[:, i]
        for j in range(6):
            y = data.iloc[:, j]
            c, p = pearsonr(x, y)
            c_val.append(np.round(c, 3))
            p_val.append(np.round(p, 3))

    c_val = np.array(c_val)
    p_val = np.array(p_val)

    c_val = np.reshape(c_val, (6, 6))
    p_val = np.reshape(p_val, (6, 6))

    fig, axs = plt.subplots(figsize=(12, 6))

    # Replace significant correlations with NaN
    c_val[p_val >= 0.05] = np.nan

    im1 = axs.imshow(c_val, interpolation='nearest', cmap='viridis', vmin=-1, vmax=1)  # Adjust vmin and vmax as needed
    #axs.set_title(f'{age} years old {sex} motion correlation map.')
    axs.set_xticks([0, 1, 2, 3, 4, 5])
    axs.set_yticks([0, 1, 2, 3, 4, 5])
    axs.set_xticklabels(['trans-X', 'trans-Y', 'trans-Z', 'rot-X', 'rot-Y', 'rot-Z'])
    axs.set_yticklabels(['trans-X', 'trans-Y', 'trans-Z', 'rot-X', 'rot-Y', 'rot-Z'])

    # Add colorbars
    cbar1 = fig.colorbar(im1, ax=axs)

    plt.tight_layout()
    plt.show()


def kolmogorov_smirnov(data: np.ndarray) -> tuple:
    """
    Perform Kolmogorov-Smirnov tests for normality on each variable in the input data.

    Parameters:
    - data (numpy.ndarray): The input data array with shape (N, 6), where N is the number of data points.

    Returns:
    tuple: A tuple containing two numpy arrays - k_val and p_val, representing the KS test statistics and p-values.

    This function performs Kolmogorov-Smirnov tests for normality on each variable in the input data and returns
    the KS test statistics and p-values.

    Example:
    >>> k_vals, p_vals = kolmogorov_smirnov(data)
    """
    k_val = []
    p_val = []

    cols = ["max_trans_x_mm", "max_trans_y_mm", "max_trans_z_mm", "max_rot_x_deg", "max_rot_y_deg", "max_rot_z_deg"]
    data = data[cols]

    for i in range(6):
        rvs = data.iloc[:, i]

        #print(rvs.shape)

        mu    = np.mean(rvs, axis=0)
        sigma = np.std(rvs, axis=0)

        # Generate a sample from the normal distribution with the same size as your data
        #norm_sample = np.random.normal(mu, sigma, len(rvs))
        rvs = (rvs - mu) / sigma

        # Perform the Kolmogorov-Smirnov test
        k, p = kstest(rvs, "norm")
        k_val.append(k)
        p_val.append(p)

    return np.array(k_val), np.array(p_val)


def plot_2dhist(data: pd.DataFrame) -> None:
    """
    Plot a 2D histogram for the given data.

    Parameters:
    - data (numpy.ndarray): The input data array with shape (N, M), where N is the number of data points and
                           M is the number of variables (at least 4 for the specified indices).

    Returns:
    None

    This function plots a 2D histogram for the specified columns (1 and 3) of the input data.

    Example:
    >>> plot_2dhist(data)
    """
    cols = ["max_trans_x_mm", "max_trans_y_mm", "max_trans_z_mm", "max_rot_x_deg", "max_rot_y_deg", "max_rot_z_deg"]
    data = data[cols]

    plt.hist2d(data.iloc[:, 1], data.iloc[:, 2], bins=225, range=[[-15, +15], [-15, +15]], density=True)

    plt.xlabel('trans-Y')
    plt.ylabel('rot-X')
    #plt.title('2D Histogram')

    cbar = plt.colorbar()
    cbar.set_label('Frequency')

    plt.show()


def plot_pie_test(data: pd.DataFrame):
    cols = ["runs"]
    data = data[cols]

    fig, axs = plt.subplots(figsize=(15, 10))

    def plot_pie_chart(ax, data, title, legend_title):
        values = list(Counter(data).values())
        labels = list(Counter(data).keys())

        # Sort both lists based on values in descending order
        sorted_indices = sorted(range(len(values)), key=lambda k: values[k], reverse=True)
        sorted_labels = [labels[i] for i in sorted_indices]
        sorted_values = [values[i] for i in sorted_indices]

        set1_cmap = plt.get_cmap("Set1", 9)

        # Extract the first 4 colors and reverse their order
        first_4_colors_flipped = set1_cmap.colors

        # Define a color map using the modified colors
        color_map = {label: first_4_colors_flipped[i] for i, label in enumerate(set(sorted_labels))}

        # Map colors based on the defined color map
        colors = [color_map[sorted_label] for sorted_label in sorted_labels]

        percentages = [(value / sum(sorted_values)) * 100 for value in sorted_values]
        percentages = [f'{v:.1f}%\n{l}' for v, l in zip(percentages, sorted_labels)]  # Add percentages to labels
        wedges, texts, autotexts = ax.pie(sorted_values, autopct='', startangle=90, colors=colors)
        ax.legend(wedges, percentages, title=legend_title, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1))
        #ax.set_title(title)

    plot_pie_chart(axs, data.iloc[:, 0], "Number of runs", "%-run")

    plt.tight_layout(rect=[0, 0, 0.9, 0.96], h_pad=2)
    plt.show()


def plot_pie(m: pd.DataFrame, f: pd.DataFrame, age: int, motion: str):
    cols = ["run_max_trans_x_mm", "run_max_trans_y_mm", "run_max_trans_z_mm", "run_max_rot_x_deg", "run_max_rot_y_deg", "run_max_rot_z_deg"]
    m = m[cols]
    f = f[cols]

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))

    def plot_pie_chart(ax, data, title, legend_title):
        values = list(Counter(data).values())
        labels = list(Counter(data).keys())

        # Sort both lists based on values in descending order
        sorted_indices = sorted(range(len(values)), key=lambda k: values[k], reverse=True)
        sorted_labels = [labels[i] for i in sorted_indices]
        sorted_values = [values[i] for i in sorted_indices]

        set1_cmap = plt.get_cmap("Set1", 9)

        # Extract the first 4 colors and reverse their order
        first_4_colors = set1_cmap.colors[:4]

        # Define a color map using the modified colors
        color_map = {label: first_4_colors[i % 4] for i, label in enumerate(set(sorted_labels))}

        # Map colors based on the defined color map
        colors = [color_map[sorted_label] for sorted_label in sorted_labels]

        percentages = [(value / sum(sorted_values)) * 100 for value in sorted_values]
        percentages = [f'{v:.1f}%\n{l}' for v, l in zip(percentages, sorted_labels)]  # Add percentages to labels
        wedges, texts, autotexts = ax.pie(sorted_values, autopct='', startangle=90, colors=colors)
        ax.legend(wedges, percentages, title=legend_title, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1))
        ax.set_title(title)


    if motion == "translation":
        plot_pie_chart(axs[0, 0], m.iloc[:, 0], "X-translation (mm) - Males", "%-run")
        plot_pie_chart(axs[0, 1], m.iloc[:, 1], "Y-translation (mm) - Males", "%-run")
        plot_pie_chart(axs[0, 2], m.iloc[:, 2], "Z-translation (mm) - Males", "%-run")
        plot_pie_chart(axs[1, 0], f.iloc[:, 0], "X-translation (mm) - Females", "%-run")
        plot_pie_chart(axs[1, 1], f.iloc[:, 1], "Y-translation (mm) - Females", "%-run")
        plot_pie_chart(axs[1, 2], f.iloc[:, 2], "Z-translation (mm) - Females", "%-run")
    elif motion == "rotation":
        plot_pie_chart(axs[0, 0], m.iloc[:, 3], "X-rotation (deg) - Males", "%-run")
        plot_pie_chart(axs[0, 1], m.iloc[:, 4], "Y-rotation (deg) - Males", "%-run")
        plot_pie_chart(axs[0, 2], m.iloc[:, 5], "Z-rotation (deg) - Males", "%-run")
        plot_pie_chart(axs[1, 0], f.iloc[:, 3], "X-rotation (deg) - Females", "%-run")
        plot_pie_chart(axs[1, 1], f.iloc[:, 4], "Y-rotation (deg) - Females", "%-run")
        plot_pie_chart(axs[1, 2], f.iloc[:, 5], "Z-rotation (deg) - Females", "%-run")

    #plt.suptitle(f'Runs where highest displacement occurred for {age} years old males and females.')
    plt.tight_layout(rect=[0, 0, 0.9, 0.96], h_pad=2)
    plt.show()


def u_test(x: pd.DataFrame, y: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Perform independent two-sample t-tests for each variable between two groups.

    Parameters:
    - x (numpy.ndarray): Data array for the first group with shape (N, M), where N is the number of data points
                         and M is the number of variables.
    - y (numpy.ndarray): Data array for the second group with the same shape (N, M).

    Returns:
    - tuple[np.ndarray, np.ndarray]: A tuple containing arrays of t-statistics and corresponding p-values for
                                     each variable.

    This function conducts independent two-sample t-tests for each variable between the two groups.

    Example:
    >>> t_vals, p_vals = t_test(group1_data, group2_data)
    """

    u_vals = []
    p_vals = []

    cols = ["max_trans_x_mm", "max_trans_y_mm", "max_trans_z_mm", "max_rot_x_deg", "max_rot_y_deg", "max_rot_z_deg"]
    x = x[cols]
    y = y[cols]

    for i in range(x.shape[1]):
        variable_x = x.iloc[:, i]
        variable_y = y.iloc[:, i]

        #print(variable_x)

        u_stat, p_value = mannwhitneyu(variable_x, variable_y)
        u_vals.append(np.round(u_stat, 4))
        p_vals.append(np.round(p_value, 4))

    return np.array(u_vals), np.array(p_vals)


def plot_u_test(u_val: np.ndarray, p_val: np.ndarray, age: int):
    """
    Plot a table showing t-statistics and p-values for each variable between two groups.

    Parameters:
    - t_val (numpy.ndarray): Array of t-statistics with shape (M,), where M is the number of variables.
    - p_val (numpy.ndarray): Array of p-values with the same shape (M,).

    This function creates a table showing t-statistics and p-values for each variable between two groups.

    Example:
    >>> plot_t_test(t_values, p_values)
    """
    df = pd.DataFrame({'u-val': u_val, 'p-val': p_val},
                      index=["trans-X", "trans-Y", "trans-Z", "rot-X", "rot-Y", "rot-Z"])

    fig, ax = plt.subplots()

    table = ax.table(cellText=df.reset_index().values,
                     colLabels=df.columns.insert(0, ''),
                     cellLoc='center',
                     loc='upper center')

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.2)

    ax.set_title(f"Statistical significance of translation along/rotation around each axis for {age} years old males and females.",
                 fontweight="bold", fontsize=12)
    ax.axis('off')

    plt.show()


def plot_fd(male, female, age, modality):
    import os
    m_i, _ = np.random.randint(0, male.shape)
    f_i, _ = np.random.randint(0, female.shape)
    patients = [male['participant_id'].iloc[m_i], female['participant_id'].iloc[f_i]]
    run_nr = 1
    subj_path = os.path.join("data", modality, "resting_state", "fmriresults01", "derivatives", "abcd-hcp-pipeline")
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    for i, patient in enumerate(patients):
        if i == 0:
            sex = "male"
        else:
            sex = "female"

        if modality == 'dti':
            path = os.path.join(subj_path, patient + "_ses-baselineYear1Arm1_confounds.tsv")
        if modality == 'fmri':
            path = os.path.join(subj_path, patient, "ses-baselineYear1Arm1", "func", patient + "_ses-baselineYear1Arm1_task-rest_run-" + str(run_nr) + "_desc-includingFD_motion.tsv")
        if modality == 'nback':
            path = os.path.join(subj_path, patient, "ses-baselineYear1Arm1", "func", patient + "_ses-baselineYear1Arm1_task-nback_run-" + str(run_nr) + "_desc-includingFD_motion.tsv")
        displacements = data_loader.load_data(path, modality)
        axes[i, 0].plot(displacements[:, 0], label="trans-X", color="red")
        axes[i, 0].plot(displacements[:, 1], label="trans-Y", color="green")
        axes[i, 0].plot(displacements[:, 2], label="trans-Z", color="blue")
        axes[i, 0].set_ylim(-2, 2)
        axes[i, 0].legend()
        #axes[i, 0].set_title(f"Random {age}-year-old {sex} translation motion curve.")
        axes[i, 0].set_xlabel("Frame number")
        axes[i, 0].set_ylabel("Translation (mm)")
        if modality == 'dti':
            axes[i, 1].plot(np.rad2deg(displacements[:, 3]), label="rot-X", color="purple")
            axes[i, 1].plot(np.rad2deg(displacements[:, 4]), label="rot-Y", color="orange")
            axes[i, 1].plot(np.rad2deg(displacements[:, 5]), label="rot-Z", color="teal")
        if modality == 'fmri' or modality == 'nback':
            axes[i, 1].plot(displacements[:, 3], label="rot-X", color="purple")
            axes[i, 1].plot(displacements[:, 4], label="rot-Y", color="orange")
            axes[i, 1].plot(displacements[:, 5], label="rot-Z", color="teal")
        axes[i, 1].set_ylim(-2, 2)
        axes[i, 1].legend()
        #axes[i, 1].set_title(f"Random {age}-year-old {sex} rotation motion curve.")
        axes[i, 1].set_xlabel("Frame number")
        axes[i, 1].set_ylabel("Rotation (deg)")
    plt.tight_layout()
    plt.show()

def plot_fd_test(data, modality):
    d_i       = np.random.randint(0, data.shape[0])
    run_nr    = 1
    
    if modality == "fmri":
        subj_path = os.path.join("data", modality, "resting_state", "fmriresults01", "derivatives", "abcd-hcp-pipeline")
    if modality == "dti":
        subj_path = os.path.join("data", modality)
    
    patient   = data['participant_id'].iloc[d_i]

    if modality == 'dti':
            path = os.path.join(subj_path, patient + "_ses-baselineYear1Arm1_confounds.tsv")
    if modality == 'fmri':
            path = os.path.join(subj_path, patient, "ses-baselineYear1Arm1", "func", patient + "_ses-baselineYear1Arm1_task-rest_run-" + str(run_nr) + "_desc-includingFD_motion.tsv")
    if modality == 'nback':
            path = os.path.join(subj_path, patient, "ses-baselineYear1Arm1", "func", patient + "_ses-baselineYear1Arm1_task-nback_run-" + str(run_nr) + "_desc-includingFD_motion.tsv")
    
    displacements = data_loader.load_data(path, modality)

    fig, axes = plt.subplots(1, 2, figsize=(15, 10))

    axes[0].plot(displacements[:, 0], label="trans-X", color="red")
    axes[0].plot(displacements[:, 1], label="trans-Y", color="green")
    axes[0].plot(displacements[:, 2], label="trans-Z", color="blue")
    axes[0].set_ylim(-2,2)
    axes[0].legend()
    axes[0].set_xlabel("Frame number")
    axes[0].set_ylabel("Translation (mm)")

    axes[1].plot(np.rad2deg(displacements[:, 3]), label="rot-X", color="purple")
    axes[1].plot(np.rad2deg(displacements[:, 4]), label="rot-Y", color="orange")
    axes[1].plot(np.rad2deg(displacements[:, 5]), label="rot-Z", color="teal")
    axes[1].set_ylim(-2,2)
    axes[1].legend()
    axes[1].set_xlabel("Frame number")
    axes[1].set_ylabel("Rotation (deg)")

    plt.tight_layout()
    plt.show()

modality = input("Enter modality: ")
data = data_loader.read_data(modality)
n_rows, n_cols = data.shape

#print(data)
#plot_pie_test(data)
plot_2dhist(data)

group_age = data.groupby('age')
for age in group_age.groups.keys():
    x = group_age.get_group(age)
    group_sex = x.groupby('sex')
    m = group_sex.get_group(1)
    f = group_sex.get_group(2)
    m_rows, _ = m.shape
    f_rows, _ = f.shape

    plot_qq_plot(m, "males", age)
    plot_qq_plot(f, "females", age)

    #plot_fd(m, f, age, modality)

    #plot_2dhist(f)
    #plot_correlation(m)

    #plot_pie(m, f, age, "rotation")

    """print("---m------------")
    print("m-t-x var: " + str(round(m['max_trans_x_mm'].std(),3)))
    print("m-t-y var: " + str(round(m['max_trans_y_mm'].std(),3)))
    print("m-t-z var: " + str(round(m['max_trans_z_mm'].std(),3)))
    print("m-r-x var: " + str(round(m['max_rot_x_deg'].std(),3)))
    print("m-r-y var: " + str(round(m['max_rot_y_deg'].std(),3)))
    print("m-r-z var: " + str(round(m['max_rot_z_deg'].std(),3)))

    print("---f------------")
    print("f-t-x var: " + str(round(f['max_trans_x_mm'].std(),3)))
    print("f-t-y var: " + str(round(f['max_trans_y_mm'].std(),3)))
    print("f-t-z var: " + str(round(f['max_trans_z_mm'].std(),3)))
    print("f-r-x var: " + str(round(f['max_rot_x_deg'].std(),3)))
    print("f-r-y var: " + str(round(f['max_rot_y_deg'].std(),3)))
    print("f-r-z var: " + str(round(f['max_rot_z_deg'].std(),3)))"""

    """u_vals, p_vals = u_test(m, f)
    plot_u_test(u_vals, p_vals, age)"""

    """m = m[(m['runs'] >= 4) & (m['runs'] <= 9)]
    f = f[(f['runs'] >= 4) & (f['runs'] <= 9)]

    plot_pie(m, f, age, "translation")"""