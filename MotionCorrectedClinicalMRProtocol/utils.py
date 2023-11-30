import numpy as np
import matplotlib.pyplot as plt
import os



def MakeBoxplot(metric, colors):
    '''
    Creating a box plot with preset colors

    Parameters
    ----------
    metric : array
        values to plot.
    colors : array or list
        colors for the different types.

    Returns
    -------
    None.

    '''
    
    small = dict(markersize=3)
    box1 = plt.boxplot(metric, flierprops=small, widths=0.5)
    for patch, patch2, color in zip(box1['boxes'], box1['medians'], colors):
            patch.set(color=color)
            patch2.set(color=color)


def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.03, fs=None, maxasterix=None, col='dimgrey'):
    """
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    :param col: color of the asterixes or text. Optional. The default is dimgrey.
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    print('center: ', center)
    print('height: ', height)
    print('num1: ', num1)
    print('num2: ', num2)
    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh-0.035*barh)

    plt.plot(barx, bary, c='dimgrey')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs, c=col)


def Show_Stars(p_cor, ind, bars, heights, arange_dh=False, dh=.1, flexible_dh=False, col='dimgrey'):
    '''
    Function to show brackets with asterisks indicating statistical 
    significance

    Parameters
    ----------
    p_cor : array or list
        corrected p-values.
    ind : list
        indices for the comparisons corresponding to the individual p-values.
    bars : list or array
        x position of boxplots.
    heights : list or array
        maximal value visualised in boxplots.
    arange_dh : float, optional
        offset above heights. The default is False.
    col : str
        color of the asterixes. Optional, the default is dimgrey.

    Returns
    -------
    int
        returns 0, when completed..

    '''
    
    star = p_cor < 0.05
    starstar = p_cor < 0.001
    
    all_nrs = []
    #dh = .1
    if arange_dh==True:
        dhs = np.array([.1,.1,.2, .1, .1, .2, .1, .1, .2,.2, .2, .2, .2, .1, .1, .1, .1, .2, .3, .1, .2, .3, .2, .2])
    if arange_dh == 'REAC':
        dhs = np.array([.1, .1, .1, .2, .4, .4, .1, .1, .2, .4,.4])
    if arange_dh == 'OLD':
        dhs = np.array([.1, .1, .1, .2, .1, .1, .2 ])
    if arange_dh == '2D':
        dhs = np.array([.1, .1, .1, .2, .3, .3, .2 ])
    if arange_dh == 'diff':
        dhs = np.array([.1, .1, .3])
    if arange_dh == 'PAPER':
        dhs = np.array([.1, .1, .15, .25, .25, .1, .1, .15, .25, .25])
    fl_dh = .03
    print('p_cor: ', p_cor)
    print('len(p_cor): ', len(p_cor))
    print('ind: ', ind)
    for i in range(0,len(p_cor)):
        nr = ind[i]
        all_nrs.append(nr)
        if arange_dh:
            dh = dhs[i]
        if flexible_dh:
            dh = fl_dh
        if starstar[i]==True:
            barplot_annotate_brackets(nr[0], nr[1], '**', bars, heights, dh=dh, col=col)
            fl_dh += 0.05
        else:
            if star[i]==True:
                barplot_annotate_brackets(nr[0], nr[1], '*', bars, heights, dh=dh, col=col)
                fl_dh += 0.05

    return 0


def SortFiles(file, no_check=False):
    '''
    sort a number of files after date and pick the most recent file.

    Parameters
    ----------
    file : list
        list of filenames.

    Returns
    -------
    file : string
        filename corresponding to the most recent file.

    '''
    
    ind = []
    for i in range(len(file)):
        base = os.path.basename(file[i])
        if no_check:
            ind.append(i)
        else:
            if len(base)<30:
                ind.append(i)
        
    file = np.array(file)
    file = file[ind]
    file = np.sort(file)[::-1]
    return file


def Scale_Axis_Lim(data):
    '''
    Adjust the limits of the y axis.

    Parameters
    ----------
    data : array
        values to be plotted.

    Returns
    -------
    y_min : float
        minimal value for y-axis.
    y_max : float
        maixmal value for y-axis.

    '''
    
    maxim = np.amax(data)
    minim = np.amin(data)
    if np.isinf(maxim):
        temp = np.sort(data)
        maxim = temp[-2]

    y_min = minim - 0.1*(maxim-minim)
    y_max = maxim + 0.05*(maxim-minim)

    if np.isnan(y_min):
        y_min = 0

    if np.isnan(y_max):
        y_max = 1

    return y_min, y_max


def DrawLines(a, b, c, d, metric, lw=0.5, col='gray'):
    '''
    Function to connect single paired data in a boxplot with thin grey lines

    Parameters
    ----------
    a : list
        Defines from which column of the metric array to choose the y-values 
        from which the lines start.
    b : list
        Defines from which column of the metric array to choose the y-values 
        at which the lines end.
    c : list
        Defines the x position of the start of the lines.
    d : list
        Defines the x position of the end of the lines.
    metric : array
        metric values which are used for the boxplot and should be connected.
    lw : float
        Linewidth of the lines. Optional, the default is 0.5.
    col : string
        Color of the lines. Optional, the default is 'gray'.

    Returns
    -------
    0

    '''
    counts = {}
    
    for A,B,C,D in zip(a,b,c,d):
        for y1, y2 in zip(metric[:,A], metric[:,B]):
            key = str(y1)+','+str(y2)
            if key in counts:
                counts[key] += 1
            else:
                counts[key] = 1
                
            linewidth = lw*(1+counts[key]/2)
            plt.plot([C,D], [y1, y2], col, lw=linewidth)

    return 0


def DrawLines2(a, b, c, d, metric, lw=0.5, col='gray'):
    '''
    Function to connect single paired data in a boxplot with thin grey lines,
    works if data in metric not as array but as list of arrays

    Parameters
    ----------
    a : list
        Defines from which column of the metric array to choose the y-values 
        from which the lines start.
    b : list
        Defines from which column of the metric array to choose the y-values 
        at which the lines end.
    c : list
        Defines the x position of the start of the lines.
    d : list
        Defines the x position of the end of the lines.
    metric : list of arrays
        metric values which are used for the boxplot and should be connected.
    lw : float
        Linewidth of the lines. Optional, the default is 0.5.
     col : string
        Color of the lines. Optional, the default is 'gray'.

    Returns
    -------
    0

    '''
    
    counts = {}
    
    for A,B,C,D in zip(a,b,c,d):
        #print('zip(metric[A], metric[B]) : ', zip(metric[A], metric[B]))
        for y1, y2 in zip(metric[A], metric[B]):
            print('y1 : ', y1, 'y2 : ', y2)
            key = str(y1)+','+str(y2)
            #print('counts : ', counts)
            if key in counts:
                counts[key] += 1
            else:
                counts[key] = 1
                
            linewidth = lw*(1+counts[key]/2)
            plt.plot([C,D], [y1, y2], col, lw=linewidth)
            

    return 0