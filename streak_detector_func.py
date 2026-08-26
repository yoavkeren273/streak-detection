import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import astropy
from astropy.io import fits
from astropy.table import Table
import time

def make_point_list(df):
    df=df.loc[:,['x','y']]
    data_array = df.to_numpy()
    # make a list of (x,y) for x and y in data['x','y']
    points= [(float(x), float(y)) for x, y in data_array[:, 0:2]]
    return points

def estimate_oriented_ccd_bounds(points): #point = list of tuples (x,y)
    """Estimate CCD bounds in sky coordinates using the minimum-area oriented bounding box."""
    if len(points) == 0: 
        return None
    if len(points) <= 2: #draws a squre based on the diagonal, returns list of four tuples representing the corners
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)
        return [(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)]

    def convex_hull(pts):
        pts = sorted(pts)

        def cross(o, a, b): #cross product
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

        lower = []
        for p in pts:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        upper = []
        for p in reversed(pts):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
                upper.pop()
            upper.append(p)
        return lower[:-1] + upper[:-1]

    hull = convex_hull(points)
    if len(hull) <= 2:
        xs = [p[0] for p in hull]
        ys = [p[1] for p in hull]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)
        return [(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)]

    best_area = None
    best_corners = None

    for i in range(len(hull)):
        p1 = hull[i]
        p2 = hull[(i + 1) % len(hull)]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        angle = math.atan2(dy, dx)
        ca = math.cos(-angle)
        sa = math.sin(-angle)
        rot = [(p[0] * ca - p[1] * sa, p[0] * sa + p[1] * ca) for p in hull]
        xs = [p[0] for p in rot]
        ys = [p[1] for p in rot]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)
        area = (maxx - minx) * (maxy - miny)
        if (best_area is None) or (area < best_area):
            best_area = area
            corners_rot = [(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)]
            ca = math.cos(angle)
            sa = math.sin(angle)
            best_corners = [(p[0] * ca - p[1] * sa, p[0] * sa + p[1] * ca) for p in corners_rot]

    return best_corners

def find_sides(corners):
    #returns two of the CCDs sides as vectors
    vectors = (corners - corners[0])[1:]
    for i in range(3):
        vec = vectors[i]
        others = np.delete(vectors, i, axis=0)
        cross1 = (others[0][0] * vec[1]) - (others[0][1] * vec[0])
        cross2 = (vec[0] * others[1][1]) - (vec[1] * others[1][0])
        if cross1*cross2 > 0:
            return (others[0], others[1])

def projection_norm(matrix, vector):
    """Calculated the projection's normal for each row vector in the matrix and returns a vector of the normals"""
    proj = (((matrix @ vector) / (vector @ vector))[:,np.newaxis]) * vector
    proj_norm = np.linalg.norm(proj, axis=1)
    return proj_norm

def find_streaks(event, bandwidth,n1,n2):
    points = make_point_list(event)
    corners = np.array(estimate_oriented_ccd_bounds(points))
    top_vector ,left_vector = find_sides(corners)
    if corners is not None:
        poly = corners + [corners[0]]
        px = [p[0] for p in poly]
        py = [p[1] for p in poly]

    #Choosing the axis parallel to chipx
    common_chipy = event['chipy'].mode()[0] 
    chip_row = event[event['chipy'] == common_chipy].sort_values('chipx')
    p1 = chip_row.iloc[0] 
    p2 = chip_row.iloc[-1]
    chipx_vec = np.array([p2['x'] - p1['x'], p2['y'] - p1['y']])
    chipx_vec = chipx_vec / np.linalg.norm(chipx_vec)
    top_norm = top_vector / np.linalg.norm(top_vector)
    left_norm = left_vector / np.linalg.norm(left_vector)
    top_score = np.abs(np.dot(top_norm, chipx_vec))
    left_score = np.abs(np.dot(left_norm, chipx_vec))
    if top_score > left_score:
        chipx_axis_vector = top_vector
        chipy_axis_vector = left_vector
    else:
        chipx_axis_vector = left_vector
        chipy_axis_vector = top_vector

    #Calculating the projections on two of the borders
    xy_matrix = np.hstack((np.array(event['x'])[:,np.newaxis] , np.array(event['y'])[:,np.newaxis]))
    xy_matrix = xy_matrix - corners[0] #Move origin to top left corner

    #Primary band division:
    event['proj_norm_v'] = projection_norm(xy_matrix,chipx_axis_vector)
    vaxis_size = np.linalg.norm(chipx_axis_vector)
    vaxis_middle = vaxis_size / 2
    v_stop = event['proj_norm_v'].values.max() + bandwidth
    vbands = np.arange(0, vaxis_size + bandwidth, bandwidth)
    event['vband_index'] = pd.cut(event.loc[:,'proj_norm_v'], bins=vbands, labels=False, include_lowest=True)
    max_vband = event['vband_index'].max()

    #Splitting the primary bands in the middle:
    event['proj_norm_h'] = projection_norm(xy_matrix,chipy_axis_vector)
    haxis_size = np.linalg.norm(chipy_axis_vector)
    haxis_middle = haxis_size / 2
    event['vband_side'] = np.where(event['proj_norm_h'] > haxis_middle, 1, 0)

    #Find absolute maxima in Primary bands:
    vband_source_count = event.groupby('vband_index')['proj_norm_v'].count()
    vband_source_count = vband_source_count.to_numpy()
    v_max_index = np.argmax(vband_source_count)
    #Distance from mean
    vband_source_count_nomax = np.delete(vband_source_count,v_max_index)
    vmean = np.mean(vband_source_count_nomax)
    # vmedian = np.median(vband_source_count)
    vsigma = np.std(vband_source_count_nomax)
    #check if photon count in maximal bin is larger than the median by n1*sigma
    condition1 = (vband_source_count[v_max_index] - vmean) > n1*vsigma
    is_streak = False
    if condition1:
        vband_source_count_0 = event[event['vband_side'] == 0].groupby('vband_index')['proj_norm_v'].count()
        vband_source_count_0 = vband_source_count_0.to_numpy()
        v_max_index_0 = np.argmax(vband_source_count_0)
        vband_source_count_1 = event[event['vband_side'] == 1].groupby('vband_index')['proj_norm_v'].count()
        vband_source_count_1 = vband_source_count_1.to_numpy()
        v_max_index_1 = np.argmax(vband_source_count_1)
        streak_band_0 = v_max_index_0
        streak_band_1 = v_max_index_1
        if v_max_index_0 == v_max_index_1 == v_max_index:
            vband_source_count_nomax_0 = np.delete(vband_source_count_0,v_max_index_0)
            vmean_0 = np.mean(vband_source_count_nomax_0)
            # vmedian_0 = np.median(vband_source_count_0)
            vsigma_0 = np.std(vband_source_count_nomax_0)
            vband_source_count_nomax_1 = np.delete(vband_source_count_1,v_max_index_1)
            vmean_1 = np.mean(vband_source_count_nomax_1)
            # vmedian_1 = np.median(vband_source_count_1)
            vsigma_1 = np.std(vband_source_count_nomax_1)
            condition2 = ((vband_source_count_0[v_max_index_0] - vmean_0) > n2*vsigma_0) & ((vband_source_count_1[v_max_index_1] - vmean_1) > n2*vsigma_1)
            if condition2:
                is_streak = True

        #--------------------------------------------------------------Plot-------------------------------------------------------------------
        #Histogram and photon counts per band plot:
        fig1,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,12), constrained_layout=True)
        ax1.plot(np.arange(0,(len(vband_source_count_0))), vband_source_count_0)
        ax1.set_title(f'photon count per band obsid: {obsid} side 1')
        ax1.set_xlabel('band index')
        ax1.set_ylabel('photon count per band')
        ax1.grid()
        ax1.set_xticks([streak_band_0], minor=True)
        ax1.set_xticklabels([str(streak_band_0)], minor=True)
        ax1.axvline(streak_band_0, color='r')
        ax1.tick_params(axis='x', which='minor', length=8, width=2, color='red', direction='out')
        ax2.plot(np.arange(0,(len(vband_source_count_1))), vband_source_count_1)
        ax2.set_title(f'photon count per band obsid: {obsid} side 2')
        ax2.set_xlabel('band index')
        ax2.set_ylabel('photon count per band')
        ax2.set_xticks([streak_band_1], minor=True)
        ax2.set_xticklabels([str(streak_band_1)], minor=True)
        ax2.axvline(streak_band_1, color='r')
        ax2.tick_params(axis='x', which='minor', length=8, width=2, color='red', direction='out')
        ax2.grid()
        ax3.hist(vband_source_count_0, bins=100)
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax3.grid()
        ax3.set_title(f'photon count histogram obsid: {obsid} side 1')
        ax3.set_xlabel('photon count')
        ax4.hist(vband_source_count_1, bins=100)
        ax4.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax4.grid()
        ax4.set_title(f'photon count histogram obsid: {obsid} side 2')
        ax4.set_xlabel('photon count')
        # plt.savefig(f".......png")

        #Map with borders:
        fig2,ax = plt.subplots(figsize=(6, 6))
        for i in range(max_vband):
            condition = (event['vband_index'] == i) & (event['vband_side'] == 0)
            points = make_point_list(event[condition])
            corners = corners = estimate_oriented_ccd_bounds(points)
            if corners is not None:
                poly = corners + [corners[0]]
                px = [p[0] for p in poly]
                py = [p[1] for p in poly]
            ax.plot(px,py,'m-',linewidth=0.5, label='oriented bounds')
        for i in range(max_vband):
            condition = (event['vband_index'] == i) & (event['vband_side'] == 1)
            points = make_point_list(event[condition])
            corners = corners = estimate_oriented_ccd_bounds(points)
            if corners is not None:
                poly = corners + [corners[0]]
                px = [p[0] for p in poly]
                py = [p[1] for p in poly]
            ax.plot(px,py,'m-',linewidth=0.5, label='oriented bounds')
        ax.scatter(event['x'],event['y'], c='k', s=1, alpha=0.07)
        condition = ((event['vband_index'] == streak_band_0) & (event['vband_side'] == 0)) | \
                    ((event['vband_index'] == streak_band_1) & (event['vband_side'] == 1))
        condition2 = ~condition
        # plt.scatter(event[condition]['x'], event[condition]['y'], c='r', s=1, alpha=0.07)
        # plt.scatter(event[condition2]['x'], event[condition2]['y'], c='k', s=1, alpha=0.07)
        #ax.plot(px,py,'m-',linewidth=2.0, label='oriented bounds')
        plt.title(f'map with band borders - obsid: {obsid}')
        plt.xlabel('x') 
        plt.ylabel('y')
        # ax.set_aspect('equal', adjustable='box')
        # plt.savefig(f"......png")
        plt.show()
    else:
        streak_band = v_max_index
        #------------------------------------------------------------Plot---------------------------------------------------------------------
        fig1,(ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        ax1.plot(np.arange(0,(len(vband_source_count))), vband_source_count)
        ax1.set_title(f'photon count per band obsid: {obsid}')
        ax1.set_xlabel('band index')
        ax1.set_ylabel('photon count per band')
        ax1.set_xticks([streak_band], minor=True)
        ax1.set_xticklabels([str(streak_band)], minor=True)
        ax1.axvline(streak_band, color='r')
        ax1.tick_params(axis='x', which='minor', length=8, width=2, color='red', direction='out')
        ax1.grid()
        ax2.hist(vband_source_count, bins=100)
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax2.set_title(f'photon count histogram obsid: {obsid}')
        ax2.set_xlabel('photon count')
        ax2.grid()
        # plt.savefig(f"C:/Users/bboyh/OneDrive/Desktop/{obsid}_histandplot.png")
        
        #Map with borders:
        fig2,ax = plt.subplots(figsize=(6, 6))
        for i in range(max_vband):
            condition = (event['vband_index'] == i)
            points = make_point_list(event[condition])
            corners = corners = estimate_oriented_ccd_bounds(points)
            if corners is not None:
                poly = corners + [corners[0]]
                px = [p[0] for p in poly]
                py = [p[1] for p in poly]
            ax.plot(px,py,'m-',linewidth=0.5, label='oriented bounds')
        # ax.scatter(event['x'],event['y'], c='k', s=1, alpha=0.07)
        condition = event['vband_index'] == streak_band
        condition2 = ~condition
        plt.scatter(event[condition]['x'], event[condition]['y'], c='r', s=1, alpha=0.07)
        plt.scatter(event[condition2]['x'], event[condition2]['y'], c='k', s=1, alpha=0.07)
        #ax.plot(px,py,'m-',linewidth=2.0, label='oriented bounds')
        plt.title(f'map with band borders - obsid: {obsid}')
        plt.xlabel('x') 
        plt.ylabel('y')
        # ax.set_aspect('equal', adjustable='box')
        # plt.savefig(f"......png")
        plt.show()


    return is_streak


obsid = 830
hdul = fits.open(f".....{obsid}.fits")
data = Table(hdul[1].data)
event_df = data.to_pandas()
event_df=event_df[event_df["ccd_id"] == 7] #When writing the function - maybe we will iterate over each CCD?
print(find_streaks(event_df, 20, 22, 8))