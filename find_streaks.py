import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
from astropy.io import fits
from pathlib import Path
import matplotlib.patches as mpatches

streak_flagged_list = None

def get_streak_flagged_list(filepath="cda_flagged_list.txt"):
  global streak_flagged_list
  if streak_flagged_list is None:
    df = pd.read_csv(filepath)
    df.columns = df.columns.str.strip().str.lower()
    obsids = df["obsid"].astype(int)
    ccds = df["chip_id"].astype(int)
    streak_flagged_list = set(zip(obsids, ccds))
  return streak_flagged_list

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
    if len(points) <= 2:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)
        return [(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)]

    def convex_hull(pts):
        pts = sorted(pts)

        def cross(o, a, b):
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

def find_streaks(event_path, plot_path, bandwidth=50, n1=10, n2=4, second_iter_factor=0.8, second=False):
    #Reading the fits
    with fits.open(event_path) as hdul:
        obsid = hdul['EVENTS'].header.get('OBS_ID')
        obsid = int(str(obsid).strip())
        data = hdul[1].data
        columns = {}
        for name in data.names:
            if data[name].ndim != 1:
                continue
            column = np.asarray(data[name])
            if column.dtype.byteorder not in ('=', '|'):
                column = column.byteswap().view(column.dtype.newbyteorder('='))
            columns[name] = column
        event_df = pd.DataFrame(columns)
    
    #setting up the path - creating a folder for the obsid:
    base_path = Path(plot_path)
    obsid_dir = Path(base_path) / f"obsid_{obsid}"
    obsid_dir.mkdir(parents=True, exist_ok=True)

    #Iterating over all CCDs:
    ccd_list = event_df['ccd_id'].dropna().unique().astype(int).tolist()
    chandra_labeled_obsid_ccd = get_streak_flagged_list()
    records = []
    for ccd in ccd_list:
        ccd_dir = obsid_dir / f"ccd_{ccd}"
        ccd_dir.mkdir(parents=True, exist_ok=True)
        obsid_in_labeled = (obsid, ccd) in chandra_labeled_obsid_ccd
        event_ccd_specific = event_df[event_df['ccd_id'] == ccd]
        is_streak, sigma_distance, sigma_distance_0, sigma_distance_1, mean, mean_0, mean_1 = find_streaks_logic(event_ccd_specific, obsid, ccd, bandwidth, n1, n2, ccd_dir, second_iter_factor, second)
        records.append({
        "ccd": ccd,
        "is_streak": bool(is_streak),
        "sigma_distance": sigma_distance,
        "sigma_distance_side1": sigma_distance_0,
        "sigma_distance_side2": sigma_distance_1,
        "mean_photon_count": mean,
        "mean_photon_count_side1": mean_0,
        "mean_photon_count_side2": mean_1,
        "labeled_as_streaked_in_cda": obsid_in_labeled,
        })

    verdict_df = pd.DataFrame(records)
    df_path = Path(obsid_dir) / f"verdict_obsid_{obsid}.csv"
    verdict_df.to_csv(df_path, index=False)
    return verdict_df

def find_streaks_logic(event_df, obsid, ccd, bandwidth, n1, n2, plot_path, second_iter_factor=0.8, second=False):
    #Pre-initialization of second check's parameters:
    event_df = event_df.copy()
    sigma_distance_0 = None
    sigma_distance_1 = None
    mean_0 = None
    mean_1 = None

    points = make_point_list(event_df)
    corners = np.array(estimate_oriented_ccd_bounds(points))
    top_vector ,left_vector = find_sides(corners)

    #Choosing the axis parallel to chipx
    common_chipy = event_df['chipy'].mode()[0] 
    chip_row = event_df[event_df['chipy'] == common_chipy].sort_values('chipx')
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
    xy_matrix = np.hstack((np.array(event_df['x'])[:,np.newaxis] , np.array(event_df['y'])[:,np.newaxis]))
    xy_matrix = xy_matrix - corners[0] #Move origin to top left corner

    #Primary band division:
    event_df['proj_norm'] = projection_norm(xy_matrix,chipx_axis_vector)
    axis_size = np.linalg.norm(chipx_axis_vector)
    bands = np.arange(0, axis_size + bandwidth, bandwidth)
    event_df['band_index'] = pd.cut(event_df.loc[:,'proj_norm'], bins=bands, labels=False, include_lowest=True)

    #Splitting the primary bands in the middle:
    event_df['proj_norm_orth'] = projection_norm(xy_matrix,chipy_axis_vector)
    orth_axis_size = np.linalg.norm(chipy_axis_vector)
    orth_axis_middle = orth_axis_size / 2
    event_df['band_side'] = np.where(event_df['proj_norm_orth'] > orth_axis_middle, 1, 0)

    #Find absolute maxima in Primary bands:
    band_source_count, _ = np.histogram(event_df['proj_norm'], bins=bands)
    max_index = np.argmax(band_source_count)

    #Distance from mean
    band_source_count_nomax = np.delete(band_source_count,max_index)
    mean = np.mean(band_source_count_nomax)
    sigma = np.std(band_source_count_nomax)
    streak_band = max_index
    #check if photon count in maximal bin is larger than the median by n1*sigma
    sigma_distance = (band_source_count[max_index] - mean) / sigma if sigma > 0 else 0.0
    condition1 = n1 < sigma_distance
    is_streak = False
    if condition1:
        band_source_count_0, _ = np.histogram(event_df[event_df['band_side'] == 0]['proj_norm'], bins=bands)
        max_index_0 = np.argmax(band_source_count_0)
        band_source_count_1, _ = np.histogram(event_df[event_df['band_side'] == 1]['proj_norm'], bins=bands)
        max_index_1 = np.argmax(band_source_count_1)
        streak_band_0 = max_index_0
        streak_band_1 = max_index_1
        band_source_count_nomax_0 = np.delete(band_source_count_0,max_index_0)
        mean_0 = np.mean(band_source_count_nomax_0)
        sigma_0 = np.std(band_source_count_nomax_0)
        band_source_count_nomax_1 = np.delete(band_source_count_1,max_index_1)
        mean_1 = np.mean(band_source_count_nomax_1)
        sigma_1 = np.std(band_source_count_nomax_1)
        sigma_distance_0 = (band_source_count_0[max_index_0] - mean_0) / sigma_0 if sigma_0 > 0 else 0.0
        sigma_distance_1 = (band_source_count_1[max_index_1] - mean_1) / sigma_1 if sigma_1 > 0 else 0.0
        if max_index_0 == max_index_1 == max_index:
            condition2 = (sigma_distance_0 > n2) and (sigma_distance_1 > n2)
            if condition2:
                is_streak = True

    #Second iteration in case of False:
    if not is_streak and not second:
        return find_streaks_logic(event_df, obsid, ccd, bandwidth * second_iter_factor, n1, n2, plot_path, second=True)
    
    #Plots for first check:
    fig, (ax_hist, ax_plot, ax_map) = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(f'main band check obsid {obsid} ccd {ccd}', fontsize=16)
    ax_plot.plot(np.arange(0,(len(band_source_count))), band_source_count)
    ax_plot.set_title(f'photon count per band')
    ax_plot.set_xlabel('band index')
    ax_plot.set_ylabel('photon count per band')
    ax_plot.grid()
    ax_plot.set_xticks([streak_band], minor=True)
    ax_plot.set_xticklabels([str(streak_band)], minor=True)
    ax_plot.axvline(streak_band, color='r')
    ax_plot.tick_params(axis='x', which='minor', length=8, width=2, color='red', direction='out')

    counts, bins, _ = ax_hist.hist(band_source_count, bins=100)
    ax_hist.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax_hist.grid()
    ax_hist.set_title(f'photon count histogram')
    ax_hist.set_xlabel('photon count')
    y_pos = counts.max() * 0.6
    x_peak = mean
    x_outlier = band_source_count[max_index]
    ax_hist.annotate('', xy=(x_outlier, y_pos), xytext=(x_peak, y_pos), 
                    arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
    ax_hist.text((x_peak + x_outlier)/2, y_pos + (counts.max() * 0.02), 
                f'{sigma_distance:.2f}', color='red', 
                ha='center', va='bottom', fontweight='bold')
    ax_hist.axvline(mean, color='gray', linestyle='--', alpha=0.7)
    ax_hist.set_xticks([mean], minor=True)

    v_hat = chipx_axis_vector / axis_size
    h_hat = chipy_axis_vector / orth_axis_size
    origin = corners[0]
    for v_val in bands:
        p_start = origin + (v_val * v_hat)
        p_end = origin + (v_val * v_hat) + (orth_axis_size * h_hat)
        ax_map.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], 'm-', linewidth=0.5)
    for h_val in [0, orth_axis_size]:
        p_start = origin + (h_val * h_hat)
        p_end = origin + (axis_size * v_hat) + (h_val * h_hat)
        ax_map.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], 'm-', linewidth=0.5)
    v_start = bands[streak_band]
    v_end = bands[streak_band + 1]
    corner_bl = origin + (v_start * v_hat)
    corner_br = origin + (v_start * v_hat) + (orth_axis_size * h_hat)
    corner_tr = origin + (v_end * v_hat) + (orth_axis_size * h_hat)
    corner_tl = origin + (v_end * v_hat)
    band_polygon = mpatches.Polygon(
        [corner_bl, corner_br, corner_tr, corner_tl],
        closed=True,
        facecolor='red',
        alpha=0.35,       
        edgecolor=None,
        zorder=1         
    )
    ax_map.add_patch(band_polygon)
    ax_map.scatter(event_df['x'], event_df['y'], c='k', s=1, alpha=0.05, zorder=2)
    ax_map.set_title(f'map with band borders')
    ax_map.set_xlabel('x') 
    ax_map.set_ylabel('y')
    plot1_path = Path(plot_path) / f"{obsid}_first_check_ccd_{ccd}.png"
    plt.savefig(plot1_path)
    plt.close(fig)

    #Plots for second check:
    if condition1:
        layout = """
            AABB
            CCDD
            .MM.
            """
        
        fig1, ax_dict = plt.subplot_mosaic(layout, figsize=(12, 18), constrained_layout=True)
        fig1.suptitle(f'split band check obsid {obsid} ccd {ccd}', fontsize=16)
        ax1 = ax_dict['A']
        ax2 = ax_dict['B']
        ax3 = ax_dict['C']
        ax4 = ax_dict['D']
        ax_map = ax_dict['M']
        ylim_plot = max(max(band_source_count_0),max(band_source_count_1))
    
        ax1.plot(np.arange(0,(len(band_source_count_0))), band_source_count_0)
        ax1.set_title(f'photon count per band side 1')
        ax1.set_xlabel('band index')
        ax1.set_ylabel('photon count per band')
        ax1.grid()
        ax1.set_xticks([streak_band_0], minor=True)
        ax1.set_xticklabels([str(streak_band_0)], minor=True)
        ax1.axvline(streak_band_0, color='r')
        ax1.set_ylim(0,ylim_plot)
        ax1.tick_params(axis='x', which='minor', length=8, width=2, color='red', direction='out')
    
        ax2.plot(np.arange(0,(len(band_source_count_1))), band_source_count_1)
        ax2.set_title(f'photon count per band side 2')
        ax2.set_xlabel('band index')
        ax2.set_ylabel('photon count per band')
        ax2.set_xticks([streak_band_1], minor=True)
        ax2.set_xticklabels([str(streak_band_1)], minor=True)
        ax2.axvline(streak_band_1, color='r')
        ax2.set_ylim(0,ylim_plot)
        ax2.tick_params(axis='x', which='minor', length=8, width=2, color='red', direction='out')
        ax2.grid()
    
        counts_0, bins_0, _ = ax3.hist(band_source_count_0, bins=100)
        counts_1, bins_1, _ = ax4.hist(band_source_count_1, bins=100)
        ylim_hist = max(max(counts_0),max(counts_1))
    
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax3.grid()
        ax3.set_title(f'photon count histogram side 1')
        ax3.set_xlabel('photon count')
        y_pos_0 = counts_0.max() * 0.6
        x_peak_0 = mean_0
        x_outlier_0 = band_source_count_0[max_index_0]
        ax3.annotate('', xy=(x_outlier_0, y_pos_0), xytext=(x_peak_0, y_pos_0), 
                        arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
        ax3.text((x_peak_0 + x_outlier_0)/2, y_pos_0 + (counts_0.max() * 0.02), 
                    f'{sigma_distance_0:.2f}', color='red', 
                    ha='center', va='bottom', fontweight='bold')
        ax3.axvline(mean_0, color='gray', linestyle='--', alpha=0.7)
        ax3.set_xticks([mean_0], minor=True)
        ax3.tick_params(axis='x', which='minor', length=8, width=2, color='gray', direction='out')
        
        ax3.set_ylim(0,ylim_hist)
    
        ax4.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax4.grid()
        ax4.set_title(f'photon count histogram side 2')
        ax4.set_xlabel('photon count')
        y_pos_1 = counts_1.max() * 0.6
        x_peak_1 = mean_1
        x_outlier_1 = band_source_count_1[max_index_1]
        ax4.annotate('', xy=(x_outlier_1, y_pos_1), xytext=(x_peak_1, y_pos_1), 
                        arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
        ax4.text((x_peak_1 + x_outlier_1)/2, y_pos_1 + (counts_1.max() * 0.02), 
                    f'{sigma_distance_1:.2f}', color='red', 
                    ha='center', va='bottom', fontweight='bold')
        ax4.axvline(mean_1, color='gray', linestyle='--', alpha=0.7)
        ax4.set_xticks([mean_1], minor=True)
        ax4.set_ylim(0,ylim_hist)
    
        v_hat = chipx_axis_vector / axis_size
        h_hat = chipy_axis_vector / orth_axis_size
        origin = corners[0]
        s0_p1 = origin
        s0_p2 = origin + (axis_size * v_hat)
        s0_p3 = origin + (axis_size * v_hat) + (orth_axis_middle * h_hat)
        s0_p4 = origin + (orth_axis_middle * h_hat)
        ax_map.fill([s0_p1[0], s0_p2[0], s0_p3[0], s0_p4[0]], 
                    [s0_p1[1], s0_p2[1], s0_p3[1], s0_p4[1]], 
                    color='cornflowerblue', alpha=0.15, zorder=1)
        s1_p1 = origin + (orth_axis_middle * h_hat)
        s1_p2 = origin + (axis_size * v_hat) + (orth_axis_middle * h_hat)
        s1_p3 = origin + (axis_size * v_hat) + (orth_axis_size * h_hat)
        s1_p4 = origin + (orth_axis_size * h_hat)
        ax_map.fill([s1_p1[0], s1_p2[0], s1_p3[0], s1_p4[0]], 
                    [s1_p1[1], s1_p2[1], s1_p3[1], s1_p4[1]], 
                    color='lightcoral', alpha=0.15, zorder=1)
        for v_val in bands:
            p_start = origin + (v_val * v_hat)
            p_end = origin + (v_val * v_hat) + (orth_axis_size * h_hat)
            ax_map.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], 'm-', linewidth=0.5)
        for h_val in [0, orth_axis_middle, orth_axis_size]:
            p_start = origin + (h_val * h_hat)
            p_end = origin + (axis_size * v_hat) + (h_val * h_hat)
            ax_map.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], 'm-', linewidth=0.5)
        v0_start = bands[streak_band_0]
        v0_end = bands[streak_band_0 + 1]
        s0_b_p1 = origin + (v0_start * v_hat)
        s0_b_p2 = origin + (v0_end * v_hat)
        s0_b_p3 = origin + (v0_end * v_hat) + (orth_axis_middle * h_hat)
        s0_b_p4 = origin + (v0_start * v_hat) + (orth_axis_middle * h_hat)
        ax_map.fill(
            [s0_b_p1[0], s0_b_p2[0], s0_b_p3[0], s0_b_p4[0]],
            [s0_b_p1[1], s0_b_p2[1], s0_b_p3[1], s0_b_p4[1]],
            color='red', alpha=0.35, zorder=2
        )
        v1_start = bands[streak_band_1]
        v1_end = bands[streak_band_1 + 1]
        s1_b_p1 = origin + (v1_start * v_hat) + (orth_axis_middle * h_hat)
        s1_b_p2 = origin + (v1_end * v_hat) + (orth_axis_middle * h_hat)
        s1_b_p3 = origin + (v1_end * v_hat) + (orth_axis_size * h_hat)
        s1_b_p4 = origin + (v1_start * v_hat) + (orth_axis_size * h_hat)
        ax_map.fill(
            [s1_b_p1[0], s1_b_p2[0], s1_b_p3[0], s1_b_p4[0]],
            [s1_b_p1[1], s1_b_p2[1], s1_b_p3[1], s1_b_p4[1]],
            color='red', alpha=0.35, zorder=2
        )
        ax_map.scatter(event_df['x'], event_df['y'], c='k', s=1, alpha=0.07, zorder=3)
        side0_patch = mpatches.Patch(color='cornflowerblue', alpha=0.4, label='Side 1')
        side1_patch = mpatches.Patch(color='lightcoral', alpha=0.4, label='Side 2')
        streak_patch = mpatches.Patch(color='red', alpha=0.35, label='Detected Band')
        ax_map.legend(handles=[side0_patch, side1_patch, streak_patch], loc='best', framealpha=0.9)
        ax_map.set_title(f'map with split band borders')
        ax_map.set_xlabel('x') 
        ax_map.set_ylabel('y')
        plot2_path = Path(plot_path) / f"{obsid}_second_check_{ccd}.png"
        plt.savefig(plot2_path)
        plt.close(fig1)

    return is_streak, sigma_distance, sigma_distance_0, sigma_distance_1, mean, mean_0, mean_1