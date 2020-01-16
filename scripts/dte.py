import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm, tqdm_pandas
import itertools
from multiprocessing import Pool

def calculate_dte(inoculum_size, proportion_of_taxa, viability=1, num_wells=96):
    '''
    Main method to calculate DTE probabilities
    :param inoculum_size: The number of cells added to each well
    :param proportion_of_taxa: The relative abundance of the target taxa in the inoculum
    :param viability: The viability of the target taxa
    :param num_wells: The number of wells being inoculated
    :return: the number of positive wells,
            the number of pure wells,
            the number of positive wells containing the target taxon,
            the number of pure wells containing the target taxon
    '''
    inocula=np.random.poisson(inoculum_size, num_wells)
    taxon_wells=np.random.binomial(inocula, proportion_of_taxa)
    positive_wells = np.sum(inocula >=1)
    single_wells = np.sum(inocula == 1)
    #This captures all the wells with greater than 0 cells where the number
    #of cells in it are all from the selected taxon
    positive_wells_with_taxon=np.sum(taxon_wells>=1)
    taxon_pure_wells = np.sum((taxon_wells == inocula) & (inocula >0))
    if viability <1:
        pure_wells_to_test =  taxon_wells[(taxon_wells == inocula) & (inocula >0)]
        viable_pure_wells=np.random.binomial(pure_wells_to_test, viability)
        taxon_pure_wells = np.sum(viable_pure_wells>=1)
        unviable_taxon_pure = np.sum(viable_pure_wells==0)
        positive_wells = positive_wells - unviable_taxon_pure
    return positive_wells, single_wells, positive_wells_with_taxon, taxon_pure_wells


def get_CI(array_of_numbers, as_string=True, ci=95):
  median_value = np.median(array_of_numbers)
  ci_frac = (100 - ci) / 2
  percentiles = np.percentile(array_of_numbers, [ci_frac, 100 - ci_frac])
  if as_string:
    return f'{median_value:.1f} wells ({percentiles[0]:.1f}-{percentiles[1]:.1f}, 95% CI)'
  else:
    return median_value, percentiles[0], percentiles[1]


def bootstrap_dte(cells_per_well, rel_abund, viability, num_wells):
    positive_wells = np.zeros(number_of_experiments)
    single_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)

    for i in range(number_of_experiments):
        positive_wells[i], single_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(cells_per_well, rel_abund, viability=viability, num_wells=num_wells)

    taxon_pure_med, taxon_pure_low, taxon_pure_high = get_CI(taxon_pure_wells, as_string=False)
    taxon_positive_med, taxon_positive_low, taxon_positive_high = get_CI(taxon_positive_wells, as_string=False)
    positive_med, positive_low, positive_high = get_CI(positive_wells, as_string=False)
    single_med, single_low, single_high = get_CI(single_wells, as_string=False)

    return pd.Series([cells_per_well, rel_abund, viability,
                      num_wells, number_of_experiments,
                      positive_med, positive_low, positive_high,
                      single_med, single_low, single_high,
                      taxon_positive_med, taxon_positive_low, taxon_positive_high,
                      taxon_pure_med, taxon_pure_low, taxon_pure_high],
                     index=['cells_per_well', 'rel_abund', 'viability',
                            'num_wells', 'number_of_experiments',
                            'positive_well_med','positive_well_95pc_low', 'positive_well_95pc_high',
                            'single_med', 'single_low', 'single_high',
                            'taxon_positive_med', 'taxon_positive_95pc_low', 'taxon_positive_95pc_high',
                            'pure_well_med', 'pure_well_95pc_low', 'pure_well_95pc_high'])


def process_row(row):
    positive_wells = np.zeros(number_of_experiments)
    single_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)


    for i in range(number_of_experiments):
        positive_wells[i], single_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(row.cells_per_well, row.rel_abund,num_wells=row.num_wells_inoculated)

    taxon_pure_med, taxon_pure_low, taxon_pure_high = get_CI(taxon_pure_wells, as_string=False)

    if row.num_isolates > taxon_pure_high:
        p_value = np.sum(taxon_pure_wells >= row.num_isolates) / number_of_experiments
        within_estimates = 'no'
        deviance = row.num_isolates - taxon_pure_high
    elif row.num_isolates < taxon_pure_low:
        p_value = np.sum(taxon_pure_wells <= row.num_isolates) / number_of_experiments
        within_estimates = 'no'
        deviance = row.num_isolates - taxon_pure_low
    else:
        p_value = ''
        within_estimates = 'yes'
        deviance = ''

    taxon_positive_med, taxon_positive_low, taxon_positive_high = get_CI(taxon_positive_wells, as_string=False)
    positive_med, positive_low, positive_high = get_CI(positive_wells, as_string=False)
    single_med, single_low, single_high = get_CI(single_wells, as_string=False)
    return pd.Series([row.ID, row.Site, row.num_wells_inoculated,
                      positive_med, positive_low, positive_high,
                      taxon_positive_med, taxon_positive_low, taxon_positive_high,
                      single_med, single_low, single_high,
                      taxon_pure_med, taxon_pure_low, taxon_pure_high, within_estimates, deviance, p_value],
                     index=['ID', 'Site', 'num_inoculated_wells',
                            'positive_well_med', 'positive_well_95pc_low','positive_well_95pc_high',
                            'taxon_positive_med', 'taxon_positive_95pc_low', 'taxon_positive_95pc_high',
                                'single_well_med', 'single_well_95pc_low', 'single_well_95pc_high',
                               'taxon_pure_med',
                               'taxon_pure_95pc_low',
                               'taxon_pure_95pc_high', 'within_estimates', 'deviance', 'p_value'])

def run_bootstrap_per_ASV(input_file):
    print(f'Reading in ASV file - {input_file} and running with {number_of_experiments} simulations')
    df = pd.read_csv(input_file)
    tqdm_pandas(tqdm())
    results = df.progress_apply(process_row, axis=1)
    joined_results = df.merge(results, how='left', on=['ID', 'Site'])
    joined_results.to_csv(f'./data/ASV_cultivation_bootstrapped_numbers_max_viability_{number_of_experiments}_simulations.joined.csv', index=False)

def process(job):
    test_proportions = np.arange(start=0.0, stop=1, step=0.001)
    for p in test_proportions:
        record = bootstrap_dte(cells_per_well=job[1],
                                                                        rel_abund=p,
                                                                        viability=1,
                                                                        num_wells=job[0])
        if (record.taxon_pure_low >= 1):
            return((job[0], job[1], p, record.taxon_pure_low))

def estimate_required_proportion_for_underepresentation():
    print('*** Testing required proportions for under-representation ***')

    number_of_wells = np.arange(start=92, stop=92*101, step=92)
    inoculum=[1,1.27, 1.96, 2,3,4,5]
    jobs = list(itertools.product(number_of_wells, inoculum))
    with Pool(16) as p:
            results = list(tqdm(p.imap(process, jobs), total=len(jobs), desc='Estimating proportion limits:'))
    results_df = pd.DataFrame(results, columns=['number_of_inoculated_wells', 'inoculum', 'p_i', 'taxon_pure_low'])
    results_df.to_csv(f'./data/estimate_min_num_wells_to_identify_viability_sensitivity_{number_of_experiments}_simulations.csv', index=False)
    print(results_df)

def calculate_interaction_rel_abundance_inoculum():
    results = []
    relative_abundance = np.arange(start=0, stop=1.01, step=0.01)
    inoculum = np.arange(start=1, stop=11, step=1)

    for cells_per_well, rel_abund in itertools.product(inoculum, relative_abundance):
        print(f'Running {cells_per_well} cells per well with a relative abundance of {rel_abund}')
        record  = bootstrap_dte(cells_per_well, rel_abund, num_wells= 9999)
        results.append((record.cells_per_well, record.rel_abund, record.taxon_pure_med, record.taxon_pure_low, record.taxon_pure_high))

    results_df = pd.DataFrame(results, columns=['cells_per_well', 'relative_abundance', 'taxon_pure_med', 'taxon_pure_low', 'taxon_pure_high'])
    results_df.to_csv(f'./data/inoculum_rel_abundance_interaction_plot_data_{number_of_experiments}_simulations.csv', index=False)

def process_viability_multiprocessor(items):
    s = bootstrap_dte(items[0], items[1], items[2], items[3])
    taxon_pure_med = s.pure_well_med
    taxon_pure_low = s.pure_well_95pc_low
    taxon_pure_high = s.pure_well_95pc_high
    if items[4] >=taxon_pure_low and items[4] <= taxon_pure_high:
        return items[2]
    else:
        return 0

def estimate_viability(cells_per_well, rel_abund, observed_pure, num_wells, test_min=0, test_max=1, initial_test_step=0.01):
    viability_range = np.arange(start=test_min, stop=test_max, step=initial_test_step)
    item_list = [(cells_per_well, rel_abund, v, num_wells, observed_pure) for v in viability_range]
    print(f'Testing viability values between {test_min} and {test_max} in increments of {initial_test_step} to see at what viability an observed value of {observed_pure} wells out of {num_wells} wells falls within 95%CI intervals, given an inoculum of {cells_per_well:.2f} cells per well')
    with Pool(16) as p:
        viability_values = list(p.imap(process_viability_multiprocessor, item_list))
        #print(viability_values)

    cleaned_viability_values = [x for x in viability_values if x !=0]
    try:
        min = np.min(cleaned_viability_values)
        max= np.max(cleaned_viability_values)

        new_min = np.max([0,min-initial_test_step])
        new_max = np.min([1,max+initial_test_step])
        new_step = initial_test_step/10
        print(f'''Minimum and maximum values are {min} and {max}, respectively - refining:\nTesting values between {new_min:.3f} and {new_max:.3f} in increments of {new_step}''')

        viability_range = np.arange(start=new_min, stop=new_max, step=new_step)
        item_list = [(cells_per_well, rel_abund, v, num_wells, observed_pure) for v in viability_range]

        with Pool(16) as p:
            viability_values = list(p.imap(process_viability_multiprocessor, item_list))
            #print(viability_values)
        cleaned_viability_values = [x for x in viability_values if x !=0]

        min = np.min(cleaned_viability_values)
        max= np.max(cleaned_viability_values)
    except ValueError:
        #To get to here, the list of viability values is full of zeros, which means
        #the range hasn't been found. This is most likely due to extremely low viability
        #So, we set new_min to 0 and new_max to 0.1 and try at steps of 0.0001
        new_step = initial_test_step
        while len(cleaned_viability_values) ==0 and new_step > 0.00001:
            new_min = 0
            new_max = new_step
            new_step=new_step/10
            print(f'Trying values between {new_min} and {new_max} with increments of {new_step}')
            viability_range = np.arange(start=new_min, stop=new_max, step=new_step)
            item_list = [(cells_per_well, rel_abund, v, num_wells, observed_pure) for v in viability_range]

            with Pool(16) as p:
                viability_values = list(p.imap(process_viability_multiprocessor, item_list))
            #print(viability_values)
            cleaned_viability_values = [x for x in viability_values if x !=0]
            #print(cleaned_viability_values)
        if len(cleaned_viability_values) > 0:
            min = np.min(cleaned_viability_values)
            max= np.max(cleaned_viability_values)
        else:
            print(f'Gave up - estimated viability is really rather low.')
            min=np.NaN
            max=np.NaN


    print(f'The final bootstrapped values for viability are {min:1.3%}-{max:1.3%} %')
    return min, max


def process_viability(row):
    positive_wells = np.zeros(number_of_experiments)
    single_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)

    viability_range = np.arange(start=0, stop=1, step=0.01)[::-1]
    max_viability = 0

    for v in viability_range:
        print(f'*** Processing ID:{row.ID} at site: {row.Site}, testing viability {v:.2f} ***')

        for i in range(number_of_experiments):
            positive_wells[i], single_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(row.cells_per_well, row.rel_abund, viability=v, num_wells=row.num_wells_inoculated)

        taxon_pure_med, taxon_pure_low, taxon_pure_high = get_CI(taxon_pure_wells, as_string=False)
#        if row.num_isolates >=taxon_pure_low and row.num_isolates <= taxon_pure_high:
#            max_viability = v
#            break
#    return pd.Series([row.ID, row.Site, row.rel_abund, row.num_wells_inoculated,
#                    row.num_isolates, row.cells_per_well,
#                    row.taxon_pure_95pc_low, min_viability, max_viability], index=['ID', 'Site',
#                                                                'rel_abund', 'num_inoculated_wells',
#                                                                'num_isolates', 'cells_per_well', 'taxon_pure_95pc_low', 'min_viability','max_viability'])

def process_viability(row):
    min, max = estimate_viability(row.cells_per_well, row.rel_abund, row.num_isolates, row.num_wells_inoculated)
    return pd.Series([row.ID, row.Site, row.rel_abund, row.num_wells_inoculated,  row.num_isolates, row.cells_per_well, min, max], index=['ID','Site','rel_abund','num_inoculated_wells','num_isolates','cells_per_well','min_viability', 'max_viability'])



def estimate_viability_range_for_underrepresented_taxa():
    print('*** Calculating maximum viability to explain observed number of wells ***')

    df = pd.read_csv(f'./data/ASV_cultivation_bootstrapped_numbers_max_viability_{number_of_experiments}_simulations.joined.csv')
    lower_than_expected_df = df[df.deviance <0]
    tqdm_pandas(tqdm())
    results = lower_than_expected_df.progress_apply(process_viability, axis=1)
    results.to_csv(f'./data/estimate_viability_range_for_underrepresented_taxa_{number_of_experiments}_simulations.csv', index=False)

def make_inoculum_rel_abundance_interaction_plot():
    fig = plt.figure(figsize=[6,6], dpi=100)
    ax = Axes3D(fig)
    ax.set_xlabel('Inoculum', fontsize=12)
    ax.set_ylabel('Relative abundance', fontsize=12)
    ax.set_zlabel('% Pure', fontsize=12)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    df = pd.read_csv(f'./data/inoculum_rel_abundance_interaction_plot_data_{number_of_experiments}_simulations.csv')

    ax = fig.gca(projection='3d')
    ax.plot_trisurf(df['cells_per_well'], df['relative_abundance'], df['taxon_pure_med']/9999*100,
                    cmap=plt.cm.jet, linewidth=0.0, antialiased=True)
    ax.view_init(azim=-10)
    plt.show()


def process_3d_plot_item(item):
    """
    Utility function for multiprocessing
    :param item: The item to be processed
    :return:
    """
    return bootstrap_dte(item[0], item[1], item[2], num_wells=number_of_wells_to_simulate)

def generate_3d_plot_data():
    """
    Generates the dataframe for simulated relative abundances vs viability vs inoculum size
    :return:
    """
    inocula = np.arange(start=1, stop=11, step=1)
    viability = np.arange(start=0, stop=1.05, step=0.05)
    rel_abund = np.arange(start=0, stop=1.01, step=0.01)
    items = list(itertools.product(inocula, rel_abund, viability))
    with Pool(16) as p:
        results = list(tqdm(p.imap(process_3d_plot_item, items), total=len(items),
                            desc='Generating data for 3d plots:'))
    df = pd.DataFrame(results)
    df.to_csv(f'./data/3d_plot_data_{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.csv', index=False)





def generate_3d_plot(x_label, y_label, z_label, x, y, z, filename):
    """
    Generates the 3d plots
    :param y_label: The label of the y axis
    :param x: The x axis values to plot
    :param y: The y axis values to plot
    :param z: The z axis values to plot
    :param filename: The filename to save it to
    :return: None
    """
    print(f'Generating plot for {x_label} vs {y_label} effect on {z_label} (may take some time).....')
    fig = plt.figure(figsize=[6, 6], dpi=100)
    ax = Axes3D(fig)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_zlabel(f' % {z_label}', fontsize=12)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    ax = fig.gca(projection='3d')

    def init():
        # Plot the surface.
        ax.plot_trisurf(x, y, z,
                        cmap=plt.cm.jet, linewidth=0.0, antialiased=True)
        return fig,

    def animate(i):
        # azimuth angle : 0 deg to 360 deg
        ax.view_init(elev=10, azim=i +1 )
        return fig,

    # Animate
    ani = animation.FuncAnimation(fig, animate, init_func=init,
                                  frames=360, interval=20, blit=True)

    ani.save(f'./plots/{filename}', writer='ffmpeg', fps=30)


def main():
    global number_of_experiments
    number_of_experiments = 9999

    global number_of_wells_to_simulate
    number_of_wells_to_simulate=460
    #run_bootstrap_per_ASV('./data/ASV_cultivation_numbers.csv')
    #estimate_required_proportion_for_underepresentation()
    #calculate_interaction_rel_abundance_inoculum()
    #estimate_viability_range_for_underrepresented_taxa()
    #make_inoculum_rel_abundance_interaction_plot()

    #generate_3d_plot_data()

    #df = pd.read_csv(f'./data/3d_plot_data_{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.csv')
    #full_viability = df[df['viability'] == 1]

    #generate_3d_plot(x_label='Inoculum',
    #                 y_label='Relative Abundance',
    #                 z_label='Taxon Positive Well',
    #                 x=full_viability['cells_per_well'],
    #                 y=full_viability['rel_abund'],
    #                 z=full_viability['taxon_positive_med']/number_of_wells_to_simulate*100,
    #                 filename=f'taxon_positive_med.viability_1.{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.3d_plot.mp4')

    #generate_3d_plot(x_label='Inoculum',
    #                 y_label='Relative Abundance',
    #                 z_label='Pure Well',
    #                 x=full_viability['cells_per_well'],
    #                 y=full_viability['rel_abund'],
    #                 z=full_viability['pure_well_med']/number_of_wells_to_simulate*100,
    #                 filename=f'pure_well_med.viability_1.{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.3d_plot.mp4')

    #generate_3d_plot(x_label='Inoculum',
    #                 y_label='Relative Abundance',
    #                 z_label='Taxon positive well / Positive Well',
    #                 x=full_viability['cells_per_well'],
    #                 y=full_viability['rel_abund'],
    #                 z=full_viability['taxon_positive_med'] / full_viability['positive_well_med'] * 100,
    #                 filename=f'taxon_positive_med_frac.viability_1.{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.3d_plot.mp4')

    #full_abundance = df[df['rel_abund'] == 1]

    #generate_3d_plot(x_label='Inoculum',
    #                 y_label='Viability',
    #                 z_label='Positive Well',
    #                 x=full_abundance['cells_per_well'],
    #                 y=full_abundance['viability'],
    #                 z=full_abundance['positive_well_med'] / number_of_wells_to_simulate * 100,
    #                 filename=f'positive_well_med.rel_abund_1.{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.3d_plot.mp4')

    #inoculum_2 = df[df['cells_per_well'] == 2]

    #generate_3d_plot(x_label='Relative Abundance',
    #                 y_label='Viability',
    #                 z_label='Pure Well',
    #                 x=inoculum_2['rel_abund'],
    #                 y=inoculum_2['viability'],
    #                 z=inoculum_2['pure_well_med'] / number_of_wells_to_simulate * 100,
    #                 filename=f'pure_well_med.inoculum_2.{number_of_wells_to_simulate}_wells.{number_of_experiments}_bootstraps.3d_plot.mp4')
    #cells_per_well, rel_abund, observed_pure, num_wells
    #ID,Site,rel_abund,num_inoculated_wells,num_isolates,cells_per_well,taxon_pure_95pc_low,max_viability
    #7471,FWC3,0.067158385,460,2,2.0,4.0,0.79
    estimate_viability_range_for_underrepresented_taxa()

    estimate_viability_range_for_




if __name__=="__main__":
    main()
