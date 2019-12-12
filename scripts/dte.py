import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm, tqdm_pandas
import itertools
from multiprocessing import Pool

def calculate_dte(inoculum_size, proportion_of_taxa, viability=1, num_wells=96):
    inocula=np.random.poisson(inoculum_size, num_wells)
    taxon_wells=np.random.binomial(inocula, proportion_of_taxa)
    positive_wells = np.sum(inocula >=1)
    pure_wells = np.sum(inocula==1)
    #This captures all the wells with greater than 0 cells where the number
    #of cells in it are all from the selected taxon
    positive_wells_with_taxon=np.sum(taxon_wells>=1)
    taxon_pure_wells = np.sum((taxon_wells == inocula) & (inocula >0))
    if viability <1:
        pure_wells_to_test =  taxon_wells[(taxon_wells == inocula) & (inocula >0)]
        viable_pure_wells=np.random.binomial(pure_wells_to_test, viability)
        taxon_pure_wells = np.sum(viable_pure_wells>=1)
    return positive_wells, pure_wells, positive_wells_with_taxon, taxon_pure_wells


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
    pure_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)
    for i in range(number_of_experiments):
        positive_wells[i], pure_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(cells_per_well, rel_abund, viability=viability, num_wells=num_wells)
    taxon_pure_med, taxon_pure_low, taxon_pure_high = get_CI(taxon_pure_wells, as_string=False)
    return taxon_pure_med, taxon_pure_low, taxon_pure_high


def process_row(row):
    positive_wells = np.zeros(number_of_experiments)
    pure_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)


    for i in range(number_of_experiments):
        positive_wells[i], pure_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(row.cells_per_well, row.rel_abund,num_wells=row.num_wells_inoculated)

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
    pure_med, pure_low, pure_high = get_CI(pure_wells, as_string=False)
    return pd.Series([row.ID, row.Site, row.num_wells_inoculated,
                      positive_med, positive_low, positive_high,
                      taxon_positive_med, taxon_positive_low, taxon_positive_high,
                      pure_med, pure_low, pure_high,
                      taxon_pure_med, taxon_pure_low, taxon_pure_high, within_estimates, deviance, p_value],
                     index=['ID', 'Site', 'num_inoculated_wells',
                            'positive_well_med', 'positive_well_95pc_low','positive_well_95pc_high',
                            'taxon_positive_med', 'taxon_positive_95pc_low', 'taxon_positive_95pc_high',
                                'pure_well_med', 'pure_well_95pc_low', 'pure_well_95pc_high',
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
        taxon_pure_med, taxon_pure_low, taxon_pure_high = bootstrap_dte(cells_per_well=job[1],
                                                                        rel_abund=p,
                                                                        viability=1,
                                                                        num_wells=job[0])
        if (taxon_pure_low >= 1):
            return((job[0], job[1], p, taxon_pure_low))

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
        taxon_pure_med, taxon_pure_low, taxon_pure_high = bootstrap_dte(cells_per_well, rel_abund, num_wells= 9999)
        results.append((cells_per_well, rel_abund, taxon_pure_med, taxon_pure_low, taxon_pure_high))

    results_df = pd.DataFrame(results, columns=['cells_per_well', 'relative_abundance', 'taxon_pure_med', 'taxon_pure_low', 'taxon_pure_high'])
    results_df.to_csv(f'./data/inoculum_rel_abundance_interaction_plot_data_{number_of_experiments}_simulations.csv', index=False)


def estimate_viability(cells_per_well, rel_abund, observed_pure, num_wells, test_min=0, test_max=1, test_step=0.1):
    viability_range = np.arange(start=test_min, stop=test_max, step=test_step)[::-1]
    viability_values = []
    for v in viability_range:
        print(f'*** Testing viability {v:.4f} ***')
        taxon_pure_med, taxon_pure_low, taxon_pure_high = bootstrap_dte(cells_per_well, rel_abund, v, num_wells)
        if observed_pure >=taxon_pure_low and observed_pure <= taxon_pure_high:
            viability_values.append(v)
    try:
        min = np.min(viability_values)
    except ValueError:
        min = 0
    try:
        max = np.max(viability_values)
    except ValueError:
        max = 1

    if test_step==0.1:
        print(f'*** Refining between {np.max((min-test_step, 0)):.4f} and {np.min((max+test_step,1)):.4f}***')
        min, max = estimate_viability(cells_per_well, rel_abund, observed_pure, num_wells,
                    test_min=np.max((np.min(viability_values)-test_step, 0)),
                    test_max=np.min((np.max(viability_values)+test_step,1)),
                    test_step=0.01)

    return min, max


def process_viability(row):
    positive_wells = np.zeros(number_of_experiments)
    pure_wells = np.zeros(number_of_experiments)
    taxon_positive_wells = np.zeros(number_of_experiments)
    taxon_pure_wells = np.zeros(number_of_experiments)

    viability_range = np.arange(start=0, stop=1, step=0.01)[::-1]
    max_viability = 0

    for v in viability_range:
        print(f'*** Processing ID:{row.ID} at site: {row.Site}, testing viability {v:.2f} ***')

        for i in range(number_of_experiments):
            positive_wells[i], pure_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(row.cells_per_well, row.rel_abund, viability=v, num_wells=row.num_wells_inoculated)

        taxon_pure_med, taxon_pure_low, taxon_pure_high = get_CI(taxon_pure_wells, as_string=False)
        if row.num_isolates >=taxon_pure_low and row.num_isolates <= taxon_pure_high:
            max_viability = v
            break
    return pd.Series([row.ID, row.Site, row.rel_abund, row.num_wells_inoculated,
                    row.num_isolates, row.cells_per_well,
                    row.taxon_pure_95pc_low, max_viability], index=['ID', 'Site',
                                                                'rel_abund', 'num_inoculated_wells',
                                                                'num_isolates', 'cells_per_well', 'taxon_pure_95pc_low', 'max_viability'])



def estimate_viability_range_for_underrepresented_taxa():
    print('*** Calculating maximum viability to explain observed number of wells ***')

    df = pd.read_csv(f'./data/ASV_cultivation_bootstrapped_numbers_max_viability_{number_of_experiments}_simulations.joined.csv')
    lower_than_expected_df = df[df.deviance <0]
    print(lower_than_expected_df.shape)
    tqdm_pandas(tqdm())
    results = lower_than_expected_df.progress_apply(process_viability, axis=1)
    results.to_csv(f'./data/estimate_max_viability_range_for_underrepresented_taxa_{number_of_experiments}_simulations.csv', index=False)

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


def main():
    global number_of_experiments
    number_of_experiments = 9999
    #run_bootstrap_per_ASV('./data/ASV_cultivation_numbers.csv')
    #estimate_required_proportion_for_underepresentation()
    #calculate_interaction_rel_abundance_inoculum()
    #estimate_viability_range_for_underrepresented_taxa()
    #make_inoculum_rel_abundance_interaction_plot()
if __name__=="__main__":
    main()
