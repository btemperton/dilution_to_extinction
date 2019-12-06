import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm, tqdm_pandas
import itertools
from multiprocessing import Pool

def calculate_dte(inoculum_size, proportion_of_taxa, viability=1, num_wells=96):
    inocula=np.random.poisson(inoculum_size, num_wells)
    taxon_wells=np.random.binomial(inocula, proportion_of_taxa)
    viable_taxon_wells=np.random.binomial(taxon_wells, viability)
    positive_wells = np.sum(inocula >=1)
    pure_wells = np.sum(inocula==1)
    #This captures all the wells with greater than 0 cells where the number
    #of cells in it are all from the selected taxon
    positive_wells_with_taxon=np.sum(viable_taxon_wells>=1)
    taxon_pure_wells = np.sum((viable_taxon_wells == inocula) & (inocula >0))
    return positive_wells, pure_wells, positive_wells_with_taxon, taxon_pure_wells

def get_CI(array_of_numbers, as_string=True, ci=95):
  median_value = np.median(array_of_numbers)
  ci_frac = (100 - ci) / 2
  percentiles = np.percentile(array_of_numbers, [ci_frac, 100 - ci_frac])
  if as_string:
    return f'{median_value:.1f} wells ({percentiles[0]:.1f}-{percentiles[1]:.1f}, 95% CI)'
  else:
    return median_value, percentiles[0], percentiles[1]


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
    tqdm.pandas()
    results = df.progress_apply(process_row, axis=1)
    joined_results = df.merge(results, how='left', on=['ID', 'Site'])
    joined_results.to_csv(f'./data/ASV_cultivation_bootstrapped_numbers_max_viability_{number_of_experiments}_simulations.joined.csv', index=False)

def process(job):
    test_proportions = np.arange(start=0.0, stop=1, step=0.001)
    for p in test_proportions:
        positive_wells = np.zeros(number_of_experiments)
        pure_wells = np.zeros(number_of_experiments)
        taxon_positive_wells = np.zeros(number_of_experiments)
        taxon_pure_wells = np.zeros(number_of_experiments)
        for i in range(number_of_experiments):
            positive_wells[i], pure_wells[i], taxon_positive_wells[i], taxon_pure_wells[i] = calculate_dte(job[1], p,num_wells=job[0])
        taxon_pure_med, taxon_pure_low, taxon_pure_high = get_CI(taxon_pure_wells, as_string=False)
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



def main():
    global number_of_experiments
    number_of_experiments = 9999
    run_bootstrap_per_ASV('./data/ASV_cultivation_numbers.csv')
    #estimate_required_proportion_for_underepresentation()

if __name__=="__main__":
    main()
