import pandas as pd
import re
import mygene
import regex_exposure
import regex_observation_code
import regex_subject_baseline

mg = mygene.MyGeneInfo()
annotations_file = 'test_annotation_POI_functional.xlsx'

results = pd.read_excel(annotations_file, sheet_name='RESULTS', index_col=None).reset_index().iloc[1:]
old_col_names = results.columns[1:].tolist()
old_col_names.insert(0, 'index')
new_col_names = results.iloc[0].values.tolist()
new_col_names[1] = 'comparison_type'
for i in range(len(new_col_names)):
    new_col_names[i] = re.sub(r"\s+", "_", new_col_names[i])
y = 1
for x in range(8,22):
    new_col_names[x] = 'series_exposure_code_'+str(y)
    y+=1
col_rename_dict = {i: j for i, j in zip(old_col_names, new_col_names)}
results.rename(columns=col_rename_dict, inplace=True)
results = results.dropna(axis=1, thresh=2)

results_grouped = results.groupby(by='PMID')

regex_exposure_compiled = []
regex_subject_baseline_compiled = []
regex_observation_code_compiled = []


def compile_regex(type):
    if type not in ('regex_exposure', 'regex_observation_code', 'regex_subject_baseline'):
        print('Error. Enter an existing regex file')
    else:
        for line in __import__(type).regex:
            compiled = re.compile(line)
            globals()[type+'_compiled'].append(compiled)


unique_pmids = results_grouped.groups.keys()
annotations_dict = {}
for pmid in unique_pmids:
    annotations_dict[pmid] = {}
    for entry in range(len(results_grouped.get_group(pmid).index)):
        annotations_dict[pmid][entry] = {}
        annotations_dict[pmid][entry] = results_grouped.get_group(pmid).iloc[entry].to_dict()

def make_dict(type):

    if type in ('binary_subject_1_exposure', 'binary_subject_2_exposure', 'single_cohort_2_exposure', 'series_exposure_code_1', 'series_exposure_code_2', 'series_exposure_code_3'):
        for pmid in unique_pmids:
            for entry in range(len(results_grouped.get_group(pmid))):
                for regex in regex_exposure_compiled:
                    match = regex.match(str(annotations_dict[pmid][entry][type]))
                    if match:
                        annotations_dict[pmid][entry][type] = match.groupdict()
    if type in ('binary_subject_1_baseline', 'binary_subject_2_baseline', 'single_cohort_1_baseline'):
        for pmid in unique_pmids:
            for entry in range(len(results_grouped.get_group(pmid))):
                for regex in regex_subject_baseline_compiled:
                    match = regex.match(str(annotations_dict[pmid][entry][type]))
                    if match:
                        annotations_dict[pmid][entry][type] = match.groupdict()
    if type in ('observation_code'):
        for pmid in unique_pmids:
            for entry in range(len(results_grouped.get_group(pmid))):
                for regex in regex_observation_code_compiled:
                    match = regex.match(str(annotations_dict[pmid][entry][type]))
                    if match:
                        annotations_dict[pmid][entry][type] = match.groupdict()

compile_regex('regex_exposure')
compile_regex('regex_observation_code')
compile_regex('regex_subject_baseline')
make_dict('binary_subject_1_baseline')
make_dict('binary_subject_2_baseline')
make_dict('binary_subject_1_exposure')
make_dict('binary_subject_2_exposure')
make_dict('observation_code')

for pmid in unique_pmids:
    for entry in range(len(results_grouped.get_group(pmid))):
        del annotations_dict[pmid][entry]['PMID']
del annotations_dict['PMID']
pmids = []
for pmid in unique_pmids:
    pmids.append(pmid)
pmids.pop(-1)

for pmid in pmids:
    for entry in range(0, len(annotations_dict[pmid])):
        for keys in annotations_dict[pmid][entry].keys():
            if type(annotations_dict[pmid][entry][keys]) not in (str, float):
                for key, value in annotations_dict[pmid][entry][keys].items():
                    if 'C05503\xa017beta-estradiol' in value:
                        annotations_dict[pmid][entry][keys][key] = 'C05503 17beta-estradiol'

functional_dict = {}

## remove from list HERE if some observations are not necessary ###
annotations_results = ['observation_code', 'direction_of_effect', 'qualitative_presence', 'qualitative_trend',
                       'qualitative_significance', 'qualitative_p_value', 'p_value_threshold', 'abnormal']

for pmid in pmids:
    functional_dict[pmid] = {}
    for entry in range(0, len(annotations_dict[pmid])):
        functional_dict[pmid][entry] = {}
        if type(annotations_dict[pmid][entry]['binary_subject_1_baseline']) not in (float,):
            if annotations_dict[pmid][entry]['binary_subject_1_baseline'] != annotations_dict[pmid][entry][
                'binary_subject_2_baseline']:
                all_keys = annotations_dict[pmid][entry]['binary_subject_1_baseline'].keys() | \
                           annotations_dict[pmid][entry]['binary_subject_2_baseline'].keys()
                sb1 = annotations_dict[pmid][entry]['binary_subject_1_baseline']
                sb2 = annotations_dict[pmid][entry]['binary_subject_2_baseline']
                functional_dict[pmid][entry]['diff_baseline'] = {}
                functional_dict[pmid][entry]['diff_baseline']['sb1'] = {}
                functional_dict[pmid][entry]['diff_baseline']['sb2'] = {}
                for key in all_keys:
                    if (key in sb1) and (key in sb2):
                        if sb1[key] != sb2[key]:
                            functional_dict[pmid][entry]['diff_baseline']['sb1'][key] = sb1[key]
                            functional_dict[pmid][entry]['diff_baseline']['sb2'][key] = sb2[key]
                    if (key in sb1) and (key not in sb2):
                        functional_dict[pmid][entry]['diff_baseline']['sb1'][key] = sb1[key]
                    if (key in sb2) and (key not in sb1):
                        functional_dict[pmid][entry]['diff_baseline']['sb2'][key] = sb2[key]
                functional_dict[pmid][entry]['same_baseline'] = {}
                functional_dict[pmid][entry]['same_baseline'].update(
                    (annotations_dict[pmid][entry]['binary_subject_2_baseline'].items()) & (
                        annotations_dict[pmid][entry]['binary_subject_1_baseline'].items()))
            if annotations_dict[pmid][entry]['binary_subject_1_baseline'] == annotations_dict[pmid][entry][
                'binary_subject_2_baseline']:
                functional_dict[pmid][entry]['identical_baseline'] = {}
                functional_dict[pmid][entry]['identical_baseline'].update(
                    annotations_dict[pmid][entry]['binary_subject_1_baseline'])
                ## EXPOSURE ##
        if type(annotations_dict[pmid][entry]['binary_subject_1_exposure']) not in (float,):
            if annotations_dict[pmid][entry]['binary_subject_1_exposure'] != annotations_dict[pmid][entry][
                'binary_subject_2_exposure']:
                all_keys = annotations_dict[pmid][entry]['binary_subject_1_exposure'].keys() | \
                           annotations_dict[pmid][entry]['binary_subject_2_exposure'].keys()
                se1 = annotations_dict[pmid][entry]['binary_subject_1_exposure']
                se2 = annotations_dict[pmid][entry]['binary_subject_2_exposure']
                functional_dict[pmid][entry]['diff_exposure'] = {}
                functional_dict[pmid][entry]['diff_exposure']['se1'] = {}
                functional_dict[pmid][entry]['diff_exposure']['se2'] = {}
                for key in all_keys:
                    if (key in se1) and (key in se2):
                        if se1[key] != se2[key]:
                            functional_dict[pmid][entry]['diff_exposure']['se1'][key] = se1[key]
                            functional_dict[pmid][entry]['diff_exposure']['se2'][key] = se2[key]
                    if (key in se1) and (key not in se2):
                        functional_dict[pmid][entry]['diff_exposure']['se1'][key] = se1[key]
                    if (key in se2) and (key not in se1):
                        functional_dict[pmid][entry]['diff_exposure']['se2'][key] = se2[key]
                functional_dict[pmid][entry]['same_exposure'] = {}
                functional_dict[pmid][entry]['same_exposure'].update(
                    (annotations_dict[pmid][entry]['binary_subject_2_exposure'].items()) & (
                        annotations_dict[pmid][entry]['binary_subject_1_exposure'].items()))
            if annotations_dict[pmid][entry]['binary_subject_1_exposure'] == annotations_dict[pmid][entry][
                'binary_subject_2_exposure']:
                functional_dict[pmid][entry]['identical_exposure'] = {}
                functional_dict[pmid][entry]['identical_exposure'].update(
                    annotations_dict[pmid][entry]['binary_subject_1_exposure'])
        for result in annotations_results:
            if type(annotations_dict[pmid][entry][result]) not in (float,):
                functional_dict[pmid][entry][result] = {}
                functional_dict[pmid][entry][result] = annotations_dict[pmid][entry][result]
        functional_dict[pmid][entry]['comparison_type'] = annotations_dict[pmid][entry]['comparison_type']

