regex =[ r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<time_between_factor_1_and_factor_2>\w+\s+\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<time_between_factor_1_and_factor_2>\w+\s+\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<non_genetic_factor_type_3>\w+\s+\w+)\s+(?P<factor_length_3>\d+)\s+(?P<factor_length_type_3>\w+)\s+(?P<factor_dose_3>\w+)\s+(?P<factor_dose_unit_3>\w+)\s+(?P<time_between_factor_3_and_4>\"\w+\"\s+\w+\s+\w+)\s+(?P<number_of_days>\d)\s+(?P<non_genetic_factor_type_4>\w+\s\w+\-\w+)\s+(?P<factor_length_4>\d+)\s+(?P<factor_length_type_4>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<time_between_factor_1_and_factor_2>\w+\s+\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<non_genetic_factor_type_3>\w+\s+\w+)\s+(?P<factor_length_3>\d+)\s+(?P<factor_length_type_3>\w+)\s+(?P<factor_dose_3>\d+\.\d+)\s+(?P<factor_dose_unit_3>\w+\/\w+)\s+(?P<non_genetic_factor_type_4>\w+\s+\w+)\s+(?P<factor_length_4>\w+)\s+(?P<factor_length_type_4>\w+)\s+(?P<factor_dose_4>\w+)\s+(?P<factor_dose_unit_4>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+(?P<other_exposure>\w+)",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<non_genetic_factor_type_3>\w+\s+\w+)\s+(?P<factor_length_3>\d+)\s+(?P<factor_length_type_3>\w+)\s+(?P<factor_dose_3>\d+)\s+(?P<factor_dose_unit_3>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<time_between_factor_1_and_factor_2>\w+\s+\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<non_genetic_factor_type_3>\w+\s+\w+)\s+(?P<factor_length_3>\d+)\s+(?P<factor_length_type_3>\w+)\s+(?P<factor_dose_3>\d+\.\d+)\s+(?P<factor_dose_unit_3>\w+\/\w+)\s+",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>\w+)\s+(?P<time_between_factor_1_and_factor_2>\w+\s+\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<no_specified_exposure>[A-Z]+\s+[A-Z]+\s+[A-Z]+)",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+.\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>.\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s+\w+\s+\w+\:\s+\d+)\s+(?P<factor_length_2>\d)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+\.\d+\/\d+\.\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<non_genetic_factor_type_3>\w+\s+\w+)\s+(?P<factor_length_3>\d+)\s+(?P<factor_length_type_3>\w+)\s+(?P<factor_dose_3>\d+)\s+(?P<factor_dose_unit_3>\w+)\s+(?P<time_between_factor_3_and_factor_4>\"\w+\"\s+\w+\s\w+)\s+(?P<number_of_days>\d+)\s+(?P<non_genetic_factor_type_4>\w+.\w+\-\w+)\s+(?P<factor_length_4>\d)\s+(?P<factor_length_type_4>\w+)",

r"\s+(?P<ovarian_stimulation_type>\w+)\s+(?P<agonist_name>\w+\s+\w+\s+\(\w+\))\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+\s+\w+\s+\w+-\w+)\s+(?P<model_of>\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+\_\w+\_\w+)\s+(?P<non_genetic_factor_type_1>\w+\-\d+\s+\(\w+\))\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+\_\w+\_\w+)\s+(?P<non_genetic_factor_type_1>\w+\-\w+\-\w+)\s+(?P<factor_length>\d+)\s+(?P<factor_length_type>\w+)\s+(?P<factor_dose>\d+)\s+(?P<factor_dose_unit>\w+\/\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+(?P<no_specified_exposure>[A-Z]+\s+[A-Z]+\s+[A-Z]+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+\_\w+\_\w+)\s+(?P<no_specified_exposure>[A-Z]+\s+[A-Z]+\s+[A-Z]+)",

r"(?P<cycle_pregnancy_phase>\w+\s+\w+\/\w+\s+\w+)\s+(?P<cycle_pregnancy_days_1>\w+\s+\d+)\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\(\w+\:\w+\))\s+(?P<factor_length>\d+)\s+(?P<factor_length_type>\w+)\s+(?P<factor_dose>\w+)\s+(?P<factor_dose_unit>(\w+|\w+\/\w+))\s+$",

r"(?P<cycle_pregnancy_phase>\w+\s+\w+\/\w+\s+\w+)\s+(?P<cycle_pregnancy_days_1>\w+\s+\d+)\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+(?P<no_specified_exposure>[A-Z]+\s+[A-Z]+\s+[A-Z]+)\s+$",

r"(?P<cycle_pregnancy_phase>\w+\s+\w+\/\w+\s+\w+)\s+(?P<cycle_pregnancy_days_1>\w+\s+\d+)\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\(\w+\:\w+\))\s+(?P<factor_length>\d+)\s+(?P<factor_length_type>\w+)\s+(?P<factor_dose>\w+)\s+(?P<factor_dose_unit>(\w+|\w+\/\w+))\s+(?P<time_between_factor_1_and_factor_2>\"\w+\"\s+\w+\s+\w+)\s+(?P<number_of_days>\d+)\s+(?P<non_genetic_factor_type_2>\w+\s+\(\w+\:\w+\))\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+$",

r"(?P<cycle_pregnancy_phase>\w+\s+\w+\/\w+\s+\w+)\s+(?P<cycle_pregnancy_days_1>\w+\s+\d+)\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\(\w+\:\w+\))\s+(?P<factor_length>\d+)\s+(?P<factor_length_type>\w+)\s+(?P<factor_dose>\w+)\s+(?P<factor_dose_unit>(\w+|\w+\/\w+))\s+(?P<time_between_factor_1_and_factor_2>\"\w+\"\s+\w+\s+\w+)\s+(?P<number_of_days>\d+)\s+(?P<non_genetic_factor_type_2>\w+\s+\(\w+\:\w+\))\s+(?P<factor_length_2>\d+)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d+)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<no_specified_exposure>[A-Z]+\s+[A-Z]+\s+[A-Z]+)",

r"(?P<cycle_pregnancy_phase>\w+\s+\w+\/\w+\s+\w+)\s+(?P<cycle_pregnancy_days_1>\w+\s+\d+)\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<established_animal_model>\w+)\s+$",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<non_genetic_factor_type_1>\w+\s+\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+\.\d+)\s+(?P<factor_dose_unit_1>\w+\/\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+(?P<no_specified_exposure>[A-Z]+\s+[A-Z]+\s+[A-Z]+)",

r"\s+(?P<genetic_exposure_type>\w+\_\w+)\s+(?P<animal_model_gene_name>\w+)\s+(?P<animal_model_modification_type>\w+)\s+(?P<het_homo>h\w+)\s+(?P<non_genetic_factor_type_1>\w+.\w+\-\w+)\s+(?P<factor_length_1>\d+)\s+(?P<factor_length_type_1>\w+)\s+(?P<factor_dose_1>\d+)\s+(?P<factor_dose_unit_1>.\w+)\s+(?P<time_between_factor_1_and_factor_2>\w+\s+\w+)\s+(?P<non_genetic_factor_type_2>\w+\;\s\w+\s\w+\:\s\d+)\s+(?P<factor_length_2>\d)\s+(?P<factor_length_type_2>\w+)\s+(?P<factor_dose_2>\d\.\d\/\d\.\d)\s+(?P<factor_dose_unit_2>\w+\/\w+)\s+(?P<ovarian_stimulation_name_mice>\w+\s+\+\s+\w+)\s+(?P<other_exposure>\w+)\s+(?P<no_specified_exposure>[A-Z]+\s[A-Z]+\s[A-Z]+)"

]