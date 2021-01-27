import csv
from pathlib import Path
import re



def process_score_file_line(line):
    '''Process a Rosetta score file line into a list.

    Args:
        line (str): Line from Rosetta score file.

    Returns:
        list: List of items in the line string.
    '''
    l = [i.strip() for i in re.split(' +', line)[1:] if i.strip()]
    for field in l:
        assert l != None
    return l


def read_score_file(filepath):
    '''Read a Rosetta score file.

    Args:
        filepath (str or Path): Path to score file.

    Returns:
        tuple: tuple of two items, the header (column names) as a list and the
        score file rows as a list of lists.
    '''
    with open(str(filepath)) as handle:
        next(handle)  # skip the first line dont need it
        header = process_score_file_line(next(handle))
        lines = []
        
        for line in handle.readlines():
            line = process_score_file_line(line)
            lines.append(line)
            

    return header, lines


def add_string_to_score_file_lines(score_file_header, score_file_lines,
                                   name, string):
    '''Add a constant string to the lines read from a Rosetta score file
    (see read_score_file) and append name to the header.

    Args:
        score_file_header (list): Header for the score_file_lines.
        score_file_lines (list): List of lists holding fields read from Rosetta.
                                 score file.
        name (str): Name to be used as the column name for the string.
        string (str): String to add to all rows.

    Returns:
        [type]: [description]
    '''
    score_file_header.append(name)
    for i, line in enumerate(score_file_lines):
        score_file_lines[i].append(string)
    # score_file_lines = [print(line) for line in score_file_lines]
    return score_file_header, score_file_lines


def add_result_pdb_to_score_file_lines(score_file_header, score_file_lines, 
                                       score_file_path):
    name = 'result_pdb'
    result_dir = score_file_path.parent
    descrip_column = score_file_header.index('description')
    for i, line in enumerate(score_file_lines):
        descrip = line[descrip_column]
        pdb_path = Path(result_dir).joinpath(Path(descrip).with_suffix('.pdb'))
        assert pdb_path.is_file()
        score_file_lines[i].append(pdb_path)
    score_file_header.append('pdb_path')

    return score_file_header, score_file_lines


def add_multi_target_protein_identifiers(target_protein_name, score_file_header, 
                                         score_file_lines, target_protein_dir=None):
    score_file_header, score_file_lines = add_string_to_score_file_lines(
        score_file_header, score_file_lines,
        'target_protein', str(target_protein_name)
    )
    if target_protein_dir:
        assert Path(target_protein_dir).is_dir()
        target_protein_path = Path(
            target_protein_dir).joinpath(target_protein_name)
        score_file_header, score_file_lines = add_string_to_score_file_lines(
            score_file_header, score_file_lines,
            'target_protein_path', str(target_protein_path)
        )

    return score_file_header, score_file_lines


def add_observation_identifiers(results_dir, score_file_path, score_file_header, score_file_lines, target_protein_dir=None):
    distance_to_parent = 0
    results_dir = Path(results_dir)
    score_file_path = Path(score_file_path)
    score_file_path_copy = score_file_path
    # print(len(score_file_lines), 'score file lines')
    while score_file_path_copy != results_dir:
        score_file_path_copy = score_file_path_copy.parent
        distance_to_parent += 1

    if distance_to_parent == 4:  # this should be a job submitted for multiple target proteins
        target_protein = list(score_file_path.parents)[2].name
        score_file_header, score_file_lines = add_multi_target_protein_identifiers(
             target_protein, score_file_header, score_file_lines, target_protein_dir
         )
        print(len(score_file_lines), 'score_file lines')
    ligand_name = score_file_path.parents[1].name

    score_file_header, score_file_lines = add_string_to_score_file_lines(
        score_file_header, score_file_lines, 'ligand_name', ligand_name)

    return score_file_header, score_file_lines


def recursive_score_file_collect(parent, score_file_name='score.sc'):
    # print('Collecting score files from {}'.format(parent))
    return Path(parent).rglob("*.sc")


def write_list_as_delim(filepath, header, nested_list, delimiter='\t'):
    with open(str(filepath), 'w') as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(header)
        writer.writerows(nested_list)

    return filepath


def write_formated_score_files(results_dir, target_protein_dir,
                               formated_filename='score.formated.tsv'):
    score_files = recursive_score_file_collect(results_dir)
    formated_score_files = []
    for sf in score_files:
        # print('Formating {}'.format(sf))
        sf_header, sf_lines = read_score_file(sf)
        # print('Read score file', 'Number lines {}'.format(len(sf_lines)))
        sf_header, sf_lines = add_observation_identifiers(
            results_dir, sf, sf_header, sf_lines, target_protein_dir
        )
        # print('added observation ids')
        formated_sf_path = Path(sf).with_name(formated_filename)
        write_list_as_delim(formated_sf_path, sf_header, sf_lines)
        formated_score_files.append(formated_sf_path)

    return formated_score_files


def concatenate_score_files(formated_score_files, target_file):
    with open(str(target_file), 'w') as target_handle:
        for i, formated_score_file in enumerate(formated_score_files):
            with open(str(formated_score_file)) as read_handle:
                header = next(read_handle)
                if i == 0:
                    target_handle.write(header)
                else:
                    for line in read_handle:
                        target_handle.write(line)
    return target_file
