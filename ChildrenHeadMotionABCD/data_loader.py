import os
import numpy as np
import pandas as pd


### NOTICE ### 
# Docstrings have been generated using ChatGPT to speed-up the book-keeping.


def filter_subj(participants: pd.DataFrame, subjects: pd.DataFrame) -> pd.DataFrame:
    """
    Filter and clean-up subjects based on specified criteria.

    Parameters:
    - participants (pd.DataFrame): DataFrame containing information about all participants.
    - subjects (pd.DataFrame): DataFrame containing a list of participant IDs.

    Returns:
    - pd.DataFrame: Filtered DataFrame containing subjects that meet the specified criteria.
    
    Criteria:
    - Keep only subjects present in the provided list.
    - Retain only data from the 'ses-baselineYear1Arm1' session.
    - Remove subjects with age code 888 (indicating missing age) and convert age to years.

    Example:
    filtered_data = filter_subj(participants_df, subjects_df)
    """

    # Clean-up by baseline year:
    subjects = participants[participants['participant_id'].isin(subjects['participant_id'])]    
    subjects = subjects[subjects['session_id'].isin(['ses-baselineYear1Arm1'])]

    # Removed sub-NDARINV66PM75JX due to the lack of sufficient scans:
    subjects = subjects[subjects['participant_id'] != 'sub-NDARINV66PM75JX']
        
    # Clean-up by age:
    subjects        = subjects[subjects["age"] != 888]
    subjects["age"] = subjects["age"] // 12

    return subjects



def read_data(modality: str) -> pd.DataFrame:
    """
    Reads and preprocess data for a specified modality.

    Parameters:
    - modality (str): Modality type, either 'dti' or 'fmri'.

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing participants and subjects DataFrames.

    Example:
    participants_df, subjects_df = read_data('dti')
    """

    # Create participants dataframe:
    part_path    = os.path.join("data", "descriptions", "participants.tsv")
    participants = pd.read_csv(part_path, delimiter='\t')

    # Create a subject DataFrame based on modality:
    if modality == 'dti':
        subj_path = os.path.join("data", "dti")
    elif modality == 'fmri':
        subj_path = os.path.join("data", "fmri", "resting_state", "fmriresults01", "derivatives", "abcd-hcp-pipeline")
    else:
        raise ValueError("Invalid modality. Supported modalities: 'dti', 'fmri'.")

    # Get a list of all subjects and put them in a DataFrame:
    subjects = os.listdir(subj_path)
    subjects = pd.DataFrame([subject.split('_')[0] for subject in subjects], columns=['participant_id'])



    # Filter and clean-up subjects based on criteria
    subjects = filter_subj(participants, subjects)
    subjects = find_greatest_displacement(subjects, subj_path, modality)

    return subjects 


def sex_split(subjects: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame) :
    """
    Split the subjects DataFrame by sex into separate groups.

    Parameters:
    - subjects (pd.DataFrame): DataFrame containing subject information.

    Returns:
    - Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames, one for male subjects and one for female subjects.

    Example:
    male_group, female_group = sex_split(subjects_df)
    """

    # Split by sex:
    group_sex = subjects.groupby('sex')

    # Define separate groups:
    male = group_sex.get_group(1)
    female = group_sex.get_group(2)

    return male, female


def load_data(path: str, modality: str) -> np.ndarray:
    """
    Load data from a file and extract relevant columns based on the specified modality.

    Parameters:
    - path (str): Path to the data file.
    - modality (str): Modality of the data ('dti' or 'fmri').

    Returns:
    - np.ndarray: Loaded data as a NumPy array.

    Example:
    dti_data = load_data('path/to/dti_data.txt', 'dti')
    fmri_data = load_data('path/to/fmri_data.txt', 'fmri')
    """

    data = []

    # Read data from file:
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if i > 0:
                data.append(line.split())

    data = np.array(data)

    # Extract relevant columns based on modality:
    if modality == 'dti':
        return np.array(data[:, 1:7], dtype=np.double)
    elif modality == 'fmri':
        return np.array(data[:, 0:6], dtype=np.double)
    else:
        raise ValueError("Invalid modality. Supported modalities: 'dti', 'fmri'.")


def find_greatest(data: np.ndarray) -> list:
    """
    Find the greatest absolute value for each column in the input data.

    Note: columns in the case of the ABCD dataset will stand for translations along
          rotations around the X, Y, and Z axis.

    Parameters:
    - data (np.ndarray): Input data as a NumPy array.

    Returns:
    - list: List of greatest values for each column.

    Example:
    greatest_values = find_greatest(my_data)
    """

    _, N = data.shape
    mins = np.min(data, axis=0)
    maxs = np.max(data, axis=0)

    # Choose the greater of the absolute values for each column:
    return [maxs[i] if abs(mins[i]) < abs(maxs[i]) else mins[i] for i in range(N)]


def find_greatest_displacement(subset: pd.DataFrame, subj_path: str, modality: str) -> np.ndarray:
    """
    Find the greatest displacement values for each subject.
    Parameters:
    - subset (pd.DataFrame): Subset of subjects containing participant_id.
    - subj_path (str): Path to the data.
    - modality (str): Modality of the data ('dti' or 'fmri').
    Returns:
    - np.ndarray: Array of greatest displacement values for each subject.
    """
    greatest_values = []
    if modality == 'dti':
        # For DTI modality, iterate through subjects and find greatest displacement from a single file
        for subject in subset['participant_id']:
            path = os.path.join(subj_path, subject + "_ses-baselineYear1Arm1_confounds.tsv")
            data = load_data(path, modality)
            greatest_values.append(find_greatest(data))
    stored_runs = []
    run_nums    = []
    if modality == 'fmri':
        # For fMRI modality, iterate through subjects and find greatest displacement across multiple runs
        for subject in subset['participant_id']:
            runs_dir = os.path.join(subj_path, subject, "ses-baselineYear1Arm1", "func")
            N = len(os.listdir(runs_dir))
            N = (N - 1) // 2
            runs = []

        

            for run_nr in range(1, N + 1):



                path = os.path.join(runs_dir, subject + "_ses-baselineYear1Arm1_task-rest_run-" + str(run_nr) + "_desc-includingFD_motion.tsv")
                data = load_data(path, modality)
                greatest_in_run = find_greatest(data)

                runs.append(greatest_in_run)

            runs = np.array(runs)
            
            # Change back to runs if you want all runs to be considered: (runs[:4])
            stored_runs.append(np.argmax(runs, axis=0) + 1)

            greatest_across_runs = find_greatest(runs)
            greatest_values.append(greatest_across_runs)
            run_nums.append(N)
    stored_runs = np.array(stored_runs)
    greatest_values = np.array(greatest_values)
    if modality == 'dti':
        greatest_values[:, 3:6] = np.rad2deg(greatest_values[:, 3:6])

    subset['max_trans_x_mm'] = greatest_values[:, 0]
    subset['max_trans_y_mm'] = greatest_values[:, 1]
    subset['max_trans_z_mm'] = greatest_values[:, 2]
    subset['max_rot_x_deg'] = greatest_values[:, 3]
    subset['max_rot_y_deg'] = greatest_values[:, 4]
    subset['max_rot_z_deg'] = greatest_values[:, 5]
    
    """subset['run_max_trans_x_mm'] = stored_runs[:, 0]
    subset['run_max_trans_y_mm'] = stored_runs[:, 1]
    subset['run_max_trans_z_mm'] = stored_runs[:, 2]
    subset['run_max_rot_x_deg'] = stored_runs[:, 3]
    subset['run_max_rot_y_deg'] = stored_runs[:, 4]
    subset['run_max_rot_z_deg'] = stored_runs[:, 5]

    subset['runs'] = run_nums"""

    return subset