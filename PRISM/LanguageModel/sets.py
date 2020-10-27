import math
import numpy as np
from sklearn.model_selection import train_test_split

def cross_validation_dataset_exact_split(dataset, k_fold):
    """
    Function splits the dataset in k_fold parts.
    :param dataset: complete dataset
    :param k_fold: number of parts the complete dataset is supposed to be split in
    :return dataset_split: dataset splitted in k_fold parts
    """
    multiplier = divmod(len(dataset), k_fold)
    multiplier_loop = multiplier[0]
    dataset_split = []
    rest = 0

    for set in range(0, k_fold):
        if set == k_fold-1:
            rest = multiplier[1]
        dataset_split.append(dataset[0+multiplier_loop*set:multiplier_loop+multiplier_loop*set+rest])

    return dataset_split

def cross_validation_dataset_split(dataset, count_protein_family_list, k_fold):
    """
    Function splits the dataset in k_fold parts with respect to borders defined by protein_families, so sequences of one protein family are not separated into different entries in dataset_split.

    :param dataset: complete dataset
    :param count_protein_family_list: protein families corresponding to entries in dataset
    :param k_fold: number of parts the complete dataset is supposed to be split in
    :return dataset_split: dataset splitted in k_fold parts with respect of lengths of protein_families defined by protein_families_length
    """

    multiplier = divmod(len(dataset), k_fold)
    multiplier_loop = multiplier[0]
    dataset_split = []
    counter = 0

    dataset_split_part = 0
    dataset_split_former_parts = 0

    for element in count_protein_family_list:
        dataset_split_part = dataset_split_part + element
        counter = counter + 1
        if counter == len(count_protein_family_list):
            dataset_split.append(
                dataset[0 + dataset_split_former_parts:dataset_split_former_parts + dataset_split_part])
        elif dataset_split_part > multiplier_loop:
            if abs(multiplier_loop - dataset_split_part) < abs(multiplier_loop - (dataset_split_part - element)):
                dataset_split.append(
                    dataset[0 + dataset_split_former_parts:dataset_split_former_parts + dataset_split_part])
                dataset_split_former_parts = dataset_split_former_parts + dataset_split_part
                dataset_split_part = 0
            else:
                dataset_split.append(
                    dataset[0 + dataset_split_former_parts:dataset_split_former_parts + dataset_split_part - element])
                dataset_split_former_parts = dataset_split_former_parts + dataset_split_part - element
                dataset_split_part = element

    return dataset_split


def cross_validation_dataset_combination(dataset_split, set):
    """
    Function combines dataset_split into the validation_set and the training_set by choosing with set the validation_set.
    
    :param dataset_split: the complete dataset splitted in equally big parts
    :param set: position of validation dataset in dataset_split
    :return validation_set: 1/10th of the dataset used
    :return training_set: rest of the dataset
    """
    validation_set = dataset_split[set]

    training_set = np.array([])

    for dataset_set in dataset_split[0:set]:
        if training_set.size == 0:
            training_set = dataset_set
        else:
            training_set = np.concatenate((training_set, dataset_set))

    for dataset_set in dataset_split[set+1:]:
        if training_set.size == 0:
            training_set = dataset_set
        else:
            training_set = np.concatenate((training_set, dataset_set))

    return training_set, validation_set

def finetuning_dataset(dataset):
    """
    Function separates the dataset in the training_set, validation_set and test_set.

    :param dataset: complete dataset
    :return training_set: 8/10th of dataset for training
    :return validation_set: 1/10th of dataset for validation
    :return test_set: 1/10th of dataset for testing
    """

    #MISCHEN?

    limit = math.ceil(len(dataset) * 0.1)

    validation_set = dataset[:limit]
    test_set = dataset[limit:limit*2]
    training_set = dataset[limit*2:]

    return training_set, validation_set, test_set

def determine_protein_family_borders(protein_families):

    protein_families_borders = []

    for element in range(0, len(protein_families) - 1):
        if protein_families[element] != protein_families[element + 1]:
            protein_families_borders.append(element)

    return protein_families_borders

def determine_protein_family_lengths(protein_families):

    count_protein_family = 1
    count_protein_family_list = []

    for element in range(0, len(protein_families)):
        if element != len(protein_families) - 1:
            if protein_families[element] == protein_families[element + 1]:
                count_protein_family = count_protein_family + 1
            else:
                count_protein_family_list.append(count_protein_family)
                count_protein_family = 1
        if element == (len(protein_families) - 1):
            count_protein_family_list.append(count_protein_family)

    return count_protein_family_list

if __name__ == "__main__":
    train_texts, val_texts, train_labels, val_labels = train_test_split(train_texts, train_labels, test_size=.2) #EXAMPLE
