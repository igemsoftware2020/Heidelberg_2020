from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import pickle as pk
import sets

def preprocessing_pca(annotation_vectors_train, annotation_vectors_val):
    '''
    Function performs PCA for Feature Reuction to 512 components for all annotation vectors.

    :param annotation_vectors_train: annotation_vectors for training
    :param annotation_vectors_val: annotation_vectors for validation
    :return annotation_vectors_train: feature-reduced annotation_vectors for training
    :return annotation_vectors_val: feature-reduced annotation_vectors for validation
    '''

    try:
        pca = pk.load(open("pca.pkl", 'rb'))
    except:
        pca = pk.load(open("./LanguageModel/pca.pkl", 'rb'))

    annotation_vectors_train = pca.transform(annotation_vectors_train)
    annotation_vectors_val = pca.transform(annotation_vectors_val)

    return annotation_vectors_train, annotation_vectors_val


def preprocessing_pca_initialization(annotation_vectors):
    '''
    Function fits PCA on all annotations vectors for Pretraining and Finetuning.

    :param annotation_vectors: annotation_vectors used in pretraining and finetuning
    '''

    pca = PCA(n_components=512)

    pca.fit(annotation_vectors)

    pk.dump(pca, open("pca.pkl", "wb"))

    columns_loadings = []

    for i in range(0, len(pca.explained_variance_ratio_)):
        columns_loadings.append("PC" + str(i + 1))

    loadings = pd.DataFrame(pca.components_.T, columns= columns_loadings)#, index=annotations_vector_description[0:7]) #annotations_vector_description would need to be imported or loaded! For real annotation vectors delete [0:7]
    print(loadings)


if __name__ == "__main__":

    size_val_set = 0.1
    protein_families_pretraining = np.load("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/training/protein_families.npy", allow_pickle=True)
    protein_families_finetuning = np.load("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/protein_families.npy", allow_pickle=True)
    count_protein_family_list_pretraining = sets.determine_protein_family_lengths(protein_families_pretraining)
    count_protein_family_list_finetuning = sets.determine_protein_family_lengths(protein_families_finetuning)

    dataset_annotation_vector_pretraining = np.load("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/training/annotation_vectors.npy",allow_pickle=True)
    dataset_annotation_vector_finetuning = np.load("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/rbp_annotation_vectors.npy",allow_pickle=True)

    dataset_annotation_vector_pretraining = sets.cross_validation_dataset_split(dataset_annotation_vector_pretraining, count_protein_family_list_pretraining, (1 / size_val_set))
    dataset_annotation_vector_finetuning = sets.cross_validation_dataset_split(dataset_annotation_vector_finetuning, count_protein_family_list_finetuning, (1 / size_val_set))

    dataset_annotation_vector_pretraining, _ = sets.cross_validation_dataset_combination(dataset_annotation_vector_pretraining, 0)
    dataset_annotation_vector_finetuning, _ = sets.cross_validation_dataset_combination(dataset_annotation_vector_finetuning, 0)

    annotation_vectors = np.concatenate((dataset_annotation_vector_pretraining, dataset_annotation_vector_finetuning), axis=0)

    preprocessing_pca_initialization(annotation_vectors)



