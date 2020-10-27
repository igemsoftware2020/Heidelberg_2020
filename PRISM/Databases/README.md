# Preprocessing Pipeline

The proposed pipeline prepares the sequences and annotations stored in several databases for usage in the Neural Network. Data from several databases, such as UniProtKB Reviewwed, BindingDB Gene Ontology and ATtRACT database  are combined and the sequences and annotation vectors are prepared. A description of entries in the annotation_vectors is provided [here](https://github.com/igemsoftware2020/Heidelberg_2020/blob/main/PRISM/Databases/data/annotations_vector_description.csv). The sequences are padded to a length of 512 amino acids and a mask is is prepared.

Please find detailed information on the [Wiki of iGEM Heidelberg 2020](https://2020.igem.org/Team:Heidelberg/Software/PRISM).

# Data preparation.

Run ```python preprocessing_pipeline.py``` to run the preprocessing pipeline. The databases are downloaded automatically, so the first initialization may take some minutes.

# Data availability

The preprocessing_pipeline.py gives you the option to redo with a newer dataset and/ or modify the data preparation. The preprocessing pipeline also offers the option to load the pre-calculated GO-Terms, GO-Terms mapped for Semantic Similarity as well as taxonomies for the generation of the annotations vectors. To access these data files, just uncomment them in the code. Please note, that the pre-calculated values are recalculated when the preprocessing_pipeline.py is run.  A text list of the [GO-Terms](https://github.com/igemsoftware2020/Heidelberg_2020/blob/main/PRISM/Databases/data/go_terms_in_neural_network.txt), as well as [taxonomies](https://github.com/igemsoftware2020/Heidelberg_2020/blob/main/PRISM/Databases/data/taxonomies_in_neural_network.txt is also accessible.

The original data prepared by the preprocessing pipeline for the Neural Network is also made available via [Zenodo](LINK).
