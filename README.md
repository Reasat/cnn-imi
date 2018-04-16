# cnn-imi
This is a repository for the code developed to produce the results in the paper ["Detection of Inferior Myocardial Infarction using Shallow Convolutional Neural Networks"](https://arxiv.org/abs/1710.01115v2)


## Abstract
Myocardial Infarction is one of the leading causes of death worldwide. This paper presents a Convolutional Neural Network (CNN) architecture which takes raw Electrocardiography (ECG) signal from lead II, III and AVF and differentiates between inferior myocardial infarction (IMI) and healthy signals. The performance of the model is evaluated on IMI and healthy signals obtained from Physikalisch-Technische Bundesanstalt (PTB) database. A subject-oriented approach is taken to comprehend the generalization capability of the model and compared with the current state of the art. In a subject-oriented approach, the network is tested on one patient and trained on rest of the patients. Our model achieved a superior metrics scores (accuracy= 84.54%, sensitivity= 85.33% and specificity= 84.09%) when compared to the benchmark. We also analyzed the discriminating strength of the features extracted by the convolutional layers by means of geometric separability index and euclidean distance and compared it with the benchmark model.

While writing the codes, files and folder was organized in the following way
     
#### Project tree
 * PTB
   * code
   * data
   * ptbdb
        
All the code files (.ipynb,.py) were placed in `code` folder. The ecg records were downloaded from [PhysioNet](https://www.physionet.org/physiobank/database/ptbdb/) to the `ptbdb` folder. The preprocessed data, extracted features were saved in `data` folder.

`preprocess_and_segment_data.ipynb` --> Preprocesses the ECG signals and segments them according to [1]

`build_train_validate_cnn.ipynb` --> Builds the convolutional network, trains on training data and evaluates the model's performance on the validation data.

`extract_features_swt.py` --> Extracts feature from ECG signals as described in [1]. These features are later used to calculate geometric separability index and Euclidean distance calculation.

[1] [Sharma, Lakhan Dev, and Ramesh Kumar Sunkaria. "Inferior myocardial infarction detection using stationary wavelet transform and machine learning approach." Signal, Image and Video Processing (2017): 1-8.](https://link.springer.com/article/10.1007/s11760-017-1146-z)
