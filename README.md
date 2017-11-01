# cnn-imi
This is a repository for the code developed to produce the results in the paper "Detection of Inferior Myocardial Infarction using Shallow Convolutional Neural Networks" (https://arxiv.org/abs/1710.01115v2)

While writing the codes, files and folder was organized in the following way
     
#### Project tree
 * PTB
   * code
   * data
   * ptbdb
        
All the code files (.ipynb,.py) were placed in code folder. The ecg records were downloaded to ptbdb folder. The preprocessed data, extracted features were saved in data folder.

preprocess_and_segment_data.ipynb --> Preprocesses the ECG signals and segments them according to [1]

build_and_train_cnn.ipynb --> Builds the convolutional network, trains on training data and evaluates the model's performance on the validation data.

extract_features_swt.py --> Extracts feature from ECG signals as described in [1]

[1] "Inferior myocardial infarction detection using stationary wavelet transform and machine learning approach" (https://link.springer.com/article/10.1007/s11760-017-1146-z)
