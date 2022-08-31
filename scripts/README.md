1. Run set_to_run_autoencoder.sh on the cgrb server in multiomics/deep_learning folder. This will run through a bunch of different parameters, in this case distance metrics, and train a neural net using each one at a bunch of different laten space sizes
2. Run getBestFiles_autoencoder.py to delete all files except the one with the highest validation accuracy
3. Copy those weight files over to local machine
4. Run getBestModel_autoencoder.py to see which model performs best when looking at the sample x sample distance matrices (horn distance) using the original data and the data encoded into latent space
6. Manually copy the file names of the best performing autoencoder per omic into the script downstream_prediction_autoencoder.py
5. Run downstream_prediction_autoencoder.py to see how that latent space compares with the original space when predicting lifestyle variables