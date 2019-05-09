# Analysis-Neuronal-Data

Scripts for the analysis of neuronal responses acquired using 2-photon Ca2+ imaging

Plotting and analysing correlations of neuronal response size across different recording sessions
and different experimental groups (statistical testing: Fisher's Z, Steiger's Z and bootstrapping): 
- compare_correlation_across_groups.py
- compare_correlation_within_group.py

Performing PCA on neuronal response traces and projecting response trajectories for each stimulus/session
into the principal component space and calculate the distance between the trajectories to investigate changes in the neuronal representation of the stimuli across sessions:
- PCA_neuronal_response_traces.m
