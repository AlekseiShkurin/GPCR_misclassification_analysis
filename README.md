# Database Curation Tool for Misclassification Analysis

The goal of this project was to create a tool capable of identifying mislabelled samples in a dataset. The project was done using G-Protein-Coupled Receptor (GPCR) database.  

Out of 7 possible classes of GPCRs, 3 have been previously shown to have a high degree of misclassification of samples between those classes.  

This project uses Random Forests (RF) to identify mislabelled samples. Due to its ensemble nature, RF is naturally suited for solving such tasks given we can analyse the collective decision and see the distributions of assignments coming from individual trees.  

As an example, the figure below shows consistency scores for samples in Class 1 that is expected to have a low degree of misclassification. As we can see, for most of the sample the trees are in agreement, assigning samples their right class.  

<img src="Plots/Consistency%20Class%201%2C%20Original.png" width=600>  

Unlike Class 1, consistency measurements of Class 5 are substantially lower. Class 5 was expected to have high misclassification rate with Classes 6 and 7.  

<img src="Plots/Consistency%20Class%205%2C%20Original.png" width=600>  
