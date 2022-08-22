# Simulation study to assess correlation between fold change and cell viability in high-throughput screening experiments

This web app is a simulation study to assess the correlation between fold change and cell viability in high-throughput screening experiments.

In the drug discovery funnel, we start with a gene of interest and a large number of compounds (oligos) and want to identify more promising oligos to move forward. High-throughput screening is used to test large number of olives to see if any of them increase expression of our gene of interest. 

Assays are used to detect increases in expression, and every experiment includes control oligos that should not affect the expression of the gene, which allows us to interpret results relative to controls. 

For typical screens, protein abundance is considered relative to cell viability to interpret screening results, which are then compared to the average of controls to calculate a fold change. Although there should ideally be no relationship between fold change and cell viability, relationships are seen in differing directions for various genes. This could be a sign of a problem with the assay, but we can use a simulation study to see what correlation to expect when we know there’s not true relationship between cell viability and fold change to help troubleshoot this behavior.

## Code flow

This is an overview of core components/functions:

![alt text](https://github.com/sheryylli/cell-viability/blob/master/images/codeflow.png)
