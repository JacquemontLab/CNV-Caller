## Entry Workflow Design

The following DAG (Directed Acyclic Graph) is meant to illustrate the main nextflow workflow logic. 



![CNV_Calling_dag](https://github.com/user-attachments/assets/9a0834aa-1ab2-4cb2-b8a5-b0ef3dcc2d4b)


#### Legend:
Arrows
  - Animated: queue channel at sample level
  - Straight: value channel
  - Line: string interpolation

Colours:
  - Green: Process/Module containing several processes
  - Grey: Nextflow functions. File manipulations operated within the workflow.
  - Blue: Directory defined either with string globbing or just the filepath
