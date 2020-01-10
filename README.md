# mtk_primer_generator
Automagically ✨ generate primers for Golden Gate cloning using the Mammalian Cloning Toolkit (MTK)


This python script is run through the command line and prompts the user to input a DNA sequence for domestication, followed by the part type and a short prefix for the output ordering/upload forms. It then outputs a list of primers that can be used to domesticate the input sequence along with the expected product size(s). The primers may be slightly longer than ideal if there is homology between the adapter sequence at the beginning of the primer and the plasmid itself, so we recommend checking primer specificity and melting temperatures using an application like Benchling and adjusting primer length as necessary.

See the publication here for details about the MTK (https://pubs.acs.org/doi/10.1021/acssynbio.9b00322)
