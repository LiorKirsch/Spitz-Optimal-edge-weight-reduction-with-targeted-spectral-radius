Download the networks from koblenz konect dataset and saves them as sparse matlab matrices.
==========================================================================================

`konect_meta.csv` - a csv file which contains the metadata for the networks in koblenz konect.

`download_network.m` - downloads the networks (that are available in 2015). This function uses the unix wget command so it might not work on windows.

`get_network_data` - get a list of networks using specific filters (edge type, number of edges ...).
