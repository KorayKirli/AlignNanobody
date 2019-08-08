### Extract Nanobody Sequences from multiple DNA sequencing results (fasta files) and align them

If you have many plasmid dna sequencing results, ie from nanobody generation, or mutagenesis experiments, this tool can help you. I set this up to work with mybinder.org, but feel free to clone the repo and use it locally. 

#### What do you need to use this?
You will need fasta files of your dna sequencing results, put them in a folder, name it as you wish, and then zip the folder. 

#### What does the tool do?
1. Extracting protein sequences
You will need to change the folder name to the name you gave to your folder. When you run the first cell:
It will look for all 6 frames, and extract the longest protein starting with Methionine. You can set a size limit to ignore small ones. All will be stored in a fasta file, which you can find under /data/
2.  Second cell will perform the multiple sequence alignment, it uses Muscle, should work for mac, windows, linux environments
Again, result it stored under /data/
3. I added a visualization part, might not work for all cases, optimized for proteins under 200 amino acids.

I am happy to make changes, add more functions, please use Github Issues for requests.


#### About mybinder

Binder is a service to take a Github repo with Jupyter Notebooks and host it (and its requirements) on a server.

Start this repo on binder, then upload your folder with fasta files under data.
In the notebook update the folder name, and run cells.

You can also share the binder version of your repo as a link. The website is https://mybinder.org/.
