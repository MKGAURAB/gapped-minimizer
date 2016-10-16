# gapped-minimizer
A novel Technique of Minimzer caculation in reference and Count the k-mers of each read found in Reference genome.

# Installation
git clone https://github.com/ddxofy/gapped-minimizer.git && cd gapped-minimizer<br />
g++ -Wall -std=c++11 -pthread Main.cpp -lz

#Parameters
./Main reference_genome.fasta reads.(fasta/fastq) K W <br />

In reference_gnome and read files file names can be provided with file path.<br />
Here K and W both are optional. If not provided K = 14 and W = 24
K = length of minimizer to calculate
W = length of window of K-mers

