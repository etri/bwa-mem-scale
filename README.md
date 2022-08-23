### BWA-MEM-SCALE

BWA-MEM-SCALE builds upon BWA-MEM2 and BWA-Mich, and includes performance improvements to entire steps of genome sequence alignments.
It adds Exact Match Filter (EMF), FM-index Accelerator (FMA), and various optimization techniques.
BWA-MEM-SCALE gives up to 3.32X speedup compared to BWA-MEM2 when the available memory capacity is over 133GB.
Also, it can reduce the memory consumption with the restricted performance improvement. 
BWA-MEM-SCALE produces identical results as BWA-MEM2 except for MAPQ scores, XS tags, and XA tags in some cases (3.14% in our experiment).

For BWA-MEM2 and BWA-Mich, refer to the following links.
https://github.com/bwa-mem2/bwa-mem2
https://github.com/bwa-mem2/bwa-mem2/tree/ert

## Preparing Alignment
```sh
# Compile from source
git clone --recursive https://github.com/etri/bwa-mem-scale.git bwa-mem-scale
cd bwa-mem-scale

# To find out vectorization features supported in your machine
cat /proc/cpuinfo

# If AVX512BW (512-bit SIMD) is supported
make clean
make -j<num_threads> scale=1 arch=avx512

# If AVX2 (256-bit SIMD) is supported
make clean
make -j<num_threads> scale=1 arch=avx2

# If SSE4.1 (128-bit SIMD) is supported (default)
make -j<num_threads> scale=1

# Build index (Takes ~2 hr for human genome with 56 threads. 1 hr for BWT, 1 hr for ERT)
./bwa-mem2.scale index -p <index prefix> <input.fasta> # Generate FM-index of BWA-MEM2. Take ~1hour.
./bwa-mem2.scale index -a ert -t <num threads> -p <index prefix> <input.fasta) # Generate ERT index. Take 4~5 hours with 8 threads
./bwa-mem2.scale smem-table <index prefix> # Generate FM-index Accelerator (FMA) indices. Take ~1min.
./bwa-mem2.scale perfect-index –l <seed length> <index prefix> # Exact Match Filter (EMF) index. Take ~20min. <seed length> is the minimum read length.

```

### Performing Alignment
```sh
# Enable hugepage (optional)
# To get enough 1GB pages,
echo 121 | sudo tee /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages
# To check the number of available 1GB pages
cat tee /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages
# If the number of available 1GB pages is less than 121, please reboot your system and retry the above.

# Load indices on the in-memory index store (optional)
# If you do not run the following command, the first alignment process will load the indices on the in-memory index store.
# If you want to use hugepage, the following command must be executed before any alignment processes.
# page size: normal (default), 2mb, 1gb
# memory capacity: 17 ~ 121 is the valid range
./bwa-mem2.scale load-shm -H <page size> -l <read length> -g <memory capacity for indices> <index prefix>

# Perform single-end alignment
./bwa-mem2.scale mem -t <num threads> -i <num pipeline> -l <read length> -o <output.sam> <index prefix> <input.fastq>

# Perform pair-end alignment (the performance improvement is restricted.)
./bwa-mem2.scale mem -t <num threads> -i <num pipeline> -l <read length> -o <output.sam> <index prefix> <input_1.fastq> <input_2.fastq>

# Remove the in-memory index store
./bwa-mem2.scale remove-shm
```

## Notes

* The minimum memory requirements of BWA-MEM-SCALE is about 20GB, same with BWA-MEM2. To enable all optimization mechanisms of BWA-MEM-SCALE, 140GB or more memory is required. The total memory consumption mainly depends on the number of threads.

* For the performance results, refer to the paper in the below.

## Citation

If you use BWA-MEM-SCALE, please cite the following [paper](https://doi.org/10.1145/3545008.3545033):

> **Changdae Kim, Kwangwon Koh, Taehoon Kim, Daegyu Han, Jiwon Seo. *BWA-MEM-SCALE: Accelerating Genome Sequence Mapping on Commodity Servers. *ICPP 2022**
