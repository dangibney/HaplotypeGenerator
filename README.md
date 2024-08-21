```bash
git clone https://github.com/dangibney/HaplotypeGenerator
cd HaplotypeGenerator
make

# test run
./haplotype_generator test_data/MHC-CHM13.0.gfa.gz generated_hap.fa

# to run multiple simulations
chmod +x multiple_simulations.sh
./multiple_simulations.sh
```

## Description
Modify the probability of recombination at the top of src/haplotype_generator.cpp to obtain different number of recombinations in the generated haplotype 
