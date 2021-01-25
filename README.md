# Clusterstates

## Paper:
[arxiv.org/abs/2011.08213](https://arxiv.org/abs/2011.08213)


## Getting Started
```
git clone https://github.com/noajshu/clusterstates.git
cd clusterstates/simulator/
make
./simulate
```

## Usage example:
```
./simulate \
    --noise-model bcc_quasistandard \
    --p1 1 --p2  1  --points 2 \
    --loss \
    --pl1 3e-4 --pl2 1e-3 \
    --noise-rate 0.001 \
    --points-loss 7 \
    --min-errors 20 \
    --max-trials 15000 \
    --Ls="5,7,9,11,13,15,17,19,21,23" \
    --fname out.csv
```

## Credits
The most important code was written by Kianna Wan (kcqyw), who contributed the syndrome decoder including supercheck decoding for heralded loss errors, as well as noise models for Protocols A & B of our [paper](https://arxiv.org/abs/2011.08213). We use the (Blossom V)[https://pub.ist.ac.at/~vnk/papers/blossom5.pdf] algorithm by Vladimir Kolmogorov for minimum weight perfect matching.
