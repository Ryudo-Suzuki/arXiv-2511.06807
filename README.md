# Overview

This repository contains the datasets and source code used in the paper:

**“Anomalous Enhancement of Yield Strength due to Static Friction”**  
<https://arxiv.org/abs/2511.06807>

**Data and code prepared by:**  
Ryudo Suzuki  
Email: suzuki.ryudo.55e@st.kyoto-u.ac.jp

---
# Repository Structure
```
.
├── figures/
│   ├── fig2a/
│   │   ├── data/
│   │   └── fig/
│   ├── fig2b/
│   │   ├── data/
│   │   └── fig/
│   ├── figS2/
│   │   ├── data/
│   │   └── fig/
│   ├── figS3a/
│   │   ├── data/
│   │   └── fig/
│   └── figS3b/
│       ├── data/
│       └── fig/
└── src/
    ├── compress/
    └── pile/

```
---

# Description of Contents

## figures/

Each directory (`fig2a` – `figS3b`) contains:

- **data/**  
  Raw numerical data used to generate the corresponding figure in the paper.

- **fig/**  
  The plotted figures produced by the Gnuplot script `load.plt`.

A Gnuplot script (`load.plt`) is included in each directory to reproduce the plots.  
The final figures in the paper were created by post-processing these base plots.

---

## src/

### src/compress/

Contains C code used to simulate the compression of three stacked cylinders and compute the resulting yield force.  
This code was used to produce the data for **Fig. 2** and **Fig. S3**.

### src/pile/

Contains:

- **pile.c**  
  Code for generating initial configurations of three stacked cylinders.  
  These configurations are used as input for the compression simulations in `compress/`.

- **animation.c**  
  Code for generating GIF animations that visualize the stacking and deformation process.
