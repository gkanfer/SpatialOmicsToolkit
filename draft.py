import STAligner
import os

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg

import scipy
import networkx

import torch

pathout = "/data/kanferg/Sptial_Omics/SpatialOmicsToolkit/out_4"
andata_76 = sc.read_h5ad(os.path.join(pathout, "andata_BreastCancer_76age_STAligner_input.h5ad"))