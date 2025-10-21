
import numpy as np


CUTOFFS = {
    "CD": 20,
    "CE": 20,
    "DE": 10,
    "angle_67": np.pi,
    "angle_111": np.pi,
    "angle_101": (5/6)*np.pi,
    "angle_57": (5/6)*np.pi,
    "angle_90": (7/9)*np.pi,
    "angle_46": (7/9)*np.pi,
    "angle_27": (2/9)*np.pi,
    "angle_31": (2/9)*np.pi,
}


DIRECTIONS = {
    "CD": True,    # <=
    "CE": True,    # <=
    "DE": True,    # <=
    "angle_67": True,    # <=
    "angle_111": True,    # <=
    "angle_101": True,    # <=
    "angle_57": True,    # <=
    "angle_90": True,    # <=
    "angle_46": True,    # <=
    "angle_27": False,    # >=
    "angle_31": False,    # >=
}

APPLIED = {
    "CD": True,     # Minimum C-alpha distance between Peptide (C) and TCR alpha (D)
    "CE": True,     # Minimum C-alpha distance between Peptide (C) and TCR beta (E)
    "DE": True,     # Minimum C-alpha distance between TCR alpha (D) and TCR beta (E)
    # -------------------------------------------------------------------------------------------------------
    # TCR β Normal vs. MHC PCA1 (V11 vs V1)
    "angle_111": True, # pca3E (11) | pca1MHC (1) | TCR β Normal vs. MHC PCA1
    # TCR α Normal vs. MHC PCA1 (V7 vs V1)
    "angle_67": True, # pca3D (7) | pca1MHC (1) | TCR α Normal vs. MHC PCA1
    # -------------------------------------------------------------------------------------------------------
    # TCR β PCA2 vs. MHC PCA2 (V10 vs V2)
    "angle_101": True, # pca2E (10) | pca2MHC (2) | TCR β PCA2 vs. MHC PCA2
    # TCR α PCA2 vs. MHC PCA2 (V6 vs V2)
    "angle_57": True, # pca2D (6) | pca2MHC (2) | TCR α PCA2 vs. MHC PCA2
    # -------------------------------------------------------------------------------------------------------
    # TCR α PCA1 vs. MHC PCA2 (V5 vs V2)
    "angle_46": True, # pca1D (5) | pca2MHC (2) | TCR α PCA1 vs. MHC PCA2
    # TCR β PCA1 vs. MHC PCA2 (V9 vs V2)
    "angle_90": True, # pca1E (9) | pca2MHC (2) | TCR β PCA1 vs. MHC PCA2
    # -------------------------------------------------------------------------------------------------------
    # MHC Normal vs. TCR α PCA1 (V3 vs V5)
    "angle_27": True, # pca3MHC (3) | pca1D (5) | MHC Normal vs. TCR α PCA1
    # MHC Normal vs. TCR β PCA1 (V3 vs V9)
    "angle_31": True, # pca3MHC (3) | pca1E (9) | MHC Normal vs. TCR β PCA1
}