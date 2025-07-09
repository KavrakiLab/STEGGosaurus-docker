---

# 🦕 STEGGosaurus-Docker

This is a minimal working Docker implementation of **STEGG** for testing purposes. The long-term goal is to make the Docker image available via Docker Hub and host a web server so others can use the tool easily—without ever needing to look at the code (my code is known by the state of California to cause cancer).

---

## Quick Start

### 1. **Download Required Files**

Before building the Docker container, you must download two sets of files:

#### MHCFlurry2.0 Frequency Matrices

Download the frequency matrix CSV from this link:
[MHCFlurry2.0 Frequencies (Box)](https://rice.box.com/s/duxshqxtkykg7u9y3b0tl79j8byc2dp2)

Place the CSV file into the following directory:

```
Ape-Gen2.0-main/helper_files/
```

#### ⚙️ KORP Binary

Download `korp6Dv1.bin` from this link:
[Download KORP](https://chaconlab.org/modeling/korp/down-korp/item/korp-linux)

Copy the binary file to **both** of the following directories:

```
Ape-Gen2.0-main/RCD_required_files/
T-RECS/RCD_required_files/
```

---

## 🐳 Building & Running the Docker Container

Once dependencies are in place:

1. **Build the Docker container**

   run `build_image.sh` and `start_image.sh`

2. **Run the model** inside the container:

   Navegate to:

   ```
   /home/STEGG_controler/
   ```

   And run

   ```bash
   python3 model_complex.py input.json
   ```

2. **Model Additional Complexes**
   To model other TCR–pMHC complexes, create a new JSON input file following the same structure as `input.json`, and run:

   ```bash
   python3 model_complex.py your_input.json
   ```

---
